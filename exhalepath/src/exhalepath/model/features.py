from __future__ import annotations

import numpy as np
import pandas as pd

from ..knowledge.loader import KnowledgeBase
from ..schemas import PathwayScore, TumorContext


FEATURE_PREFIX_PATHWAY = "pw_"


def stage_ordinal(stage: str | None) -> float:
    if not stage:
        return 1.5
    s = stage.upper().replace("STAGE", "").strip()
    for roman, num in [("IV", 4), ("III", 3), ("II", 2), ("I", 1), ("0", 0)]:
        if s.startswith(roman) or roman in s.split():
            return float(num)
    return 1.5


def build_feature_vector(
    *,
    kb: KnowledgeBase,
    disease: dict,
    pathway_scores: list[PathwayScore],
    tumor: TumorContext | None,
    voc_id: str,
    age_years: float | None = None,
    sex: str | None = None,
    smoking_status: str | None = None,
) -> dict[str, float]:
    score_map = {p.pathway_id: p.score for p in pathway_scores}
    feats: dict[str, float] = {}

    for pid in kb.pathways:
        feats[f"{FEATURE_PREFIX_PATHWAY}{pid}"] = float(score_map.get(pid, 0.0))
        coef = float(kb.pathways[pid].get("voc_effects", {}).get(voc_id, 0.0))
        feats[f"emit_{pid}"] = feats[f"{FEATURE_PREFIX_PATHWAY}{pid}"] * coef

    feats["disease_prior_log2fc"] = float(disease.get("voc_log2fc_prior", {}).get(voc_id, 0.0))
    feats["stage_num"] = stage_ordinal(tumor.stage if tumor else None)
    feats["metastatic"] = 1.0 if (tumor and tumor.metastatic) else 0.0
    feats["tumor_burden"] = (
        float(tumor.tumor_burden_proxy) if tumor and tumor.tumor_burden_proxy is not None else 0.3
    )
    feats["age_years"] = float(age_years) if age_years is not None else 60.0
    feats["sex_male"] = 1.0 if sex == "male" else 0.0
    feats["smoke_current"] = 1.0 if smoking_status == "current" else 0.0
    feats["smoke_former"] = 1.0 if smoking_status == "former" else 0.0
    feats["is_cancer"] = 1.0 if disease.get("category") == "cancer" else 0.0
    return feats


def feature_frame_for_training(
    case_features: pd.DataFrame,
    voc_targets: pd.DataFrame,
    kb: KnowledgeBase,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Join case pathway scores with VOC targets into an ML design matrix (vectorized)."""
    score_cols = [c for c in case_features.columns if c.startswith("score_")]
    base = case_features[["case_id", "project_id", "stage_num", "age_years", "gender"] + score_cols].copy()
    base["stage_num"] = base["stage_num"].fillna(1.5)
    base["age_years"] = base["age_years"].fillna(60.0)
    base["sex_male"] = (base["gender"].astype(str).str.lower() == "male").astype(float)

    voc_cols = ["case_id", "project_id", "voc_id", "log2_fold_change", "target_ppb", "healthy_ppb"]
    voc_cols = [c for c in voc_cols if c in voc_targets.columns]
    merged = voc_targets[voc_cols].merge(base, on=["case_id", "project_id"], how="inner")
    pathway_ids = list(kb.pathways.keys())

    # Disease prior lookup table
    prior_rows = []
    for project_id in merged["project_id"].dropna().unique():
        disease = kb.resolve_disease(str(project_id))
        for voc_id in kb.vocs:
            prior_rows.append(
                {
                    "project_id": project_id,
                    "voc_id": voc_id,
                    "disease_prior_log2fc": float(
                        disease.get("voc_log2fc_prior", {}).get(voc_id, 0.0)
                    ),
                }
            )
    priors = pd.DataFrame(prior_rows)
    merged = merged.merge(priors, on=["project_id", "voc_id"], how="left")
    merged["disease_prior_log2fc"] = merged["disease_prior_log2fc"].fillna(0.0)

    feats = pd.DataFrame(index=merged.index)
    feats["disease_prior_log2fc"] = merged["disease_prior_log2fc"].astype(float)
    feats["stage_num"] = merged["stage_num"].astype(float)
    feats["age_years"] = merged["age_years"].astype(float)
    feats["sex_male"] = merged["sex_male"].astype(float)
    feats["metastatic"] = (merged["stage_num"] >= 4).astype(float)
    feats["tumor_burden"] = ((merged["stage_num"] - 0.5) / 4.0).clip(0, 1)
    feats["smoke_current"] = 0.0
    feats["smoke_former"] = 0.0
    feats["is_cancer"] = 1.0

    coef_lookup = {
        (pid, voc_id): float(kb.pathways[pid].get("voc_effects", {}).get(voc_id, 0.0))
        for pid in pathway_ids
        for voc_id in kb.vocs
    }

    for pid in pathway_ids:
        score = merged.get(f"score_{pid}", pd.Series(0.0, index=merged.index)).astype(float)
        feats[f"{FEATURE_PREFIX_PATHWAY}{pid}"] = score
        coefs = merged["voc_id"].map(lambda v, pid=pid: coef_lookup.get((pid, v), 0.0)).astype(float)
        feats[f"emit_{pid}"] = score * coefs

    y = merged["log2_fold_change"].astype(float)
    voc_ids = merged["voc_id"].astype(str)
    return feats.fillna(0.0), y.reset_index(drop=True), voc_ids.reset_index(drop=True)


def mechanistic_log2fc(
    *,
    kb: KnowledgeBase,
    disease: dict,
    pathway_scores: list[PathwayScore],
    voc_id: str,
    tumor: TumorContext | None,
) -> tuple[float, list[str]]:
    """Pure pathway prior (no ML) for interpretability / cold-start."""
    log2fc = float(disease.get("voc_log2fc_prior", {}).get(voc_id, 0.0))
    drivers: list[tuple[str, float]] = []
    for ps in pathway_scores:
        coef = float(kb.pathways[ps.pathway_id].get("voc_effects", {}).get(voc_id, 0.0))
        contrib = coef * ps.score
        log2fc += contrib
        if abs(contrib) > 1e-6:
            drivers.append((ps.pathway_id, contrib))

    if tumor:
        log2fc += 0.08 * (stage_ordinal(tumor.stage) - 1.5)
        if tumor.metastatic:
            log2fc += 0.12
        if tumor.tumor_burden_proxy is not None:
            log2fc += 0.25 * (tumor.tumor_burden_proxy - 0.3)

    drivers.sort(key=lambda x: abs(x[1]), reverse=True)
    return float(log2fc), [d[0] for d in drivers[:5]]
