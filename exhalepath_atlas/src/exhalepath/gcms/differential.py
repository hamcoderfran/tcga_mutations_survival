"""Cross-disease differential diagnosis (DDx) from breath VOC vectors.

Field context (2025–2026): Sci Reports baseline on Sci Data GC-MS reports
fit-on-cohort XGBoost macro-AUC ≈0.998 under nested CV — impressive but
non-transferable. ExhalePath's differentiator is *mechanism-signature DDx*:
rank diseases by template match without fitting labels on the cohort.
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler

from .score import disease_signature, score_observed_vector


# Sci Data pulmonary triad used for multiclass demos
SCIDATA_DISEASE_MAP = {
    "asthma": "asthma",
    "copd": "copd",
    "bronchiectasis": "bronchiectasis",
}


def rank_diseases_for_vector(
    observed: dict[str, float],
    disease_ids: list[str],
    *,
    source: str = "hybrid",
    top_n: int = 40,
    method: str = "cosine",
) -> list[dict[str, Any]]:
    """Rank atlas diseases by signature match to one observed VOC vector."""
    rows = []
    for did in disease_ids:
        try:
            sig = disease_signature(did, source=source, top_n=top_n)  # type: ignore[arg-type]
        except Exception as exc:  # noqa: BLE001
            rows.append({"disease_id": did, "score": None, "error": str(exc)})
            continue
        out = score_observed_vector(observed, sig, method=method)  # type: ignore[arg-type]
        rows.append(
            {
                "disease_id": did,
                "score": out.get("score"),
                "n_overlap": out.get("n_overlap"),
                "method": method,
                "source": source,
            }
        )
    rows.sort(key=lambda r: (-1e9 if r.get("score") is None else -float(r["score"])))
    for i, r in enumerate(rows):
        r["rank"] = i + 1
    return rows


def multiclass_nested_auroc(
    X: pd.DataFrame,
    y: pd.Series,
    *,
    n_splits: int = 5,
    seed: int = 42,
) -> dict[str, Any]:
    """Fit-on-cohort nested multinomial logistic — ceiling, not transferable."""
    le = LabelEncoder()
    y_enc = le.fit_transform(y.astype(str))
    labels = list(le.classes_)
    if len(labels) < 2:
        return {"status": "skipped", "reason": "need_≥2_classes"}
    skf = StratifiedKFold(n_splits=min(n_splits, min(np.bincount(y_enc))), shuffle=True, random_state=seed)
    oof = np.zeros((len(y_enc), len(labels)), dtype=float)
    fold_macro = []
    for tr, te in skf.split(X, y_enc):
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=2000,
                        class_weight="balanced",
                        random_state=seed,
                    ),
                ),
            ]
        )
        pipe.fit(X.iloc[tr].to_numpy(), y_enc[tr])
        proba = pipe.predict_proba(X.iloc[te].to_numpy())
        # align columns to le.classes_
        full = np.zeros((len(te), len(labels)))
        for j, cls_idx in enumerate(pipe.named_steps["clf"].classes_):
            full[:, int(cls_idx)] = proba[:, j]
        oof[te] = full
        try:
            fold_macro.append(
                float(roc_auc_score(y_enc[te], full, multi_class="ovr", average="macro"))
            )
        except Exception:  # noqa: BLE001
            pass
    try:
        macro = float(roc_auc_score(y_enc, oof, multi_class="ovr", average="macro"))
    except Exception:  # noqa: BLE001
        macro = None
    per_class = {}
    for j, lab in enumerate(labels):
        y_bin = (y_enc == j).astype(int)
        if len(np.unique(y_bin)) < 2:
            per_class[lab] = None
            continue
        try:
            per_class[lab] = float(roc_auc_score(y_bin, oof[:, j]))
        except Exception:  # noqa: BLE001
            per_class[lab] = None
    return {
        "status": "ok",
        "kind": "fit_on_cohort_nested_multinomial",
        "classes": labels,
        "macro_ovr_auroc": macro,
        "mean_fold_macro_auroc": float(np.mean(fold_macro)) if fold_macro else None,
        "per_class_ovr_auroc": per_class,
        "n_samples": int(len(y_enc)),
        "n_features": int(X.shape[1]),
        "caveat": (
            "Fit-on-cohort ceiling (labels used in training). Compare to mechanism "
            "signature DDx below — that is the transferable scientific claim."
        ),
    }


def mechanism_signature_ddx(
    X: pd.DataFrame,
    y: pd.Series,
    *,
    disease_map: dict[str, str] | None = None,
    source: str = "hybrid",
    top_n: int = 40,
) -> dict[str, Any]:
    """One-vs-rest AUROC of mechanism signatures without fitting cohort labels.

    For each true class, score patients with that disease's signature (relative
    to cohort mean) and compute OVR AUROC.
    """
    disease_map = disease_map or SCIDATA_DISEASE_MAP
    # relative to cohort mean (no class labels in reference)
    Xp = X.clip(lower=1e-12)
    ref = Xp.mean(axis=0).clip(lower=1e-12)
    rel = np.log2(Xp.div(ref, axis=1))

    classes = sorted(y.astype(str).unique())
    per_class: dict[str, Any] = {}
    scores_for_rank: dict[str, pd.Series] = {}
    for lab in classes:
        did = disease_map.get(lab, lab)
        try:
            sig = disease_signature(did, source=source, top_n=top_n)  # type: ignore[arg-type]
        except Exception as exc:  # noqa: BLE001
            per_class[lab] = {"auroc": None, "error": str(exc), "disease_id": did}
            continue
        sc = []
        for sid in rel.index:
            vec = rel.loc[sid].to_dict()
            out = score_observed_vector(vec, sig, method="cosine")
            sc.append(float(out["score"]) if out.get("score") is not None else 0.0)
        series = pd.Series(sc, index=rel.index)
        scores_for_rank[lab] = series
        y_bin = (y.astype(str) == lab).astype(int)
        auc = None
        if len(np.unique(y_bin)) >= 2:
            try:
                auc = float(roc_auc_score(y_bin, series))
            except Exception:  # noqa: BLE001
                auc = None
        per_class[lab] = {
            "disease_id": did,
            "auroc": auc,
            "n_overlap_median": int(
                np.median(
                    [
                        score_observed_vector(rel.loc[sid].to_dict(), sig).get("n_overlap") or 0
                        for sid in rel.index[:5]
                    ]
                )
            ),
            "source": source,
        }

    # Top-1 accuracy: predicted class = argmax signature score
    if scores_for_rank:
        score_mat = pd.DataFrame(scores_for_rank)
        pred = score_mat.idxmax(axis=1)
        top1 = float((pred.astype(str) == y.astype(str).loc[pred.index]).mean())
    else:
        top1 = None

    aucs = [v["auroc"] for v in per_class.values() if v.get("auroc") is not None]
    return {
        "status": "ok",
        "kind": "mechanism_signature_ddx",
        "source": source,
        "per_class_ovr_auroc": per_class,
        "macro_ovr_auroc": float(np.mean(aucs)) if aucs else None,
        "top1_accuracy": top1,
        "n_samples": int(len(y)),
        "caveat": (
            "Transferable: signatures from ExhalePath hybrid/stack, not fit on these labels. "
            "Expect lower AUROC than fit-on-cohort XGBoost/logistic — that gap is the point."
        ),
    }


def evaluate_scidata_differential(
    *,
    mapped_vocs_only: bool = True,
    seed: int = 42,
    residualize_age_sex: bool = True,
) -> dict[str, Any]:
    """Full Sci Data multiclass pack: fit ceiling vs mechanism DDx ± residualization."""
    from .residualize import residualize_matrix
    from .scidata_samples import load_scidata_intensity_matrix

    X, cov, y = load_scidata_intensity_matrix(mapped_vocs_only=mapped_vocs_only)
    cov = cov.set_index("sample_id")
    cov.index = cov.index.astype(str)
    X.index = X.index.astype(str)
    y.index = y.index.astype(str)
    # log1p intensities
    X_log = np.log1p(X.clip(lower=0))

    ceiling = multiclass_nested_auroc(X_log, y, seed=seed)
    mech = mechanism_signature_ddx(X_log, y, source="hybrid")
    try:
        mech_stack = mechanism_signature_ddx(X_log, y, source="stack")
    except Exception as exc:  # noqa: BLE001
        mech_stack = {"status": "error", "error": str(exc)}

    resid_block: Optional[dict[str, Any]] = None
    if residualize_age_sex:
        try:
            r = residualize_matrix(
                X_log, cov, covariate_cols=["age", "sex"]
            )
            X_res = r["residual_matrix"]
            y_res = y.loc[X_res.index]
            resid_block = {
                "mean_covariate_r2": r["mean_covariate_r2"],
                "covariate_names": r["covariate_names"],
                "fit_on_cohort_ceiling": multiclass_nested_auroc(X_res, y_res, seed=seed),
                "mechanism_hybrid": mechanism_signature_ddx(X_res, y_res, source="hybrid"),
            }
        except Exception as exc:  # noqa: BLE001
            resid_block = {"error": str(exc)}

    return {
        "title": "Sci Data cross-disease differential (mechanism DDx vs fit-on-cohort)",
        "dataset": "Sci Data 2024 clinical breathomics (Kuo et al.)",
        "n_samples": int(len(y)),
        "n_features": int(X_log.shape[1]),
        "class_counts": y.value_counts().to_dict(),
        "field_baseline_note": (
            "Sci Reports 2025 (doi:10.1038/s41598-025-28143-x) reports nested XGBoost "
            "macro-AUC ≈0.998 on this dataset — a fit-on-cohort ceiling, not a "
            "mechanism-transferable claim. Our multinomial logistic nested AUROC is the "
            "open comparable ceiling; hybrid/stack DDx is the transferable score."
        ),
        "fit_on_cohort_ceiling": ceiling,
        "mechanism_ddx_hybrid": mech,
        "mechanism_ddx_stack": mech_stack,
        "age_sex_residualized": resid_block,
    }


__all__ = [
    "evaluate_scidata_differential",
    "mechanism_signature_ddx",
    "multiclass_nested_auroc",
    "rank_diseases_for_vector",
]
