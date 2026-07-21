from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import KNOWLEDGE_DIR
from ..knowledge.loader import clear_knowledge_cache
from ..model.predict import ExhalePathPredictor
from ..schemas import DiseaseQuery, TumorContext


def _load_spec(path: Path | None = None) -> dict[str, Any]:
    p = Path(path or KNOWLEDGE_DIR / "multisite_literature_eval.json")
    return json.loads(p.read_text())


def run_multisite_literature_eval(
    *,
    out_dir: Path | None = None,
    mode: str = "hybrid",
) -> dict[str, Any]:
    """
    Evaluate 20 diseases × 5 sites against literature directional VOC expectations.

    Metrics
    -------
    - directional_accuracy: fraction of expect_elevated / expect_suppressed hits
    - min_fold_accuracy: fraction of min_fold thresholds met (primary site)
    - site_sensitivity: fraction of site-boost checks where primary site ≥ alternate mean
    - disease_pass_rate: diseases with directional_accuracy ≥ 0.7 on primary site
    """
    clear_knowledge_cache()
    spec = _load_spec()
    predictor = ExhalePathPredictor(use_opentargets=False)
    out_dir = Path(out_dir or Path("runs/multisite_eval"))
    out_dir.mkdir(parents=True, exist_ok=True)

    case_rows: list[dict[str, Any]] = []
    voc_rows: list[dict[str, Any]] = []
    disease_summaries: list[dict[str, Any]] = []

    for dspec in spec["diseases"]:
        disease = dspec["disease"]
        genes = list(dspec.get("genes") or [])
        sites = list(dspec["sites"])
        primary = sites[0]
        expect_elevated = list(dspec.get("expect_elevated") or [])
        expect_suppressed = list(dspec.get("expect_suppressed") or [])
        min_fold = dict(dspec.get("min_fold") or {})
        site_boost = dict(dspec.get("expect_site_boost") or {})

        site_pred: dict[str, pd.DataFrame] = {}
        for site in sites:
            tumor = TumorContext(primary_site=site, stage="III" if genes else None)
            q = DiseaseQuery(
                disease=disease,
                tumor=tumor,
                mutated_genes=genes,
                mode=mode,  # type: ignore[arg-type]
            )
            res = predictor.predict(q)
            df = res.to_dataframe().set_index("voc_id")
            site_pred[site] = df
            for voc_id, row in df.iterrows():
                voc_rows.append(
                    {
                        "disease": disease,
                        "resolved_disease_id": res.bundle.disease_id,
                        "site": site,
                        "is_primary_site": site == primary,
                        "voc_id": voc_id,
                        "fold_change": float(row["fold_change"]),
                        "predicted_ppb": float(row["predicted_ppb"]),
                        "healthy_ppb": float(row["healthy_ppb"]),
                        "confidence": float(row["confidence"]),
                    }
                )

        pdf = site_pred[primary]
        elev_hits = elev_tot = 0
        for voc_id in expect_elevated:
            elev_tot += 1
            ok = voc_id in pdf.index and float(pdf.loc[voc_id, "fold_change"]) > 1.0
            elev_hits += int(ok)
            case_rows.append(
                {
                    "disease": disease,
                    "site": primary,
                    "check": "elevated",
                    "voc_id": voc_id,
                    "fold_change": float(pdf.loc[voc_id, "fold_change"]) if voc_id in pdf.index else None,
                    "pass": ok,
                }
            )

        supp_hits = supp_tot = 0
        for voc_id in expect_suppressed:
            supp_tot += 1
            ok = voc_id in pdf.index and float(pdf.loc[voc_id, "fold_change"]) < 1.0
            supp_hits += int(ok)
            case_rows.append(
                {
                    "disease": disease,
                    "site": primary,
                    "check": "suppressed",
                    "voc_id": voc_id,
                    "fold_change": float(pdf.loc[voc_id, "fold_change"]) if voc_id in pdf.index else None,
                    "pass": ok,
                }
            )

        fold_hits = fold_tot = 0
        for voc_id, thr in min_fold.items():
            fold_tot += 1
            ok = voc_id in pdf.index and float(pdf.loc[voc_id, "fold_change"]) >= float(thr)
            fold_hits += int(ok)
            case_rows.append(
                {
                    "disease": disease,
                    "site": primary,
                    "check": "min_fold",
                    "voc_id": voc_id,
                    "fold_change": float(pdf.loc[voc_id, "fold_change"]) if voc_id in pdf.index else None,
                    "threshold": float(thr),
                    "pass": ok,
                }
            )

        # Site sensitivity: for boosted VOCs, primary site fold ≥ mean of other sites
        boost_hits = boost_tot = 0
        for site_key, vocs in site_boost.items():
            if site_key not in site_pred:
                # use primary if alias missing
                site_key = primary
            for voc_id in vocs:
                boost_tot += 1
                primary_fold = float(site_pred[primary].loc[voc_id, "fold_change"])
                others = [
                    float(site_pred[s].loc[voc_id, "fold_change"])
                    for s in sites
                    if s != primary and voc_id in site_pred[s].index
                ]
                other_mean = float(np.mean(others)) if others else primary_fold
                # Pass if primary is at least as high as the mean of alternate sites (tolerance 2%)
                ok = primary_fold >= other_mean * 0.98
                boost_hits += int(ok)
                case_rows.append(
                    {
                        "disease": disease,
                        "site": primary,
                        "check": "site_boost",
                        "voc_id": voc_id,
                        "fold_change": primary_fold,
                        "other_sites_mean_fold": other_mean,
                        "pass": ok,
                    }
                )

        dir_tot = elev_tot + supp_tot
        dir_hits = elev_hits + supp_hits
        dir_acc = dir_hits / dir_tot if dir_tot else 1.0
        fold_acc = fold_hits / fold_tot if fold_tot else None
        site_acc = boost_hits / boost_tot if boost_tot else None

        disease_summaries.append(
            {
                "disease": disease,
                "n_sites": len(sites),
                "primary_site": primary,
                "directional_accuracy": dir_acc,
                "directional_hits": dir_hits,
                "directional_total": dir_tot,
                "min_fold_accuracy": fold_acc,
                "min_fold_hits": fold_hits,
                "min_fold_total": fold_tot,
                "site_sensitivity": site_acc,
                "site_boost_hits": boost_hits,
                "site_boost_total": boost_tot,
                "passed": dir_acc >= 0.7 and (fold_acc is None or fold_acc >= 0.5),
                "refs": dspec.get("refs") or [],
            }
        )

    disease_df = pd.DataFrame(disease_summaries)
    case_df = pd.DataFrame(case_rows)
    voc_df = pd.DataFrame(voc_rows)

    dir_hits = int(case_df[case_df["check"].isin(["elevated", "suppressed"])]["pass"].sum())
    dir_tot = int(case_df["check"].isin(["elevated", "suppressed"]).sum())
    fold_sub = case_df[case_df["check"] == "min_fold"]
    site_sub = case_df[case_df["check"] == "site_boost"]

    overall = {
        "n_diseases": int(len(disease_df)),
        "n_sites_per_disease": 5,
        "n_predictions": int(len(voc_df)),
        "mode": mode,
        "directional_accuracy": float(dir_hits / dir_tot) if dir_tot else None,
        "min_fold_accuracy": float(fold_sub["pass"].mean()) if len(fold_sub) else None,
        "site_sensitivity": float(site_sub["pass"].mean()) if len(site_sub) else None,
        "disease_pass_rate": float(disease_df["passed"].mean()),
        "n_diseases_passed": int(disease_df["passed"].sum()),
        "mean_disease_directional_accuracy": float(disease_df["directional_accuracy"].mean()),
        "grading": None,
    }

    # Composite accuracy score (weighted)
    comps = []
    weights = []
    if overall["directional_accuracy"] is not None:
        comps.append(overall["directional_accuracy"])
        weights.append(0.5)
    if overall["min_fold_accuracy"] is not None:
        comps.append(overall["min_fold_accuracy"])
        weights.append(0.3)
    if overall["site_sensitivity"] is not None:
        comps.append(overall["site_sensitivity"])
        weights.append(0.2)
    w = np.array(weights, dtype=float)
    w = w / w.sum()
    composite = float(np.dot(np.array(comps), w)) if comps else 0.0
    overall["composite_accuracy"] = composite

    if composite >= 0.85:
        grade = "A — strong literature agreement for directional VOC modeling"
    elif composite >= 0.75:
        grade = "B — good agreement; site/magnitude nuances imperfect"
    elif composite >= 0.65:
        grade = "C — moderate agreement; useful for hypothesis generation"
    elif composite >= 0.5:
        grade = "D — weak agreement; major gaps vs literature"
    else:
        grade = "F — poor agreement with literature expectations"
    overall["grading"] = grade

    # Persist
    disease_df.to_csv(out_dir / "disease_summary.csv", index=False)
    case_df.to_csv(out_dir / "check_details.csv", index=False)
    voc_df.to_csv(out_dir / "all_predictions.csv", index=False)
    report = {
        "overall": overall,
        "diseases": disease_summaries,
        "failures": case_df[~case_df["pass"]].to_dict(orient="records"),
    }
    (out_dir / "multisite_accuracy_report.json").write_text(json.dumps(report, indent=2))
    return report
