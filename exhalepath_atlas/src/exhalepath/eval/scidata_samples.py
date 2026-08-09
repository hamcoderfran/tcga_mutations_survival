"""Per-sample Sci Data 2024 evaluation (one-vs-rest + stratified AUCs)."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from ..gcms.metrics import compute_diagnostic_metrics, stratified_auroc
from ..gcms.paper_pack import export_paper_pack, write_methods_stub
from ..gcms.preprocess import apply_feature_filter, blank_ratio_filter
from ..gcms.report import write_diagnostic_report
from ..gcms.scidata_samples import (
    COHORT_FILES,
    list_scidata_cohorts,
    load_scidata_intensity_matrix,
    load_scidata_ovr_matrix,
)
from ..gcms.score import disease_signature, score_patients


def _nested_logistic(
    X: np.ndarray, y: np.ndarray, *, n_splits: int = 5, seed: int = 42
) -> dict[str, Any]:
    if len(np.unique(y)) < 2 or len(y) < n_splits * 2:
        return {
            "mean_test_auroc": None,
            "folds": [],
            "note": "insufficient class balance for nested CV",
        }
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    fold_rows = []
    oof = np.full(len(y), np.nan)
    for fold_id, (tr, te) in enumerate(skf.split(X, y)):
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=800, class_weight="balanced", solver="lbfgs"
                    ),
                ),
            ]
        )
        pipe.fit(X[tr], y[tr])
        proba = pipe.predict_proba(X[te])[:, 1]
        oof[te] = proba
        auc = (
            float(roc_auc_score(y[te], proba))
            if len(np.unique(y[te])) >= 2
            else None
        )
        fold_rows.append({"fold_id": fold_id, "n_test": int(len(te)), "auroc": auc})
    aucs = [r["auroc"] for r in fold_rows if r["auroc"] is not None]
    # non-nested optimistic
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(max_iter=800, class_weight="balanced", solver="lbfgs"),
            ),
        ]
    )
    pipe.fit(X, y)
    proba_all = pipe.predict_proba(X)[:, 1]
    non_nested = float(roc_auc_score(y, proba_all))
    nested = float(np.mean(aucs)) if aucs else None
    return {
        "mean_test_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": (non_nested - nested) if nested is not None else None,
        "folds": fold_rows,
        "oof_scores": oof,
    }


def evaluate_scidata_samples(
    *,
    positive_cohort: str = "asthma",
    out_dir: Path | None = None,
    n_splits: int = 5,
    seed: int = 42,
    min_detect_frac: float = 0.3,
    min_blank_ratio: float = 2.0,
    signature_source: str = "hybrid",
    mapped_vocs_only: bool = True,
    make_paper_pack: bool = True,
) -> dict[str, Any]:
    """Run one-vs-rest Sci Data per-sample diagnostic research eval."""
    matrix = load_scidata_ovr_matrix(
        positive_cohort, mapped_vocs_only=mapped_vocs_only
    )
    # blank / detection filter on intensity table
    filt = blank_ratio_filter(
        matrix.matrix,
        min_ratio=min_blank_ratio,
        min_detect_frac=min_detect_frac,
    )
    filtered = apply_feature_filter(matrix.matrix, filt["kept_features"])
    if filtered.shape[1] < 2:
        # fall back to unfiltered if filter too aggressive
        filtered = matrix.matrix
        filt["note"] = (filt.get("note") or "") + " | fallback_unfiltered"

    X = np.log1p(filtered.clip(lower=0).to_numpy(dtype=float))
    y = matrix.labels.loc[filtered.index].to_numpy(dtype=int)
    subject_ids = list(filtered.index.astype(str))

    logistic = _nested_logistic(X, y, n_splits=n_splits, seed=seed)
    # full-fit scores for ROC figure / stratified
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(max_iter=800, class_weight="balanced", solver="lbfgs"),
            ),
        ]
    )
    pipe.fit(X, y)
    scores = pipe.predict_proba(X)[:, 1]
    metrics = compute_diagnostic_metrics(y, scores, seed=seed)

    # covariates for stratification
    age = []
    sex = []
    smoking = []
    sid_to_sample = {s.subject_id: s for s in matrix.samples}
    for sid in subject_ids:
        rec = sid_to_sample.get(sid)
        age.append(rec.age if rec else None)
        sex.append(rec.sex if rec else None)
        smoking.append(rec.smoking_status if rec else None)

    strata_payload: dict[str, Any] = {}
    strata_in = {"age": age, "sex": sex}
    if any(s is not None for s in smoking):
        strata_in["smoking"] = smoking
    strata_payload = stratified_auroc(y, scores, strata_in)

    # optional mechanism signature scores on mapped VOC columns only
    sig_metrics = None
    sig_used: list[str] = []
    try:
        voc_cols = [c for c in filtered.columns if not str(c).startswith("cid:")]
        if len(voc_cols) >= 2:
            from ..gcms.patient_matrix import PatientVOCMatrix

            sub = PatientVOCMatrix(
                study_id=matrix.study_id,
                disease_id=matrix.disease_id,
                disease_name=matrix.disease_name,
                modality=matrix.modality,
                unit=matrix.unit,
                matrix=filtered[voc_cols],
                labels=matrix.labels.loc[filtered.index],
                samples=matrix.samples,
                metadata=matrix.metadata,
            )
            sig = disease_signature(
                matrix.disease_id,
                source=signature_source,  # type: ignore[arg-type]
                top_n=50,
                voc_ids=voc_cols,
            )
            sig_used = sorted(sig.keys())
            if len(sig) >= 2:
                sig_scores = score_patients(
                    sub.log1p(), sig, method="cosine", reference="control_mean"
                )
                y_sig = sub.labels.loc[sig_scores.index].astype(int)
                sig_metrics = compute_diagnostic_metrics(
                    y_sig.tolist(), sig_scores.tolist(), seed=seed
                ).model_dump()
    except Exception as exc:  # noqa: BLE001
        sig_metrics = {"error": str(exc)}

    report: dict[str, Any] = {
        "title": f"Sci Data 2024 per-sample · {positive_cohort} one-vs-rest",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": matrix.study_id,
        "disease_id": matrix.disease_id,
        "disease_name": matrix.disease_name,
        "positive_cohort": positive_cohort,
        "modality": matrix.modality,
        "signature_source": signature_source,
        "score_method": "logistic_ovr",
        "n_voc_features": int(filtered.shape[1]),
        "signature_vocs_used": sig_used,
        "matrix_metadata": matrix.metadata,
        "filter": {k: v for k, v in filt.items() if k != "ratios"},
        "metrics": metrics.model_dump(),
        "signature_metrics": sig_metrics,
        "stratified_auroc": strata_payload,
        "nested": {
            "strategy": "stratified_kfold",
            "seed": seed,
            "content_sha256": None,
            "mean_test_auroc": logistic["mean_test_auroc"],
            "non_nested_auroc": logistic["non_nested_auroc"],
            "optimism_gap": logistic["optimism_gap"],
            "folds": logistic["folds"],
            "logistic_baseline": {
                "mean_test_auroc": logistic["mean_test_auroc"],
                "non_nested_auroc": logistic["non_nested_auroc"],
                "optimism_gap": logistic["optimism_gap"],
                "folds": logistic["folds"],
            },
        },
        "patient_scores": [
            {
                "subject_id": sid,
                "label": int(y[i]),
                "score": float(scores[i]),
                "oof_logistic_score": (
                    float(logistic["oof_scores"][i])
                    if logistic.get("oof_scores") is not None
                    and not np.isnan(logistic["oof_scores"][i])
                    else None
                ),
                "age": age[i],
                "sex": sex[i],
                "smoking": smoking[i],
            }
            for i, sid in enumerate(subject_ids)
        ],
        "cohorts_available": list_scidata_cohorts(),
        "caveats": [
            "One-vs-rest across Asthma/COPD/Bronchiectasis — **no healthy controls**.",
            "Smoking metadata absent in Sci Data 2024 CBD_metadata — smoking strata N/A.",
            "Blank samples absent — blank-ratio filter reduces to detection-fraction.",
            "Bronchiectasis signature uses COPD atlas prior as documented proxy.",
            "Research enablement only — not clinical validation.",
        ],
    }

    pack_info = None
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        write_diagnostic_report(report, out_dir)
        pd.DataFrame(report["patient_scores"]).to_csv(
            out_dir / "patient_scores.csv", index=False
        )
        (out_dir / "filter_report.json").write_text(
            json.dumps({k: v for k, v in filt.items() if k != "ratios"}, indent=2)
        )
        (out_dir / "stratified_auroc.json").write_text(
            json.dumps(strata_payload, indent=2)
        )
        write_methods_stub(
            out_dir / "METHODS.md",
            study_id=matrix.study_id,
            disease_id=matrix.disease_id,
            n_subjects=metrics.n_subjects,
            split_sha256=None,
            signature_source=signature_source,
            extra_notes=report["caveats"],
        )
        md = [
            f"# Sci Data per-sample eval — {positive_cohort} OVR",
            "",
            f"- n={metrics.n_subjects} (pos={metrics.n_positive}, neg={metrics.n_negative})",
            f"- Logistic AUROC (full)={metrics.auroc}",
            f"- Nested logistic AUROC={logistic['mean_test_auroc']}",
            f"- Features kept after filter={filtered.shape[1]} ({filt.get('note')})",
            f"- Stratified AUROC keys: {list(strata_payload.get('strata', {}))}",
            "",
            "> One-vs-rest pulmonary cohorts; not disease-vs-healthy.",
            "",
        ]
        (out_dir / "SCIDATA_SAMPLE_EVAL.md").write_text("\n".join(md))
        (out_dir / "SCIDATA_SAMPLE_EVAL.json").write_text(
            json.dumps(report, indent=2, default=str)
        )
        if make_paper_pack:
            pack_info = export_paper_pack(
                out_dir,
                study_id=matrix.study_id,
                disease_id=matrix.disease_id,
                n_subjects=metrics.n_subjects,
                signature_source=signature_source,
                methods_notes=report["caveats"],
            )
            report["paper_pack"] = pack_info

    return report


def evaluate_all_scidata_cohorts(
    *,
    out_dir: Path = Path("runs/scidata_samples"),
    **kwargs: Any,
) -> dict[str, Any]:
    out_dir = Path(out_dir)
    rows = []
    for cohort in COHORT_FILES:
        try:
            r = evaluate_scidata_samples(
                positive_cohort=cohort, out_dir=out_dir / cohort, **kwargs
            )
            rows.append(
                {
                    "cohort": cohort,
                    "ok": True,
                    "auroc": (r.get("metrics") or {}).get("auroc"),
                    "nested_auroc": (r.get("nested") or {}).get("mean_test_auroc"),
                    "n": (r.get("metrics") or {}).get("n_subjects"),
                    "n_features": r.get("n_voc_features"),
                    "paper_pack": (r.get("paper_pack") or {}).get("zip_path"),
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append({"cohort": cohort, "ok": False, "error": str(exc)})

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohorts": rows,
        "label_scheme": "one_vs_rest_pulmonary_cohorts",
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "SCIDATA_ALL_SUMMARY.json").write_text(json.dumps(summary, indent=2))
    lines = [
        "# Sci Data per-sample summary (one-vs-rest)",
        "",
        "| Cohort | n | Features | Logistic AUROC | Nested AUROC |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        if not r.get("ok"):
            lines.append(f"| {r['cohort']} | — | — | ERROR | {r.get('error')} |")
        else:
            lines.append(
                f"| {r['cohort']} | {r['n']} | {r['n_features']} | "
                f"{r.get('auroc')} | {r.get('nested_auroc')} |"
            )
    lines += [
        "",
        "Negatives = other pulmonary cohorts (no healthy arm).",
        "",
    ]
    (out_dir / "SCIDATA_ALL_SUMMARY.md").write_text("\n".join(lines))
    return summary


def intensity_filter_demo(
    *,
    min_detect_frac: float = 0.5,
    min_blank_ratio: float = 2.0,
) -> dict[str, Any]:
    """Run blank/detection filter on full (mapped+unmapped) Sci Data intensity table."""
    intensity, _, _ = load_scidata_intensity_matrix(mapped_vocs_only=False)
    return blank_ratio_filter(
        intensity, min_ratio=min_blank_ratio, min_detect_frac=min_detect_frac
    )


__all__ = [
    "evaluate_all_scidata_cohorts",
    "evaluate_scidata_samples",
    "intensity_filter_demo",
]
