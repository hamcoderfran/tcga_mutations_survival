"""Patient-level GC-MS diagnostic evaluation harness.

Closes the research gap named in SOTA.md: no ROC/AUC on locked patient holdouts.
Uses bundled Metabolomics Workbench studies (ST000883 malaria, ST000587 HF).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from ..gcms.locked_split import lock_split
from ..gcms.metrics import compute_diagnostic_metrics, stratified_auroc
from ..gcms.patient_matrix import (
    export_patient_matrix_csv,
    list_bundled_diagnostic_studies,
    load_mw_patient_matrix,
)
from ..gcms.paper_pack import export_paper_pack, write_methods_stub
from ..gcms.report import write_diagnostic_report
from ..gcms.score import disease_signature, score_patients


def _fold_control_relative_scores(
    matrix,
    signature: dict[str, float],
    train_ids: list[str],
    test_ids: list[str],
    *,
    method: str,
) -> pd.Series:
    """Score test patients using control mean estimated on train only (leakage-safe)."""
    X = matrix.matrix.clip(lower=1e-12)
    train_ctrl = [i for i in train_ids if int(matrix.labels.loc[i]) == 0]
    if len(train_ctrl) < 2:
        ref = X.loc[train_ids].mean(axis=0)
    else:
        ref = X.loc[train_ctrl].mean(axis=0)
    ref = ref.clip(lower=1e-12)
    profiles = np.log2(X.loc[test_ids].div(ref, axis=1))
    common = [c for c in profiles.columns if c in signature and abs(signature[c]) > 1e-9]
    if len(common) < 2:
        return pd.Series(dtype=float)
    P = profiles[common].to_numpy(dtype=float)
    s = np.asarray([signature[c] for c in common], dtype=float)
    out = []
    for i in range(P.shape[0]):
        x = P[i]
        if method == "dot":
            out.append(float(np.dot(x, s)))
        else:
            out.append(float(np.dot(x, s) / (np.linalg.norm(x) * np.linalg.norm(s) + 1e-12)))
    return pd.Series(out, index=test_ids)


def _nested_logistic_auroc(matrix, manifest) -> dict[str, Any]:
    """Data-driven nested logistic baseline (fit on train VOC matrix only)."""
    fold_rows = []
    oof_scores = pd.Series(index=matrix.subject_ids, dtype=float)
    X_all = matrix.matrix.to_numpy(dtype=float)
    y_all = matrix.labels.to_numpy(dtype=int)
    id_to_i = {sid: i for i, sid in enumerate(matrix.subject_ids)}

    for fold in manifest.folds:
        tr = [id_to_i[s] for s in fold.train_subject_ids]
        te = [id_to_i[s] for s in fold.test_subject_ids]
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
        pipe.fit(X_all[tr], y_all[tr])
        proba = pipe.predict_proba(X_all[te])[:, 1]
        for sid, p in zip(fold.test_subject_ids, proba):
            oof_scores.loc[sid] = float(p)
        if len(np.unique(y_all[te])) >= 2 and len(te) >= 2:
            auc = float(roc_auc_score(y_all[te], proba))
        else:
            auc = None
        fold_rows.append(
            {
                "fold_id": fold.fold_id,
                "n_test": len(te),
                "auroc": auc,
                "model": "logistic_voc_matrix",
            }
        )

    aucs = [r["auroc"] for r in fold_rows if r["auroc"] is not None]
    # non-nested: fit on all, score all (optimistic)
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(max_iter=800, class_weight="balanced", solver="lbfgs"),
            ),
        ]
    )
    pipe.fit(X_all, y_all)
    proba_all = pipe.predict_proba(X_all)[:, 1]
    non_nested = float(roc_auc_score(y_all, proba_all)) if len(np.unique(y_all)) >= 2 else None
    nested = float(np.mean(aucs)) if aucs else None
    return {
        "folds": fold_rows,
        "mean_test_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": (non_nested - nested) if (non_nested is not None and nested is not None) else None,
        "oof_scores": oof_scores,
    }


def evaluate_patient_diagnostic(
    *,
    study_id: str = "ST000883",
    signature_source: Literal["hybrid", "stack", "literature"] = "hybrid",
    score_method: Literal["cosine", "dot"] = "cosine",
    split_strategy: str = "stratified_kfold",
    n_splits: int = 5,
    seed: int = 42,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    matrix = load_mw_patient_matrix(study_id)
    matrix_l = matrix.log1p()
    sig = disease_signature(
        matrix.disease_id,
        source=signature_source,
        top_n=50,
        voc_ids=list(matrix.matrix.columns),
    )

    # Full-cohort signature scores (mechanism template — not fit on these labels)
    full_scores = score_patients(
        matrix_l, sig, method=score_method, reference="control_mean"
    )
    y = matrix.labels.loc[full_scores.index].astype(int)
    metrics = compute_diagnostic_metrics(y.tolist(), full_scores.tolist(), seed=seed)

    # Locked splits
    manifest = lock_split(
        matrix,
        strategy=split_strategy,
        n_splits=n_splits,
        seed=seed,
        out_path=(Path(out_dir) / f"{study_id}_split_manifest.json") if out_dir else None,
    )

    # Nested signature match with train-only control means
    fold_rows = []
    oof = pd.Series(index=matrix.subject_ids, dtype=float)
    for fold in manifest.folds:
        sc = _fold_control_relative_scores(
            matrix_l,
            sig,
            fold.train_subject_ids,
            fold.test_subject_ids,
            method=score_method,
        )
        for sid, val in sc.items():
            oof.loc[sid] = float(val)
        y_te = matrix.labels.loc[fold.test_subject_ids].to_numpy(dtype=int)
        if len(np.unique(y_te)) >= 2 and len(y_te) >= 2:
            auc = float(roc_auc_score(y_te, sc.to_numpy()))
        else:
            auc = None
        fold_rows.append(
            {
                "fold_id": fold.fold_id,
                "n_test": len(fold.test_subject_ids),
                "auroc": auc,
                "model": f"signature_{signature_source}_{score_method}",
            }
        )
    sig_aucs = [r["auroc"] for r in fold_rows if r["auroc"] is not None]
    sig_nested = float(np.mean(sig_aucs)) if sig_aucs else None
    sig_non_nested = metrics.auroc
    logistic = _nested_logistic_auroc(matrix_l, manifest)

    # Smoking / age / sex stratified AUCs when sample metadata exist
    sid_to_rec = {s.subject_id: s for s in matrix.samples}
    age_vals = [getattr(sid_to_rec.get(sid), "age", None) for sid in full_scores.index]
    sex_vals = [getattr(sid_to_rec.get(sid), "sex", None) for sid in full_scores.index]
    smoke_vals = [
        getattr(sid_to_rec.get(sid), "smoking_status", None) for sid in full_scores.index
    ]
    strata_in: dict[str, list] = {}
    if any(v is not None for v in age_vals):
        strata_in["age"] = age_vals
    if any(v is not None for v in sex_vals):
        strata_in["sex"] = sex_vals
    if any(v is not None for v in smoke_vals):
        strata_in["smoking"] = smoke_vals
    strata = (
        stratified_auroc(y.tolist(), full_scores.tolist(), strata_in)
        if strata_in
        else {
            "overall": metrics.auroc,
            "strata": {},
            "note": "age/sex/smoking metadata unavailable for this study",
        }
    )

    report: dict[str, Any] = {
        "title": f"{matrix.disease_name} patient-level GC-MS diagnostic ({study_id})",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": study_id,
        "disease_id": matrix.disease_id,
        "disease_name": matrix.disease_name,
        "modality": matrix.modality,
        "signature_source": signature_source,
        "score_method": score_method,
        "n_voc_features": matrix.metadata.get("n_voc_features"),
        "signature_vocs_used": sorted(sig.keys()),
        "matrix_metadata": matrix.metadata,
        "metrics": metrics.model_dump(),
        "stratified_auroc": strata,
        "nested": {
            "strategy": manifest.strategy,
            "seed": manifest.seed,
            "content_sha256": manifest.content_sha256,
            "mean_test_auroc": sig_nested,
            "non_nested_auroc": sig_non_nested,
            "optimism_gap": (sig_non_nested - sig_nested)
            if (sig_non_nested is not None and sig_nested is not None)
            else None,
            "folds": fold_rows,
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
                "label": int(matrix.labels.loc[sid]),
                "score": float(full_scores.loc[sid]),
                "oof_signature_score": float(oof.loc[sid])
                if sid in oof.index and not np.isnan(oof.loc[sid])
                else None,
                "oof_logistic_score": float(logistic["oof_scores"].loc[sid])
                if sid in logistic["oof_scores"].index
                and not np.isnan(logistic["oof_scores"].loc[sid])
                else None,
                "age": getattr(sid_to_rec.get(sid), "age", None),
                "sex": getattr(sid_to_rec.get(sid), "sex", None),
                "smoking": getattr(sid_to_rec.get(sid), "smoking_status", None),
            }
            for sid in matrix.subject_ids
        ],
        "bundled_studies": list_bundled_diagnostic_studies(),
        "caveats": [
            "Mapped atlas VOC subset only — many GC-MS peaks are unmapped and dropped.",
            "Smoking/age stratified AUCs reported only when SampleRecord metadata exist.",
            "Mechanism signature is independent of these patient labels, but VOC name mapping can still introduce circularity with literature priors.",
            "Small n (≈35) → wide bootstrap CIs; treat AUROC as research enablement evidence, not clinical validation.",
            "Not a medical device. No clinical diagnostic claim.",
        ],
        "impact": {
            "research_question": (
                "Can an open mechanism-aware exhaled-VOC signature discriminate "
                "disease vs control on public patient-level GC-MS intensities?"
            ),
            "enables": [
                "Preregisterable locked splits (SHA256)",
                "Paper-ready ROC / sens / spec / confusion + TRIPOD+AI checklist stub",
                "Optimism-gap reporting (nested vs non-nested)",
                "Comparison of mechanism signature vs data-fit logistic baseline",
                "One-zip paper pack (figures + Methods + overlay + split hash)",
            ],
        },
    }

    if out_dir is not None:
        out_dir = Path(out_dir)
        export_patient_matrix_csv(matrix, out_dir / "matrices")
        write_diagnostic_report(report, out_dir)
        # also persist nested oof table
        pd.DataFrame(report["patient_scores"]).to_csv(
            out_dir / "patient_scores.csv", index=False
        )
        (out_dir / "stratified_auroc.json").write_text(
            __import__("json").dumps(strata, indent=2, default=str)
        )
        write_methods_stub(
            out_dir / "METHODS.md",
            study_id=study_id,
            disease_id=matrix.disease_id,
            n_subjects=metrics.n_subjects,
            split_sha256=manifest.content_sha256,
            signature_source=signature_source,
            extra_notes=report["caveats"],
        )
        pack = export_paper_pack(
            out_dir,
            study_id=study_id,
            disease_id=matrix.disease_id,
            n_subjects=metrics.n_subjects,
            split_sha256=manifest.content_sha256,
            signature_source=signature_source,
            methods_notes=report["caveats"],
        )
        report["paper_pack"] = pack

    return report


def run_multi_study_diagnostic(
    *,
    studies: Optional[list[str]] = None,
    signature_source: Literal["hybrid", "stack", "literature"] = "hybrid",
    out_dir: Path = Path("runs/patient_diagnostic"),
) -> dict[str, Any]:
    studies = studies or ["ST000883", "ST000587"]
    out_dir = Path(out_dir)
    summaries = []
    for sid in studies:
        try:
            r = evaluate_patient_diagnostic(
                study_id=sid,
                signature_source=signature_source,
                out_dir=out_dir / sid,
            )
            summaries.append(
                {
                    "study_id": sid,
                    "disease_id": r["disease_id"],
                    "auroc": (r.get("metrics") or {}).get("auroc"),
                    "nested_auroc": (r.get("nested") or {}).get("mean_test_auroc"),
                    "logistic_nested_auroc": (
                        (r.get("nested") or {}).get("logistic_baseline") or {}
                    ).get("mean_test_auroc"),
                    "n": (r.get("metrics") or {}).get("n_subjects"),
                    "ok": True,
                }
            )
        except Exception as exc:  # noqa: BLE001
            summaries.append({"study_id": sid, "ok": False, "error": str(exc)})

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "signature_source": signature_source,
        "studies": summaries,
        "notes": [
            "Multi-study research enablement summary — not a pooled clinical validation.",
        ],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "MULTI_STUDY_SUMMARY.json").write_text(
        __import__("json").dumps(summary, indent=2)
    )
    lines = [
        "# Multi-study patient GC-MS diagnostic summary",
        "",
        f"Signature source: `{signature_source}`",
        "",
        "| Study | Disease | n | Signature AUROC | Nested sig AUROC | Nested logistic AUROC |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for s in summaries:
        if not s.get("ok"):
            lines.append(f"| {s['study_id']} | — | — | ERROR | {s.get('error')} | — |")
            continue
        lines.append(
            f"| {s['study_id']} | {s['disease_id']} | {s['n']} | "
            f"{_fmt(s.get('auroc'))} | {_fmt(s.get('nested_auroc'))} | "
            f"{_fmt(s.get('logistic_nested_auroc'))} |"
        )
    lines.append("")
    (out_dir / "MULTI_STUDY_SUMMARY.md").write_text("\n".join(lines))
    return summary


def run_signature_benchmark(
    *,
    study_id: str = "ST000883",
    out_dir: Path = Path("runs/patient_diagnostic/benchmark"),
) -> dict[str, Any]:
    """Compare hybrid / stack / literature × cosine / dot on one study."""
    out_dir = Path(out_dir)
    rows = []
    for src in ("hybrid", "stack", "literature"):
        for method in ("cosine", "dot"):
            try:
                r = evaluate_patient_diagnostic(
                    study_id=study_id,
                    signature_source=src,  # type: ignore[arg-type]
                    score_method=method,  # type: ignore[arg-type]
                    out_dir=out_dir / f"{src}_{method}",
                )
                rows.append(
                    {
                        "signature": src,
                        "method": method,
                        "auroc": r["metrics"]["auroc"],
                        "nested_auroc": r["nested"]["mean_test_auroc"],
                        "logistic_nested": r["nested"]["logistic_baseline"][
                            "mean_test_auroc"
                        ],
                        "auprc": r["metrics"]["auprc"],
                        "ok": True,
                    }
                )
            except Exception as exc:  # noqa: BLE001
                rows.append(
                    {"signature": src, "method": method, "ok": False, "error": str(exc)}
                )

    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": study_id,
        "rows": rows,
        "takeaway": (
            "Data-fit nested logistic is an upper reference on the same mapped VOC "
            "features. Mechanism signatures are transferable templates; literature "
            "panels can score higher when they encode the same cohort's published "
            "directions (partly circular). Report nested AUROC + optimism gap."
        ),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "SIGNATURE_BENCHMARK.json").write_text(
        __import__("json").dumps(payload, indent=2)
    )
    lines = [
        f"# Signature benchmark — {study_id}",
        "",
        "| Signature | Method | AUROC | Nested AUROC | Nested logistic | AUPRC |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        if not r.get("ok"):
            lines.append(f"| {r['signature']} | {r['method']} | ERROR | — | — | — |")
            continue
        lines.append(
            f"| {r['signature']} | {r['method']} | {_fmt(r['auroc'])} | "
            f"{_fmt(r['nested_auroc'])} | {_fmt(r['logistic_nested'])} | {_fmt(r['auprc'])} |"
        )
    lines += ["", payload["takeaway"], ""]
    (out_dir / "SIGNATURE_BENCHMARK.md").write_text("\n".join(lines))
    return payload


def _fmt(x: Any) -> str:
    if x is None:
        return "—"
    return f"{100 * float(x):.1f}%"
