"""Malaria breath diagnostic upgrade pack (ST000883 + literature remaps)."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..gcms.external_malaria import write_external_catalog
from ..gcms.locked_split import lock_split
from ..gcms.malaria_remap import MALARIA_LIT_SIGNATURE
from ..gcms.metrics import (
    bootstrap_metric_ci,
    compute_diagnostic_metrics,
    fixed_sensitivity_metrics,
)
from ..gcms.nested_sparse import learning_curve_nested_auroc, nested_sparse_logistic
from ..gcms.patient_matrix import export_patient_matrix_csv, load_mw_patient_matrix
from ..gcms.paper_pack import export_paper_pack, write_methods_stub
from ..gcms.report import write_diagnostic_report
from ..gcms.score import disease_signature, score_patients


def _hybrid_plus_lit_signature(voc_ids: list[str]) -> dict[str, float]:
    """Transferable signature: hybrid where available, else malaria lit directions."""
    try:
        sig = disease_signature("malaria", source="hybrid", top_n=50, voc_ids=voc_ids)
    except Exception:  # noqa: BLE001
        sig = {}
    for k, v in MALARIA_LIT_SIGNATURE.items():
        if k in voc_ids:
            sig.setdefault(k, float(v))
    return sig


def evaluate_malaria_diagnostic(
    *,
    out_dir: Path | None = None,
    feature_map: str = "malaria_lit",
    n_splits: int = 5,
    seed: int = 42,
    max_features: int = 8,
    external_matrix: Path | None = None,
    external_labels: Path | None = None,
    loso: bool = False,
    pool_external: bool = False,
) -> dict[str, Any]:
    """Full malaria upgrade: remap + nested sparse + fixed-sens + learning curve.

    When ``external_matrix`` (+ labels) is provided and ``loso`` is True, run
    leave-one-study-out vs ST000883 instead of (or in addition to) pooling.
    CSIRO CHMI labels are always exported for readiness; intensity still needs
    a user peak table.
    """
    matrix = load_mw_patient_matrix("ST000883", feature_map=feature_map)
    matrix_l = matrix.log1p()
    primary = matrix  # keep unpooled for LOSO

    # Optional external CSV
    pooled_note = None
    ext = None
    if external_matrix and external_labels:
        from ..gcms.external_malaria import combine_cohorts, load_external_patient_csv

        ext = load_external_patient_csv(
            Path(external_matrix), Path(external_labels), study_id="EXTERNAL_MALARIA"
        )
        if pool_external or not loso:
            matrix = combine_cohorts([primary, ext], study_id="MALARIA_POOLED")
            matrix_l = matrix.log1p()
            pooled_note = (
                f"Pooled ST000883 + EXTERNAL ({ext.metadata.get('n_subjects')} subjects)"
            )

    manifest = lock_split(
        matrix,
        strategy="stratified_kfold",
        n_splits=n_splits,
        seed=seed,
        out_path=(Path(out_dir) / f"{matrix.study_id}_split_manifest.json")
        if out_dir
        else None,
    )

    # --- Transferable: hybrid+lit signature ---
    sig = _hybrid_plus_lit_signature(list(matrix.matrix.columns.astype(str)))
    full_scores = score_patients(
        matrix_l, sig, method="cosine", reference="control_mean"
    )
    y = matrix.labels.loc[full_scores.index].astype(int)
    metrics_sig = compute_diagnostic_metrics(y.tolist(), full_scores.tolist(), seed=seed)
    # nested signature with train-only control means
    from .patient_diagnostic import _fold_control_relative_scores
    from sklearn.metrics import roc_auc_score

    fold_rows = []
    oof_sig = pd.Series(index=matrix.subject_ids, dtype=float)
    for fold in manifest.folds:
        sc = _fold_control_relative_scores(
            matrix_l, sig, fold.train_subject_ids, fold.test_subject_ids, method="cosine"
        )
        for sid, val in sc.items():
            oof_sig.loc[sid] = float(val)
        y_te = matrix.labels.loc[fold.test_subject_ids].to_numpy(dtype=int)
        auc = (
            float(roc_auc_score(y_te, sc.to_numpy()))
            if len(np.unique(y_te)) >= 2
            else None
        )
        fold_rows.append({"fold_id": fold.fold_id, "n_test": len(y_te), "auroc": auc})
    sig_nested = float(np.mean([r["auroc"] for r in fold_rows if r["auroc"] is not None]))

    # --- Fit-on-cohort: nested sparse ---
    sparse = nested_sparse_logistic(
        matrix_l, manifest, max_features=max_features, selector="kbest", seed=seed
    )
    sparse_l1 = nested_sparse_logistic(
        matrix_l, manifest, max_features=max_features, selector="l1", seed=seed
    )

    # Fixed-sens on transferable scores + sparse OOF
    fixed_sig = fixed_sensitivity_metrics(y.tolist(), full_scores.tolist(), seed=seed)
    oof_sp = sparse["oof_scores"].dropna()
    y_sp = matrix.labels.loc[oof_sp.index].astype(int)
    fixed_sparse = fixed_sensitivity_metrics(
        y_sp.tolist(), oof_sp.tolist(), seed=seed
    )
    auprc_ci = bootstrap_metric_ci(y.tolist(), full_scores.tolist(), metric="auprc", seed=seed)

    # Learning curve (n-limitation)
    curve = learning_curve_nested_auroc(matrix_l, seed=seed, max_features=max_features)
    catalog_path = write_external_catalog(
        Path(out_dir) / "external" if out_dir else None
    )

    # CSIRO CHMI label readiness + optional LOSO
    from ..gcms.csiro_chmi import csiro_readiness, export_csiro_label_bundle
    from ..gcms.external_malaria import leave_one_study_out

    csiro = csiro_readiness()
    loso_report: dict[str, Any]
    if loso and ext is not None:
        loso_report = leave_one_study_out(
            [primary.log1p(), ext.log1p() if hasattr(ext, "log1p") else ext],
            signature=sig,
            max_features=max_features,
            seed=seed,
        )
    elif loso:
        loso_report = {
            "status": "skipped",
            "reason": "external_intensity_matrix_required",
            "csiro": csiro,
            "hint": (
                "Export CSIRO MassHunter peaks indexed by Sample name, then: "
                "voc eval-malaria-diagnostic --external-matrix peaks.csv "
                "--external-labels <bundled csiro_chmi_labels.csv> --loso"
            ),
        }
    else:
        loso_report = {
            "status": "not_requested",
            "csiro": csiro,
            "hint": "Pass --loso with --external-matrix/--external-labels for LOSO AUROC",
        }

    report: dict[str, Any] = {
        "title": "Malaria breath diagnostic upgrade (literature remap + nested sparse)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": matrix.study_id,
        "disease_id": "malaria",
        "disease_name": "Malaria",
        "feature_map": feature_map,
        "n_voc_features": matrix.metadata.get("n_voc_features"),
        "voc_ids": matrix.metadata.get("voc_ids"),
        "library_mapping": {
            "n_library_metabolites": matrix.metadata.get("n_library_metabolites"),
            "n_library_metabolites_mapped": matrix.metadata.get(
                "n_library_metabolites_mapped"
            ),
        },
        "pooled_note": pooled_note,
        "transferable_signature": {
            "source": "hybrid+malaria_lit",
            "n_vocs": len(sig),
            "vocs": sorted(sig.keys()),
            "metrics": metrics_sig.model_dump(),
            "nested_auroc": sig_nested,
            "non_nested_auroc": metrics_sig.auroc,
            "optimism_gap": (metrics_sig.auroc - sig_nested)
            if (metrics_sig.auroc is not None and sig_nested is not None)
            else None,
            "auprc_ci95": auprc_ci,
            "folds": fold_rows,
            "fixed_sensitivity": fixed_sig,
        },
        "fit_on_cohort_sparse_kbest": {
            "mean_test_auroc": sparse["mean_test_auroc"],
            "non_nested_auroc": sparse["non_nested_auroc"],
            "optimism_gap": sparse["optimism_gap"],
            "consensus_features": sparse["consensus_features"],
            "folds": sparse["folds"],
            "fixed_sensitivity_oof": fixed_sparse,
        },
        "fit_on_cohort_sparse_l1": {
            "mean_test_auroc": sparse_l1["mean_test_auroc"],
            "non_nested_auroc": sparse_l1["non_nested_auroc"],
            "optimism_gap": sparse_l1["optimism_gap"],
            "consensus_features": sparse_l1["consensus_features"],
            "folds": sparse_l1["folds"],
        },
        "nested": {
            "strategy": manifest.strategy,
            "seed": manifest.seed,
            "content_sha256": manifest.content_sha256,
            "mean_test_auroc": sig_nested,
            "logistic_baseline": {
                "mean_test_auroc": sparse["mean_test_auroc"],
                "non_nested_auroc": sparse["non_nested_auroc"],
                "optimism_gap": sparse["optimism_gap"],
            },
        },
        "metrics": metrics_sig.model_dump(),
        "learning_curve": curve,
        "external_catalog": str(catalog_path),
        "csiro_chmi": csiro,
        "leave_one_study_out": loso_report,
        "story": {
            "transferable_nested_auroc": sig_nested,
            "fit_on_cohort_sparse_ceiling": sparse["mean_test_auroc"],
            "do_not_chase": (
                "Stop tuning ST000883 AUROC. Report transferable nested AUROC vs "
                "fit-on-cohort sparse ceiling; widen n via external intensity / LOSO."
            ),
        },
        "patient_scores": [
            {
                "subject_id": sid,
                "label": int(matrix.labels.loc[sid]),
                "score": float(full_scores.loc[sid]),
                "oof_signature_score": float(oof_sig.loc[sid])
                if sid in oof_sig.index and not np.isnan(oof_sig.loc[sid])
                else None,
                "oof_sparse_score": float(sparse["oof_scores"].loc[sid])
                if sid in sparse["oof_scores"].index
                and not np.isnan(sparse["oof_scores"].loc[sid])
                else None,
            }
            for sid in matrix.subject_ids
        ],
        "caveats": [
            "ST000883 n≈35 — AUROC CIs remain wide (JBR: curves flatten near n≈50).",
            "Schaber 2018: thioethers largely absent; terpenes (pinene/carene) are the remap unlock.",
            "CSIRO CHMI: baseline/peak labels bundled; intensity still needs MassHunter peak export.",
            "JID 2024 Malawi: no open intensity matrix yet.",
            "Fit-on-cohort sparse AUROC is a ceiling, not a transferable clinical claim.",
            "Do not chase further ST000883 AUROC gains — see STOP_CHASING_ST000883.md.",
            "Research enablement only — not a diagnostic device claim.",
        ],
        "comparison_to_baseline_atlas_map": None,
    }

    # Baseline atlas-only for delta reporting
    try:
        base = load_mw_patient_matrix("ST000883", feature_map="atlas")
        report["comparison_to_baseline_atlas_map"] = {
            "n_voc_features_atlas": base.metadata.get("n_voc_features"),
            "n_voc_features_malaria_lit": matrix.metadata.get("n_voc_features")
            if feature_map == "malaria_lit"
            else None,
            "newly_mapped_vocs": sorted(
                set(matrix.matrix.columns) - set(base.matrix.columns)
            )
            if feature_map == "malaria_lit"
            else [],
        }
    except Exception as exc:  # noqa: BLE001
        report["comparison_to_baseline_atlas_map"] = {"error": str(exc)}

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        export_patient_matrix_csv(matrix, out_dir / "matrices")
        # Adapt report shape for write_diagnostic_report
        diag = {
            **report,
            "signature_source": "hybrid+malaria_lit",
            "score_method": "cosine",
            "modality": matrix.modality,
        }
        write_diagnostic_report(diag, out_dir)
        pd.DataFrame(report["patient_scores"]).to_csv(
            out_dir / "patient_scores.csv", index=False
        )
        (out_dir / "MALARIA_DIAGNOSTIC_UPGRADE.json").write_text(
            json.dumps(report, indent=2, default=str)
        )
        (out_dir / "MALARIA_DIAGNOSTIC_UPGRADE.md").write_text(_md(report))
        (out_dir / "fixed_sensitivity.json").write_text(
            json.dumps(
                {"signature": fixed_sig, "sparse_oof": fixed_sparse}, indent=2, default=str
            )
        )
        (out_dir / "learning_curve.json").write_text(json.dumps(curve, indent=2))
        (out_dir / "leave_one_study_out.json").write_text(
            json.dumps(loso_report, indent=2, default=str)
        )
        export_csiro_label_bundle(out_dir / "csiro_chmi")
        (out_dir / "STOP_CHASING_ST000883.md").write_text(_stop_chasing_md(report))
        write_methods_stub(
            out_dir / "METHODS.md",
            study_id=matrix.study_id,
            disease_id="malaria",
            n_subjects=int(matrix.metadata.get("n_subjects") or 0),
            split_sha256=manifest.content_sha256,
            signature_source="hybrid+malaria_lit",
            extra_notes=report["caveats"],
        )
        pack = export_paper_pack(
            out_dir,
            study_id=matrix.study_id,
            disease_id="malaria",
            n_subjects=int(matrix.metadata.get("n_subjects") or 0),
            split_sha256=manifest.content_sha256,
            signature_source="hybrid+malaria_lit",
            methods_notes=report["caveats"],
        )
        report["paper_pack"] = pack

    return report


def _md(report: dict[str, Any]) -> str:
    t = report.get("transferable_signature") or {}
    sk = report.get("fit_on_cohort_sparse_kbest") or {}
    sl = report.get("fit_on_cohort_sparse_l1") or {}
    cmp_ = report.get("comparison_to_baseline_atlas_map") or {}
    lines = [
        "# Malaria breath diagnostic upgrade",
        "",
        f"Generated: {report.get('generated_utc')}",
        f"Study: `{report.get('study_id')}` · feature_map=`{report.get('feature_map')}`",
        "",
        "> Transferable signature vs fit-on-cohort sparse — research enablement only.",
        "",
        "## Feature remapping",
        "",
        f"- Library metabolites mapped: {report.get('library_mapping')}",
        f"- Atlas-only VOC count: {cmp_.get('n_voc_features_atlas')}",
        f"- Malaria-lit VOC count: {cmp_.get('n_voc_features_malaria_lit')}",
        f"- Newly mapped: {', '.join(cmp_.get('newly_mapped_vocs') or []) or '—'}",
        "",
        "## Transferable (hybrid + malaria literature signature)",
        "",
        f"- Nested AUROC: **{_pct(t.get('nested_auroc'))}**",
        f"- Non-nested AUROC: **{_pct(t.get('non_nested_auroc'))}**",
        f"- Optimism gap: **{_pct(t.get('optimism_gap'))}**",
        f"- AUPRC: **{_pct((t.get('metrics') or {}).get('auprc'))}** CI95={t.get('auprc_ci95')}",
        "",
        "### Fixed-sensitivity operating points (signature)",
        "",
    ]
    for p in (t.get("fixed_sensitivity") or {}).get("points") or []:
        lines.append(
            f"- sens≥{p.get('target_sensitivity')}: "
            f"achieved sens={_pct(p.get('sensitivity'))}, "
            f"spec={_pct(p.get('specificity'))} CI95={p.get('specificity_ci95')}"
        )
    lines += [
        "",
        "## Fit-on-cohort nested sparse (ceiling)",
        "",
        f"- SelectKBest nested AUROC: **{_pct(sk.get('mean_test_auroc'))}** "
        f"(gap={_pct(sk.get('optimism_gap'))})",
        f"- L1 nested AUROC: **{_pct(sl.get('mean_test_auroc'))}** "
        f"(gap={_pct(sl.get('optimism_gap'))})",
        f"- Consensus features (kbest): {', '.join(sk.get('consensus_features') or [])}",
        "",
        "## Learning curve (n-limitation)",
        "",
        "| n | mean nested-ish AUROC | std |",
        "|---:|---:|---:|",
    ]
    for r in (report.get("learning_curve") or {}).get("rows") or []:
        lines.append(
            f"| {r.get('n')} | {_pct(r.get('mean_auroc'))} | {r.get('std_auroc')} |"
        )
    story = report.get("story") or {}
    loso = report.get("leave_one_study_out") or {}
    csiro = report.get("csiro_chmi") or {}
    lines += [
        "",
        "## Transferable vs ceiling (the story — stop chasing ST000883 AUROC)",
        "",
        f"- Transferable nested AUROC: **{_pct(story.get('transferable_nested_auroc'))}**",
        f"- Fit-on-cohort sparse ceiling: **{_pct(story.get('fit_on_cohort_sparse_ceiling'))}**",
        f"- {story.get('do_not_chase')}",
        "",
        "## External cohorts / LOSO",
        "",
        f"- Catalog: `{report.get('external_catalog')}`",
        f"- Pooled: {report.get('pooled_note') or 'ST000883 only (no open second intensity table)'}",
        f"- CSIRO labels bundled: {csiro.get('labels_bundled')} "
        f"(n_subjects={csiro.get('n_subjects')}, "
        f"labeled baseline/peak={csiro.get('n_labeled_baseline_peak')})",
        f"- Intensity status: {csiro.get('intensity_status')}",
        f"- LOSO status: {loso.get('status')} · "
        f"mean signature AUROC={_pct(loso.get('mean_auroc_signature'))} · "
        f"mean sparse AUROC={_pct(loso.get('mean_auroc_sparse'))}",
        "",
        "## Caveats",
        "",
    ]
    for c in report.get("caveats") or []:
        lines.append(f"- {c}")
    lines.append("")
    return "\n".join(lines)


def _stop_chasing_md(report: dict[str, Any]) -> str:
    story = report.get("story") or {}
    t = report.get("transferable_signature") or {}
    sk = report.get("fit_on_cohort_sparse_kbest") or {}
    return "\n".join(
        [
            "# Stop chasing ST000883 AUROC",
            "",
            "ST000883 (n≈35) is a **methods harness**, not a clinical benchmark to hill-climb.",
            "",
            "## Report these two numbers",
            "",
            f"1. **Transferable nested AUROC** (hybrid+lit signature): {_pct(story.get('transferable_nested_auroc') or t.get('nested_auroc'))}",
            f"2. **Fit-on-cohort sparse ceiling** (nested SelectKBest): {_pct(story.get('fit_on_cohort_sparse_ceiling') or sk.get('mean_test_auroc'))}",
            "",
            "The gap is expected: signatures were not trained on these labels; sparse logistic was.",
            "",
            "## What to do instead",
            "",
            "- Add a second **intensity** cohort and run `--loso` (CSIRO labels are bundled; export peaks).",
            "- Widen n (JBR: learning curves flatten near n≈50).",
            "- Report fixed-sensitivity operating points + AUPRC CI95.",
            "- Use ST003200 / Sci Data for confounder (age/sex/smoking) diligence — not malaria AUROC.",
            "",
            "## Do not",
            "",
            "- Retune atlas priors or remaps solely to raise ST000883 nested AUROC.",
            "- Market sparse ceiling (~70%+) as transferable diagnostic performance.",
            "- Claim FDA/clinical readiness from this n.",
            "",
        ]
    )


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


__all__ = ["evaluate_malaria_diagnostic"]
