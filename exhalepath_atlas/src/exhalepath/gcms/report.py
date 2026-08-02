"""Paper-ready research packs for GC-MS diagnostic evaluations."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def write_diagnostic_report(
    report: dict[str, Any],
    out_dir: Path,
) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    jp = out_dir / "PATIENT_DIAGNOSTIC_REPORT.json"
    jp.write_text(json.dumps(report, indent=2, default=str))
    paths["json"] = jp

    md = _markdown(report)
    mp = out_dir / "PATIENT_DIAGNOSTIC_REPORT.md"
    mp.write_text(md)
    paths["markdown"] = mp

    # ROC figure if scores present
    scores = report.get("patient_scores") or []
    if scores:
        y = np.asarray([r["label"] for r in scores], dtype=int)
        s = np.asarray([r["score"] for r in scores], dtype=float)
        fig_path = out_dir / "roc_curve.png"
        _plot_roc(y, s, fig_path, title=report.get("title") or "Patient diagnostic ROC")
        paths["roc"] = fig_path

    checklist = out_dir / "TRIPOD_AI_BREATHVOC_CHECKLIST.md"
    checklist.write_text(_tripod_checklist(report))
    paths["checklist"] = checklist

    return paths


def _plot_roc(y: np.ndarray, scores: np.ndarray, path: Path, *, title: str) -> None:
    from sklearn.metrics import RocCurveDisplay, roc_auc_score

    fig, ax = plt.subplots(figsize=(6, 5))
    if len(np.unique(y)) >= 2:
        RocCurveDisplay.from_predictions(y, scores, ax=ax, name="signature match")
        auc = roc_auc_score(y, scores)
        ax.set_title(f"{title}\nAUROC={auc:.3f}")
    else:
        ax.text(0.5, 0.5, "AUROC undefined (single class)", ha="center")
    ax.plot([0, 1], [0, 1], "--", color="#888", lw=1)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _markdown(report: dict[str, Any]) -> str:
    m = report.get("metrics") or {}
    nested = report.get("nested") or {}
    lines = [
        "# Patient-level GC-MS diagnostic research report",
        "",
        f"**Generated:** {report.get('generated_utc')}",
        f"**Study:** `{report.get('study_id')}` · disease `{report.get('disease_id')}`",
        f"**Signature source:** {report.get('signature_source')} · method `{report.get('score_method')}`",
        "",
        "> Research enablement report — not a clinical validation claim.",
        "",
        "## Primary metrics (full-cohort signature match)",
        "",
        f"- n subjects: **{m.get('n_subjects')}** (pos={m.get('n_positive')}, neg={m.get('n_negative')})",
        f"- AUROC: **{_pct(m.get('auroc'))}**  CI95={m.get('auroc_ci95')}",
        f"- AUPRC: **{_pct(m.get('auprc'))}**",
        f"- Sensitivity / Specificity (Youden): **{_pct(m.get('sensitivity'))}** / **{_pct(m.get('specificity'))}**",
        f"- PPV / NPV: **{_pct(m.get('ppv'))}** / **{_pct(m.get('npv'))}**",
        f"- Confusion: `{m.get('confusion')}`",
        f"- Brier (rank-scaled): {m.get('brier')}",
        "",
        "## Nested / locked-split evaluation",
        "",
        f"- Strategy: `{nested.get('strategy')}` · seed={nested.get('seed')} · sha256=`{(nested.get('content_sha256') or '')[:16]}…`",
        f"- Mean test AUROC across folds: **{_pct(nested.get('mean_test_auroc'))}**",
        f"- Non-nested (full-fit) AUROC: **{_pct(nested.get('non_nested_auroc'))}**",
        f"- Optimism gap (non-nested − nested): **{_pct(nested.get('optimism_gap'))}**",
        "",
        "### Per-fold",
        "",
        "| Fold | n_test | AUROC |",
        "|---|---:|---:|",
    ]
    for row in nested.get("folds") or []:
        lines.append(
            f"| {row.get('fold_id')} | {row.get('n_test')} | {_pct(row.get('auroc'))} |"
        )
    lines += [
        "",
        "## Why this matters for science",
        "",
        "- Patient-level holdout (not peak-level) — the unit of clinical inference",
        "- Locked split hash supports preregistration-style reporting",
        "- Mechanism signature (ExhalePath/stack) tested as a transferable template",
        "- Optimism gap surfaces leakage / overfit risk for papers",
        "",
        "## Caveats",
        "",
    ]
    for c in report.get("caveats") or []:
        lines.append(f"- {c}")
    lines.append("")
    return "\n".join(lines)


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


def _tripod_checklist(report: dict[str, Any]) -> str:
    nested = report.get("nested") or {}
    return f"""# TRIPOD+AI / BreathVOC-1.0 lightweight checklist

Auto-filled fields from `voc eval-patient-diagnostic`. Complete remaining items before publication.

| Item | Status / value |
|---|---|
| Title identifies study aim | Patient-level GC-MS VOC diagnostic research enablement |
| Study / data source | {report.get('study_id')} ({report.get('disease_id')}) |
| Eligibility / labels | Disease vs control from MW factors |
| Predictors | Mapped atlas VOC intensities; signature={report.get('signature_source')} |
| Outcome | Binary infection/disease status |
| Sample size | n={ (report.get('metrics') or {}).get('n_subjects') } |
| Missing data handling | Column median fill for sparse mapped VOCs (documented) |
| Internal validation | Locked {nested.get('strategy')} · sha256={(nested.get('content_sha256') or '')[:24]}… |
| Nested vs non-nested | nested={nested.get('mean_test_auroc')} non-nested={nested.get('non_nested_auroc')} gap={nested.get('optimism_gap')} |
| Discrimination | AUROC + AUPRC reported |
| Calibration | Brier on rank-scaled scores (proxy) |
| Fairness / confounders | Smoking/age often unavailable in MW ST000883 — flagged as gap |
| Open science | Code path: `exhalepath.gcms` / `voc eval-patient-diagnostic` |
| Clinical use claim | **None** — research hypothesis / enablement only |

Generated: {datetime.now(timezone.utc).isoformat()}
"""
