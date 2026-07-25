"""Probe predictions for novel diseases with no VOC prior background.

These queries must resolve as ``custom::*`` (empty ``voc_log2fc_prior``).
Expected behavior: near-healthy VOC panel, lower confidence, explicit unresolved note.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ..biomarker import ExhaleBiomarkerEngine
from ..knowledge.loader import clear_knowledge_cache

# Deliberately outside the curated atlas (no VOC priors / pathway bias).
NOVEL_DISEASES: list[dict[str, Any]] = [
    {"disease": "completely novel ailment QZZ-42", "location": None},
    {"disease": "Xylophage syndrome type Zeta", "location": "lung"},
    {"disease": "hyperblue mitochondrial spark disease", "location": "brain"},
    {"disease": "quantum itch disorder", "location": "skin"},
    {"disease": "fibroquartz encephalopathy", "location": "brain"},
    {"disease": "neon teal cholangioflux", "location": "liver"},
    {"disease": "sporadic purple glomerulopathy XYZ", "location": "kidney"},
    {"disease": "astral cartilage liquefaction", "location": "bone"},
    {"disease": "cryptic umbra pancreatitis variant 9", "location": "pancreas"},
    {"disease": "nonexistent pathogen Omega-7 breath plague", "location": "lung"},
]

# Atlas diseases WITH VOC background for contrast
KNOWN_CONTRAST: list[dict[str, Any]] = [
    {"disease": "type 2 diabetes", "location": None, "expect_voc": "acetone"},
    {"disease": "asthma", "location": "lung", "expect_voc": "pentane"},
    {"disease": "COPD", "location": "lung", "expect_voc": "hexanal"},
    {"disease": "lung adenocarcinoma", "location": "lung", "expect_voc": "hexanal"},
]


def _summary_from_report(report, *, query: dict[str, Any]) -> dict[str, Any]:
    preds = list(report.result.bundle.predictions)
    log2 = {p.voc_id: float(p.log2_fold_change) for p in preds}
    delta = {p.voc_id: float(p.delta_ppb) for p in preds}
    conf = {p.voc_id: float(getattr(p, "confidence", 0.0) or 0.0) for p in preds}
    abs_log2 = np.array([abs(v) for v in log2.values()], dtype=float)
    top = sorted(preds, key=lambda p: abs(p.delta_ppb), reverse=True)[:8]
    did = report.disease_id
    unresolved = did.startswith("custom::") or bool(
        getattr(report, "notes", None)
        and any("not in the curated" in str(n) for n in (report.notes or []))
    )
    notes = list(getattr(report, "notes", None) or [])
    return {
        "query": query,
        "resolved_disease_id": did,
        "resolved_disease_name": report.disease_name,
        "unresolved": unresolved,
        "n_vocs": len(preds),
        "peak_abs_log2fc": float(abs_log2.max()) if len(abs_log2) else 0.0,
        "mean_abs_log2fc": float(abs_log2.mean()) if len(abs_log2) else 0.0,
        "median_abs_log2fc": float(np.median(abs_log2)) if len(abs_log2) else 0.0,
        "mean_confidence": float(np.mean(list(conf.values()))) if conf else None,
        "top_vocs": [
            {
                "voc_id": p.voc_id,
                "log2fc": float(p.log2_fold_change),
                "delta_ppb": float(p.delta_ppb),
                "predicted_ppb": float(p.predicted_ppb),
                "confidence": float(getattr(p, "confidence", 0.0) or 0.0),
            }
            for p in top
        ],
        "acetone_log2fc": log2.get("acetone"),
        "acetone_delta_ppb": delta.get("acetone"),
        "notes": notes[:6],
        "near_healthy": float(abs_log2.max()) < 0.45 if len(abs_log2) else True,
    }


def run_novel_disease_eval(
    *,
    out_dir: Path | None = None,
    mode: str = "hybrid",
) -> dict[str, Any]:
    out_dir = Path(out_dir or Path("runs/novel_diseases"))
    out_dir.mkdir(parents=True, exist_ok=True)
    clear_knowledge_cache()
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)

    novel_rows = []
    for q in NOVEL_DISEASES:
        report = engine.predict(
            q["disease"],
            location=q.get("location"),
            top_n=50,
            mode=mode,
            explain=False,
        )
        novel_rows.append(_summary_from_report(report, query=q))

    known_rows = []
    for q in KNOWN_CONTRAST:
        report = engine.predict(
            q["disease"],
            location=q.get("location"),
            top_n=50,
            mode=mode,
            explain=False,
        )
        row = _summary_from_report(report, query=q)
        row["expect_voc"] = q.get("expect_voc")
        known_rows.append(row)

    n_novel = len(novel_rows)
    n_unresolved = sum(1 for r in novel_rows if r["unresolved"])
    n_near = sum(1 for r in novel_rows if r["near_healthy"])
    mean_peak = float(np.mean([r["peak_abs_log2fc"] for r in novel_rows]))
    mean_conf = float(
        np.mean([r["mean_confidence"] for r in novel_rows if r["mean_confidence"] is not None])
    )

    # Separation: known T2D acetone vs novel mean acetone
    t2d = next(r for r in known_rows if "diabetes" in r["query"]["disease"])
    novel_acetone = float(np.mean([r["acetone_delta_ppb"] or 0.0 for r in novel_rows]))
    acetone_sep = float((t2d.get("acetone_delta_ppb") or 0.0) - novel_acetone)

    gates = {
        "all_novel_unresolved": n_unresolved == n_novel,
        "most_near_healthy": n_near >= int(0.8 * n_novel),
        "mean_peak_abs_log2fc_lt_0_6": mean_peak < 0.6,
        "acetone_sep_vs_t2d_gt_50ppb": acetone_sep > 50.0,
        "known_not_unresolved": all(not r["unresolved"] for r in known_rows),
    }

    report = {
        "version": "1.0.0",
        "description": (
            "Novel disease probe: queries with no curated VOC priors / pathway bias. "
            "Should stay near-healthy vs atlas diseases with VOC background."
        ),
        "n_novel": n_novel,
        "n_unresolved": n_unresolved,
        "n_near_healthy": n_near,
        "mean_peak_abs_log2fc": mean_peak,
        "mean_confidence": mean_conf,
        "novel_mean_acetone_delta_ppb": novel_acetone,
        "t2d_acetone_delta_ppb": t2d.get("acetone_delta_ppb"),
        "acetone_sep_t2d_minus_novel_ppb": acetone_sep,
        "gates": gates,
        "passed": all(gates.values()),
        "novel_cases": novel_rows,
        "known_contrast": known_rows,
    }

    (out_dir / "novel_disease_eval.json").write_text(json.dumps(report, indent=2) + "\n")
    (out_dir / "NOVEL_DISEASES.md").write_text(_md(report))
    return report


def _md(report: dict[str, Any]) -> str:
    lines = [
        "# Novel diseases (no VOC background)",
        "",
        f"**Passed: {report['passed']}**",
        "",
        f"- novel queries: {report['n_novel']}",
        f"- unresolved (custom::): {report['n_unresolved']}/{report['n_novel']}",
        f"- near-healthy (peak |log2fc| < 0.45): {report['n_near_healthy']}/{report['n_novel']}",
        f"- mean peak |log2fc|: {report['mean_peak_abs_log2fc']:.3f}",
        f"- mean confidence: {report['mean_confidence']:.3f}",
        f"- acetone Δppb (novel mean): {report['novel_mean_acetone_delta_ppb']:.2f}",
        f"- acetone Δppb (T2D): {report['t2d_acetone_delta_ppb']}",
        f"- separation (T2D − novel): {report['acetone_sep_t2d_minus_novel_ppb']:.2f} ppb",
        "",
        "## Gates",
        "",
    ]
    for k, v in (report.get("gates") or {}).items():
        lines.append(f"- {'✓' if v else '✗'} `{k}`")
    lines += ["", "## Novel cases (top VOCs)", ""]
    for r in report.get("novel_cases") or []:
        tops = ", ".join(
            f"{t['voc_id']}({t['log2fc']:+.2f})" for t in (r.get("top_vocs") or [])[:4]
        )
        lines.append(
            f"- **{r['query']['disease']}** → `{r['resolved_disease_id']}` "
            f"unresolved={r['unresolved']} peak|log2fc|={r['peak_abs_log2fc']:.3f} "
            f"near_healthy={r['near_healthy']}"
        )
        lines.append(f"  - top: {tops}")
    lines += ["", "## Known contrast", ""]
    for r in report.get("known_contrast") or []:
        lines.append(
            f"- **{r['query']['disease']}** → `{r['resolved_disease_id']}` "
            f"peak|log2fc|={r['peak_abs_log2fc']:.3f} "
            f"acetone_Δppb={r.get('acetone_delta_ppb')}"
        )
    lines.append("")
    return "\n".join(lines)
