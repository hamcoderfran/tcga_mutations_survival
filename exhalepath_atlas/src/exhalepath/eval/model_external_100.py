"""Improve model from secured external evidence + re-run 100-disease validation."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR, PACKAGE_ROOT
from ..knowledge.external_evidence import (
    build_external_expect_map,
    fuse_external_into_priors,
)
from ..knowledge.loader import clear_knowledge_cache
from .lit_compare import (
    _direction_score,
    _voc_vector,
    literature_concordance,
    run_disease100_suite,
    run_lit_demo_disease100,
)


def _md(report: dict[str, Any]) -> str:
    fuse = report.get("fuse") or {}
    ext = report.get("external_validation") or {}
    held = report.get("held_out_literature_benchmarks") or {}
    lit = report.get("literature_concordance") or {}
    d100 = report.get("disease100") or {}
    audit = report.get("literature_benchmarks_audit") or {}
    lines = [
        "# Model improve + external-validated 100-disease suite",
        "",
        report.get("honesty") or "",
        "",
        f"**Held-out literature_benchmarks concordance: {held.get('mean_concordance_pct')}%** "
        f"({held.get('n_diseases_with_external_gt')} diseases)",
        f"**In-sample external concordance: {ext.get('mean_concordance_pct')}%** "
        f"({ext.get('n_diseases_with_external_gt')} diseases with external GT)",
        f"**Literature panel concordance: {lit.get('mean_concordance_pct')}%** "
        f"({lit.get('n_diseases')} diseases)",
        f"**Disease connection suite: {d100.get('n_diseases')} diseases**",
        f"**Audit directional accuracy: {audit.get('directional_accuracy')}** "
        f"(min-fold pass {audit.get('min_fold_pass_rate')})",
        "",
        "## Fuse summary",
        "",
        f"- Diseases updated: {fuse.get('n_diseases_updated')}",
        f"- VOC prior updates: {fuse.get('n_voc_prior_updates')}",
        f"- Evidence edges: {(fuse.get('stats') or {}).get('n_voc_edges')}",
        f"- Europe PMC DOIs (provenance): {(fuse.get('stats') or {}).get('europepmc_dois')}",
        f"- Expanded MW studies cataloged: {(fuse.get('stats') or {}).get('expanded_mw_studies')}",
        f"- VOC identity enrichments: {(fuse.get('voc_identity_enrichment') or {}).get('updated')}",
        "",
        "## External validation (vs open data)",
        "",
        f"- Held-out lit-bench mean concordance: {held.get('mean_concordance_pct')}%",
        f"- In-sample mean concordance: {ext.get('mean_concordance_pct')}%",
        f"- Diseases with external GT (in-sample): {ext.get('n_diseases_with_external_gt')}",
        f"- VOC checks: {ext.get('n_voc_checks')} (hits {ext.get('n_voc_hits')})",
        "",
        "Top disease concordances (in-sample):",
        "",
    ]
    by = sorted(
        ext.get("by_disease") or [],
        key=lambda r: (-(r.get("concordance") or -1), r.get("disease_id") or ""),
    )
    for row in by[:15]:
        lines.append(
            f"- `{row.get('disease_id')}`: "
            f"{None if row.get('concordance') is None else round(100 * row['concordance'], 1)}% "
            f"({row.get('n_hit')}/{row.get('n_checked')})"
        )
    lines += [
        "",
        "## Disease-100 connections (sample)",
        "",
    ]
    for row in (d100.get("interesting_connections") or [])[:8]:
        lines.append(
            f"- {row.get('disease_a')} ↔ {row.get('disease_b')} "
            f"(cos={row.get('cosine', float('nan')):.3f})"
        )
    lines += [
        "",
        "## Artifacts",
        "",
        "- `data/knowledge/EXTERNAL_EVIDENCE_FUSE.json`",
        "- `data/knowledge/MODEL_EXTERNAL_100DISEASE.json`",
        "- `data/knowledge/MODEL_EXTERNAL_100DISEASE.md`",
        "- `data/knowledge/lit_compare/` (connection matrices / figures)",
        "",
        f"Generated: {report.get('generated_at')}",
        "",
    ]
    return "\n".join(lines)


def external_validation(
    *,
    n_diseases: int = 100,
    include_lit_bench: bool = True,
    include_public_breath: bool = True,
    include_panels: bool = True,
    include_hbdb: bool = True,
) -> dict[str, Any]:
    """Score predictions against fused external elevate/suppress GT for up to N diseases."""
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    atlas_ids = set(eng.kb.diseases.keys())
    atlas_vocs = set(eng.kb.vocs.keys())
    expect = build_external_expect_map(
        atlas_ids,
        atlas_vocs,
        include_lit_bench=include_lit_bench,
        include_public_breath=include_public_breath,
        include_panels=include_panels,
        include_hbdb=include_hbdb,
    )

    ordered = sorted(eng.kb.diseases.keys())[:n_diseases]
    results = []
    for did in ordered:
        exp = expect.get(did)
        if not exp:
            continue
        loc = (eng.kb.diseases.get(did) or {}).get("default_site") or "systemic"
        report = eng.predict(
            disease=did,
            location=loc,
            top_n=50,
            explain=False,
            smoking_status="never",
            age_years=60,
            sex="male",
        )
        vec = _voc_vector(report)
        scored = _direction_score(vec, exp.get("elevate") or [], exp.get("suppress") or [])
        results.append(
            {
                "disease_id": report.disease_id,
                "concordance": scored["concordance"],
                "n_hit": scored["n_hit"],
                "n_checked": scored["n_checked"],
                "n_elevate": len(exp.get("elevate") or []),
                "n_suppress": len(exp.get("suppress") or []),
                "refs": (exp.get("refs") or [])[:8],
            }
        )

    conc = [r["concordance"] for r in results if r["concordance"] is not None]
    return {
        "n_diseases_with_external_gt": len(results),
        "n_external_expect_diseases_available": len(expect),
        "mean_concordance": float(np.mean(conc)) if conc else None,
        "mean_concordance_pct": round(100 * float(np.mean(conc)), 2) if conc else None,
        "n_voc_checks": int(sum(r.get("n_checked") or 0 for r in results)),
        "n_voc_hits": int(sum(r.get("n_hit") or 0 for r in results)),
        "by_disease": results,
        "gt_sources": {
            "include_lit_bench": include_lit_bench,
            "include_public_breath": include_public_breath,
            "include_panels": include_panels,
            "include_hbdb": include_hbdb,
        },
    }


def run_model_improve_and_disease100(
    *,
    out_dir: Path,
    n_diseases: int = 100,
    demo_max_patients: int | None = 200,
    skip_fuse: bool = False,
) -> dict[str, Any]:
    """
    1) Fuse secured external evidence into priors (panels/public/HBDB; lit_bench held out of fuse)
    2) Held-out validation vs literature_benchmarks
    3) In-sample external validation vs all open GT
    4) Final fuse including lit_bench for production priors
    5) Literature concordance + 100-disease connection suite
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fuse_holdout = {"skipped": True}
    fuse = {"skipped": True}
    held_out = None
    if not skip_fuse:
        # Stage A: fuse without literature_benchmarks so they stay held-out
        fuse_holdout = fuse_external_into_priors(dry_run=False, include_lit_bench=False)
        held_out = external_validation(
            n_diseases=n_diseases,
            include_lit_bench=True,
            include_public_breath=False,
            include_panels=False,
            include_hbdb=False,
        )
        # Stage B: full fuse for production priors
        fuse = fuse_external_into_priors(dry_run=False, include_lit_bench=True)

    ext = external_validation(n_diseases=n_diseases)
    # audit literature benchmarks (case-level)
    from .audit import evaluate_literature_benchmarks

    lit_audit = evaluate_literature_benchmarks()

    full = run_lit_demo_disease100(
        out_dir=out_dir / "lit_compare",
        demo_max_patients=demo_max_patients,
        n_diseases=n_diseases,
    )
    lit = full.get("literature") or literature_concordance()

    report = {
        "version": "1.1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "In-sample external concordance uses open corpora that also inform prior fusion "
            "(panels / public breath / HBDB proxies) and can look near-perfect. "
            "Held-out literature_benchmarks concordance is the less circular check. "
            "HBDB SQL / Owlstone still blocked until user dumps arrive. "
            "Calibrator *.joblib untouched."
        ),
        "fuse_holdout_stage": fuse_holdout,
        "fuse": fuse,
        "held_out_literature_benchmarks": held_out,
        "literature_benchmarks_audit": {
            "n_cases": lit_audit.get("n_cases"),
            "directional_accuracy": lit_audit.get("directional_accuracy"),
            "min_fold_pass_rate": lit_audit.get("min_fold_accuracy"),
            "case_pass_rate": lit_audit.get("case_pass_rate"),
        },
        "external_validation": ext,
        "literature_concordance": lit,
        "disease100": full.get("disease100"),
        "demographics_pca": {
            "silhouette": (full.get("demographics_pca") or {}).get("silhouette"),
        },
        "bugs_smoke": full.get("bugs"),
    }

    md = _md(report)
    for know in {KNOWLEDGE_DIR, PACKAGE_ROOT / "data" / "knowledge"}:
        if know.exists():
            (know / "MODEL_EXTERNAL_100DISEASE.json").write_text(
                json.dumps(report, indent=2) + "\n"
            )
            (know / "MODEL_EXTERNAL_100DISEASE.md").write_text(md)
    (out_dir / "MODEL_EXTERNAL_100DISEASE.json").write_text(json.dumps(report, indent=2) + "\n")
    (out_dir / "MODEL_EXTERNAL_100DISEASE.md").write_text(md)
    return report
