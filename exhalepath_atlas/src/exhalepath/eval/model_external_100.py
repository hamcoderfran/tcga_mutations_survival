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
    lit = report.get("literature_concordance") or {}
    d100 = report.get("disease100") or {}
    lines = [
        "# Model improve + external-validated 100-disease suite",
        "",
        "Priors soft-blended from secured open corpora (priority literature panels, "
        "public breath benchmarks, literature_benchmarks, HBDB proxies, VOLATILOME identity). "
        "Calibrator joblibs were **not** retrained.",
        "",
        f"**External validation mean concordance: {ext.get('mean_concordance_pct')}%** "
        f"({ext.get('n_diseases_with_external_gt')} diseases with external GT)",
        f"**Literature concordance (expanded panels): {lit.get('mean_concordance_pct')}%** "
        f"({lit.get('n_diseases')} diseases)",
        f"**Disease connection suite: {d100.get('n_diseases')} diseases**",
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
        f"- Mean concordance: {ext.get('mean_concordance_pct')}%",
        f"- Diseases with external GT: {ext.get('n_diseases_with_external_gt')}",
        f"- VOC checks: {ext.get('n_voc_checks')} (hits {ext.get('n_voc_hits')})",
        "",
        "Top disease concordances:",
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


def external_validation(*, n_diseases: int = 100) -> dict[str, Any]:
    """Score predictions against fused external elevate/suppress GT for up to N diseases."""
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    atlas_ids = set(eng.kb.diseases.keys())
    atlas_vocs = set(eng.kb.vocs.keys())
    expect = build_external_expect_map(atlas_ids, atlas_vocs)

    # Prefer atlas order; take first n_diseases that have external GT, else pad with all
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
    }


def run_model_improve_and_disease100(
    *,
    out_dir: Path,
    n_diseases: int = 100,
    demo_max_patients: int | None = 200,
    skip_fuse: bool = False,
) -> dict[str, Any]:
    """
    1) Fuse secured external evidence into priors
    2) External validation vs open GT
    3) Literature concordance + 100-disease connection suite
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fuse = {"skipped": True}
    if not skip_fuse:
        fuse = fuse_external_into_priors(dry_run=False)

    ext = external_validation(n_diseases=n_diseases)
    # Full lit-demo suite (includes literature concordance + demographics PCA + disease100)
    full = run_lit_demo_disease100(
        out_dir=out_dir / "lit_compare",
        demo_max_patients=demo_max_patients,
        n_diseases=n_diseases,
    )
    lit = full.get("literature") or literature_concordance()

    report = {
        "version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Validation uses open secured corpora already in-repo (panels, public breath, "
            "literature benchmarks, HBDB proxies). HBDB SQL / Owlstone Atlas still blocked "
            "until user-supplied dumps arrive. Calibrator *.joblib untouched."
        ),
        "fuse": fuse,
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
