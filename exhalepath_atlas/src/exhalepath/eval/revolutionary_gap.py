"""Full field evaluation: what's left to make ExhalePath one-of-a-kind.

Researches the 2024–2026 breathomics landscape, scores capabilities we already
ship vs competitors, runs the new differentiators (DDx, residualization, claim
ledger, BreathVOC interchange, mechanism↔cohort loop), and writes a roadmap.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import DATA_DIR, PACKAGE_ROOT
from ..gcms.claim_ledger import export_claim_ledger
from ..gcms.differential import evaluate_scidata_differential, rank_diseases_for_vector
from ..gcms.interchange import (
    BREATHVOC_SCHEMA_VERSION,
    export_breathvoc,
    load_and_validate_breathvoc,
    validate_breathvoc_bundle,
    matrix_to_bundle,
)
from ..gcms.patient_matrix import load_mw_patient_matrix
from .mechanism_cohort_loop import evaluate_mechanism_cohort_loop


# Capability scorecard: 0=absent, 1=partial, 2=shipped, 3=field-leading open
FIELD_CAPABILITIES: list[dict[str, Any]] = [
    {
        "id": "mechanism_fusion",
        "name": "Mechanism → VOC multi-model fusion",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 1,
        "scireports_xgboost": 0,
        "note": "20-head Great Disease Stack + hybrid physiology — unique open",
    },
    {
        "id": "patient_nested_cv",
        "name": "Patient-level nested CV + locked splits",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 1,
        "owlstone_omni": 2,
        "scireports_xgboost": 2,
        "note": "SHA256 manifests + optimism gaps; Sci Reports has nested CV but not locked hashes",
    },
    {
        "id": "raw_ms_preprocess",
        "name": "Raw SESI/PTR/GC-MS peak picking",
        "us": 0,
        "breathxplorer": 3,
        "ptairms": 3,
        "owlstone_omni": 3,
        "scireports_xgboost": 0,
        "note": "Largest instrument gap — we consume peak tables, not .mzML/.D",
    },
    {
        "id": "on_breath_blank",
        "name": "On-breath vs blank / MSI identification",
        "us": 1,
        "breathxplorer": 1,
        "ptairms": 2,
        "owlstone_omni": 3,
        "scireports_xgboost": 0,
        "note": "Schema + filters exist; Owlstone VOC Atlas is the commercial gold standard",
    },
    {
        "id": "multiclass_ddx",
        "name": "Cross-disease differential diagnosis",
        "us": 2,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 1,
        "scireports_xgboost": 3,
        "note": "NEW: mechanism DDx + nested multinomial on Sci Data; XGBoost still wins fit-on-cohort",
    },
    {
        "id": "confounder_residual",
        "name": "Confounder residualization (age/sex/smoking)",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 1,
        "scireports_xgboost": 0,
        "note": "NEW: residualize + re-score; ST003200 smoking probes already shipped",
    },
    {
        "id": "claim_ledger",
        "name": "Evidence-graded VOC↔disease claim ledger",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 2,
        "scireports_xgboost": 0,
        "note": "NEW: quantified/directional/prior grades with DOIs — open citeable graph",
    },
    {
        "id": "interchange_schema",
        "name": "Open BreathVOC interchange format",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 1,
        "owlstone_omni": 1,
        "scireports_xgboost": 0,
        "note": "NEW: BreathVOC-1.1 JSON + validator; MetaboLights scaffold already shipped",
    },
    {
        "id": "mechanism_cohort_loop",
        "name": "Mechanism ↔ cohort closed-loop eval",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 0,
        "scireports_xgboost": 0,
        "note": "NEW: unique — packs scored on locked patient matrices with optimism gaps",
    },
    {
        "id": "external_loso",
        "name": "Multi-cohort leave-one-study-out",
        "us": 2,
        "breathxplorer": 0,
        "ptairms": 1,
        "owlstone_omni": 2,
        "scireports_xgboost": 0,
        "note": "Harness + CSIRO labels shipped; blocked on second open intensity table",
    },
    {
        "id": "paper_tripod",
        "name": "TRIPOD+AI / paper-pack automation",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 1,
        "scireports_xgboost": 1,
        "note": "One-zip paper pack + Methods stubs — field-leading open",
    },
    {
        "id": "nl_clinical",
        "name": "Natural-language clinical → VOC prediction",
        "us": 3,
        "breathxplorer": 0,
        "ptairms": 0,
        "owlstone_omni": 0,
        "scireports_xgboost": 0,
        "note": "voc ask / patient / stack — unique product surface",
    },
]


def _scorecard_table() -> dict[str, Any]:
    us_total = sum(c["us"] for c in FIELD_CAPABILITIES)
    max_total = 3 * len(FIELD_CAPABILITIES)
    competitors = {
        "breathxplorer": sum(c["breathxplorer"] for c in FIELD_CAPABILITIES),
        "ptairms": sum(c["ptairms"] for c in FIELD_CAPABILITIES),
        "owlstone_omni": sum(c["owlstone_omni"] for c in FIELD_CAPABILITIES),
        "scireports_xgboost_baseline": sum(
            c["scireports_xgboost"] for c in FIELD_CAPABILITIES
        ),
    }
    return {
        "scale": "0=absent 1=partial 2=shipped 3=field-leading open",
        "n_capabilities": len(FIELD_CAPABILITIES),
        "exhalepath_total": us_total,
        "exhalepath_pct": round(100 * us_total / max_total, 1),
        "competitors": competitors,
        "rows": FIELD_CAPABILITIES,
        "verdict": (
            "ExhalePath already leads the *open research-reasoning* axis "
            "(mechanism fusion, NL, paper packs, claim ledger, closed loop). "
            "It does not lead instrument preprocessing or commercial on-breath "
            "identification (BreathXplorer / ptairMS / Owlstone). Revolutionary "
            "path = deepen the reasoning axis + bridge to peak tables, not chase "
            "fit-on-cohort AUROC."
        ),
    }


def _roadmap() -> list[dict[str, Any]]:
    return [
        {
            "priority": 1,
            "theme": "Second open intensity cohort (malaria or pulmonary)",
            "why": "LOSO / external validation is the #1 credibility gap vs closed multi-site e-nose trials",
            "effort": "data + adapter",
            "blocks": "CSIRO intensity (MassHunter export) or Malawi JID deposit",
            "status": "labels_ready_intensity_pending",
        },
        {
            "priority": 2,
            "theme": "Raw-spectrum bridge (mzML / Agilent .D → BreathVOC)",
            "why": "Without peak picking we cannot ingest the majority of public deposits",
            "effort": "large engineering; partner BreathXplorer/ptairMS rather than reimplement",
            "status": "not_started",
        },
        {
            "priority": 3,
            "theme": "Claim-ledger-grounded NL answers",
            "why": "voc ask should cite quantified DOIs, not only atlas priors",
            "effort": "medium; extend nl/ with ledger retrieval",
            "status": "ledger_shipped_nl_ungrounded",
        },
        {
            "priority": 4,
            "theme": "On-breath MSI enrichment vs Owlstone/NIST",
            "why": "Identification confidence is what converts research → translational trust",
            "effort": "medium; FeatureAnnotation already has MSI fields",
            "status": "schema_ready",
        },
        {
            "priority": 5,
            "theme": "Prospective multi-site protocol template",
            "why": "Field failures (e-nose CRC external validation) are protocol failures",
            "effort": "docs + locked split federation; not more ST000883 tuning",
            "status": "partial_multisite_literature_only",
        },
        {
            "priority": 6,
            "theme": "Do not: hill-climb ST000883 or Sci Data fit-on-cohort AUROC",
            "why": "Sci Reports already owns the black-box ceiling (~0.998); our wedge is transferable mechanism DDx",
            "effort": "discipline",
            "status": "documented_stop_chasing",
        },
    ]


def evaluate_revolutionary_gap(
    *,
    out_dir: Path | None = None,
    quick: bool = False,
) -> dict[str, Any]:
    """Run scorecard + new differentiators; archive under knowledge/revolutionary."""
    out_dir = Path(out_dir or "runs/revolutionary_gap")
    out_dir.mkdir(parents=True, exist_ok=True)

    scorecard = _scorecard_table()
    roadmap = _roadmap()

    # 1) Claim ledger
    ledger_path = export_claim_ledger(out_dir / "CLAIM_LEDGER.json")
    ledger = json.loads(ledger_path.read_text())

    # 2) BreathVOC interchange roundtrip on ST000883
    matrix = load_mw_patient_matrix("ST000883", feature_map="malaria_lit")
    bv_path = export_breathvoc(matrix, out_dir / "breathvoc" / "ST000883.breathvoc.json")
    loaded, _ = load_and_validate_breathvoc(bv_path)
    interchange = {
        "schema_version": BREATHVOC_SCHEMA_VERSION,
        "exported": str(bv_path),
        "roundtrip_ok": loaded.matrix.shape == matrix.matrix.shape
        and list(loaded.matrix.columns) == list(matrix.matrix.columns),
        "n_subjects": int(matrix.matrix.shape[0]),
        "n_vocs": int(matrix.matrix.shape[1]),
        "validation_errors": validate_breathvoc_bundle(matrix_to_bundle(matrix)),
    }

    # 3) Sci Data differential DDx
    if quick:
        ddx = {"status": "skipped_quick"}
    else:
        ddx = evaluate_scidata_differential(seed=42)

    # 4) Example DDx ranking for a synthetic malaria-ish vector
    demo_vec = {
        "alpha_pinene": 0.8,
        "delta_3_carene": 0.6,
        "cyclohexanone": 0.5,
        "tridecane": 0.4,
        "acetone": -0.2,
        "isoprene": -0.3,
        "pentane": 0.3,
        "hexanal": 0.2,
    }
    ranking = rank_diseases_for_vector(
        demo_vec,
        ["malaria", "asthma", "copd", "heart_failure", "lung_adenocarcinoma"],
        source="hybrid",
    )

    # 5) Mechanism ↔ cohort loop
    if quick:
        loop = {"status": "skipped_quick"}
    else:
        loop = evaluate_mechanism_cohort_loop(out_dir=out_dir / "mechanism_loop", seed=42)

    report: dict[str, Any] = {
        "title": "Revolutionary gap evaluation — ExhalePath vs breathomics field",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty_bar": (
            "Research / hypothesis-generation enablement. Not clinical diagnostic SOTA, "
            "not FDA claims, not a substitute for prospective multi-site trials."
        ),
        "field_landscape_2026": {
            "open_preprocess": [
                "BreathXplorer (SESI-HRMS peak recognition)",
                "ptairMS (PTR-TOF-MS Bioconductor)",
            ],
            "commercial_closed": [
                "Owlstone Breath Biopsy OMNI + VOC Atlas (on-breath MSI)",
            ],
            "fit_on_cohort_ml": [
                "Sci Reports 2025 XGBoost on Sci Data GC-MS (macro-AUC≈0.998 nested)",
                "PTR multiclass lung-disease ensembles (IJMS 2026)",
                "e-nose multi-site CRC — external validation often fails",
            ],
            "reporting_standard": "TRIPOD+AI (2024); nested CV consensus in JBR breath-ML",
            "exhalepath_wedge": (
                "Only open stack that joins mechanism-fused VOC reasoning + patient-level "
                "diagnostic harness + evidence ledger + NL clinical front-door."
            ),
        },
        "scorecard": scorecard,
        "live_differentiators": {
            "claim_ledger": {
                "n_claims": ledger.get("n_claims"),
                "by_evidence_grade": ledger.get("by_evidence_grade"),
                "path": str(ledger_path),
            },
            "breathvoc_interchange": interchange,
            "scidata_differential": ddx,
            "demo_ddx_ranking": ranking,
            "mechanism_cohort_loop": loop,
        },
        "roadmap": roadmap,
        "what_revolutionary_means_here": [
            "Be the open *reasoning layer* of breathomics — not another peak picker",
            "Make mechanism hypotheses falsifiable on locked patient matrices",
            "Make every VOC↔disease edge citeable with an evidence grade",
            "Make BreathVOC the interchange labs actually ship between instruments and papers",
            "Refuse AUROC theater on n≈35 / fit-on-cohort ceilings",
        ],
    }

    (out_dir / "REVOLUTIONARY_GAP.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    (out_dir / "REVOLUTIONARY_GAP.md").write_text(_md(report))

    # Archive into knowledge (repo + package mirrors)
    for root in (
        DATA_DIR / "knowledge" / "revolutionary",
        PACKAGE_ROOT / "data" / "knowledge" / "revolutionary",
        Path("src/exhalepath/data/knowledge/revolutionary"),
    ):
        try:
            root.mkdir(parents=True, exist_ok=True)
            for name in (
                "REVOLUTIONARY_GAP.md",
                "REVOLUTIONARY_GAP.json",
                "CLAIM_LEDGER.json",
                "CLAIM_LEDGER.md",
            ):
                src = out_dir / name
                if src.exists():
                    shutil.copy2(src, root / name)
            bv_src = out_dir / "breathvoc"
            if bv_src.exists():
                dest = root / "breathvoc"
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(bv_src, dest)
            loop_src = out_dir / "mechanism_loop"
            if loop_src.exists():
                dest = root / "mechanism_loop"
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(loop_src, dest)
        except Exception:  # noqa: BLE001
            continue

    report["archived_under"] = "data/knowledge/revolutionary/"
    return report


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


def _md(report: dict[str, Any]) -> str:
    sc = report.get("scorecard") or {}
    ddx = (report.get("live_differentiators") or {}).get("scidata_differential") or {}
    loop = (report.get("live_differentiators") or {}).get("mechanism_cohort_loop") or {}
    ledger = (report.get("live_differentiators") or {}).get("claim_ledger") or {}
    lines = [
        "# Revolutionary gap evaluation",
        "",
        f"Generated: {report.get('generated_utc')}",
        "",
        f"> {report.get('honesty_bar')}",
        "",
        "## Verdict",
        "",
        sc.get("verdict") or "",
        "",
        f"**ExhalePath scorecard:** {sc.get('exhalepath_total')} / "
        f"{3 * (sc.get('n_capabilities') or 12)} ({sc.get('exhalepath_pct')}%)",
        "",
        "### Competitor totals (same scale)",
        "",
    ]
    for k, v in (sc.get("competitors") or {}).items():
        lines.append(f"- `{k}`: {v}")
    lines += ["", "## Capability matrix", "", "| Capability | Us | BreathXplorer | ptairMS | Owlstone | SciRep XGB |", "|---|---:|---:|---:|---:|---:|"]
    for r in sc.get("rows") or []:
        lines.append(
            f"| {r['name']} | {r['us']} | {r['breathxplorer']} | {r['ptairms']} | "
            f"{r['owlstone_omni']} | {r['scireports_xgboost']} |"
        )
    lines += ["", "## Live differentiators (this run)", ""]
    lines.append(
        f"- Claim ledger: **{ledger.get('n_claims')}** claims — {ledger.get('by_evidence_grade')}"
    )
    bv = (report.get("live_differentiators") or {}).get("breathvoc_interchange") or {}
    lines.append(
        f"- BreathVOC interchange: `{bv.get('schema_version')}` "
        f"roundtrip_ok={bv.get('roundtrip_ok')}"
    )
    if ddx.get("fit_on_cohort_ceiling"):
        ceil = ddx["fit_on_cohort_ceiling"]
        mech = ddx.get("mechanism_ddx_hybrid") or {}
        resid = ddx.get("age_sex_residualized") or {}
        lines.append(
            f"- Sci Data DDx: fit-on-cohort macro AUROC **{_pct(ceil.get('macro_ovr_auroc'))}** "
            f"vs mechanism hybrid **{_pct(mech.get('macro_ovr_auroc'))}** "
            f"(top1={_pct(mech.get('top1_accuracy'))})"
        )
        if isinstance(resid, dict) and resid.get("mean_covariate_r2") is not None:
            lines.append(
                f"- Age/sex residualization: mean VOC R²={resid.get('mean_covariate_r2'):.3f}; "
                f"residual mechanism macro="
                f"{_pct((resid.get('mechanism_hybrid') or {}).get('macro_ovr_auroc'))}"
            )
    if loop.get("studies"):
        lines.append("- Mechanism↔cohort loop:")
        for st in loop["studies"]:
            hyb = (st.get("sources") or {}).get("hybrid") or {}
            lines.append(
                f"  - `{st.get('study_id')}` hybrid nested={_pct(hyb.get('nested_auroc'))}"
            )
    lines += ["", "## Roadmap (what remains)", ""]
    for item in report.get("roadmap") or []:
        lines.append(
            f"{item.get('priority')}. **{item.get('theme')}** — {item.get('why')} "
            f"(`{item.get('status')}`)"
        )
    lines += ["", "## What 'revolutionary' means here", ""]
    for w in report.get("what_revolutionary_means_here") or []:
        lines.append(f"- {w}")
    lines.append("")
    return "\n".join(lines)


__all__ = ["evaluate_revolutionary_gap", "FIELD_CAPABILITIES"]
