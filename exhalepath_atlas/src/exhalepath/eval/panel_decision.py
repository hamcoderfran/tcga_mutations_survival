"""Panel go/no-go decision memo (Pillar 3).

Consumes existing nested AUROC, optimism gap, claim grades, outside-prior
overlay, and confounder strata presence → research decision enum.

Honesty: proceed_research | hold | stop for R&D diligence only — never a
clinical clearance recommendation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from ..config import DATA_DIR, PACKAGE_ROOT
from ..gcms.claim_ledger import build_claim_ledger
from ..oem.kit import cite_vocs_for_disease
from ..viz.analysis_export import literature_overlay

Decision = Literal["proceed_research", "hold", "stop"]


def _load_patient_diag(study_id: str = "ST000883") -> dict[str, Any] | None:
    """Load the best available patient-diagnostic report for a study.

    Prefer ``*_hybrid`` archives, then highest nested AUROC across DATA_DIR and
    repo ``data/`` mirrors (packaged wheel vs editable checkout).
    """
    roots = [
        DATA_DIR / "knowledge" / "gcms_diagnostic",
        PACKAGE_ROOT / "data" / "knowledge" / "gcms_diagnostic",
        PACKAGE_ROOT / "src" / "exhalepath" / "data" / "knowledge" / "gcms_diagnostic",
    ]
    candidates: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for suffix in (f"{study_id}_hybrid", f"{study_id}_literature", study_id, f"{study_id}_stack"):
            p = root / suffix / "PATIENT_DIAGNOSTIC_REPORT.json"
            if p.exists():
                candidates.append(p)
        # also catch any other study_id_* dirs
        for p in sorted(root.glob(f"{study_id}*/PATIENT_DIAGNOSTIC_REPORT.json")):
            if p not in candidates:
                candidates.append(p)

    best: dict[str, Any] | None = None
    best_score = -1.0
    for p in candidates:
        try:
            doc = json.loads(p.read_text())
        except Exception:
            continue
        nested = doc.get("nested") or {}
        score = nested.get("mean_test_auroc")
        try:
            score_f = float(score) if score is not None else -1.0
        except (TypeError, ValueError):
            score_f = -1.0
        # Prefer hybrid path on ties
        tie_bonus = 0.001 if "hybrid" in str(p.parent.name) else 0.0
        ranked = score_f + tie_bonus
        if best is None or ranked > best_score:
            best = doc
            best_score = ranked
    return best


def _confounder_available() -> dict[str, Any]:
    path = DATA_DIR / "knowledge" / "confounder_ptr" / "CONFOUNDER_PTR.json"
    if not path.exists():
        path = PACKAGE_ROOT / "data" / "knowledge" / "confounder_ptr" / "CONFOUNDER_PTR.json"
    if not path.exists():
        return {"available": False}
    doc = json.loads(path.read_text())
    return {
        "available": True,
        "path": str(path),
        "has_smoking": "smoking" in json.dumps(doc).lower(),
        "has_sex_or_age": ("sex" in json.dumps(doc).lower()) or ("age" in json.dumps(doc).lower()),
    }


def decide_panel(
    *,
    disease_id: str,
    study_id: str | None = None,
    predicted_log2fc: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Build a research go/no-go memo for a disease panel."""
    study_id = study_id or ("ST000883" if disease_id == "malaria" else "ST000883")
    diag = _load_patient_diag(study_id) if disease_id in {"malaria", "heart_failure"} or study_id else None
    if disease_id == "heart_failure":
        diag = _load_patient_diag("ST000587")
        study_id = "ST000587"

    nested_auroc = None
    optimism_gap = None
    n_subjects = None
    split_sha = None
    if diag:
        nested = diag.get("nested") or {}
        nested_auroc = nested.get("mean_test_auroc")
        optimism_gap = nested.get("optimism_gap")
        n_subjects = (diag.get("metrics") or {}).get("n_subjects")
        split_sha = nested.get("content_sha256")

    ledger = build_claim_ledger(disease_ids=[disease_id], include_atlas_priors=True)
    claims = [c for c in (ledger.get("claims") or []) if c.get("disease_id") == disease_id]
    grade_counts: dict[str, int] = {}
    for c in claims:
        g = str(c.get("evidence_grade") or "unknown")
        grade_counts[g] = grade_counts.get(g, 0) + 1
    n_quant = grade_counts.get("quantified", 0)
    n_dir = grade_counts.get("directional_only", 0)
    n_prior = grade_counts.get("atlas_prior", 0)

    overlay = None
    outside_prior_acc = None
    if predicted_log2fc is not None:
        overlay = literature_overlay(disease_id, predicted_log2fc)
        outside_prior_acc = overlay.get("directional_accuracy_outside_prior")
    else:
        # Use prior signs as a stand-in only to surface circularity risk in memo
        from ..knowledge.loader import clear_knowledge_cache, default_knowledge

        clear_knowledge_cache()
        kb = default_knowledge()
        prior = dict((kb.diseases.get(disease_id) or {}).get("voc_log2fc_prior") or {})
        if prior:
            overlay = literature_overlay(disease_id, prior, prior_log2fc=prior)
            outside_prior_acc = overlay.get("directional_accuracy_outside_prior")

    confounders = _confounder_available()

    reasons: list[str] = []
    decision: Decision = "proceed_research"

    if n_quant == 0 and n_dir == 0:
        decision = "stop"
        reasons.append("No quantified or directional literature claims in ledger — thin prior only")
    elif nested_auroc is not None and float(nested_auroc) < 0.55 and n_quant < 2:
        decision = "hold"
        reasons.append(
            f"Nested AUROC={float(nested_auroc):.3f} near chance with sparse quantified claims"
        )
    elif optimism_gap is not None and float(optimism_gap) > 0.15:
        decision = "hold"
        reasons.append(f"Optimism gap={float(optimism_gap):.3f} suggests leakage risk")
    elif n_prior > 0 and n_quant == 0 and outside_prior_acc is None:
        decision = "hold"
        reasons.append("Claims are atlas_prior-heavy; outside-prior concordance unavailable")
    else:
        reasons.append("Enough literature grades and/or nested metrics to proceed as research diligence")

    if not confounders.get("available"):
        reasons.append("Confounder PTR report missing — report smoking/age/sex strata when possible")

    top_vocs = []
    for c in claims:
        if c.get("evidence_grade") in {"quantified", "directional_only", "mixed"}:
            top_vocs.append(c.get("voc_id"))
    top_vocs = [v for v in top_vocs if v][:12]
    citations = cite_vocs_for_disease(disease_id, top_vocs) if top_vocs else {}

    return {
        "schema_version": "PanelDecision-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "disease_id": disease_id,
        "study_id": study_id if diag else None,
        "decision": decision,
        "reasons": reasons,
        "metrics": {
            "nested_auroc": nested_auroc,
            "optimism_gap": optimism_gap,
            "n_subjects": n_subjects,
            "split_sha256": split_sha,
            "claim_grade_counts": grade_counts,
            "n_quantified": n_quant,
            "n_directional_only": n_dir,
            "n_atlas_prior": n_prior,
            "directional_accuracy_outside_prior": outside_prior_acc,
            "confounders": confounders,
        },
        "cited_vocs": citations,
        "honesty": (
            "Research diligence go/no-go only (proceed_research | hold | stop). "
            "NOT clinical clearance, FDA advice, or a diagnostic recommendation."
        ),
        "clinical_claim": False,
        "overlay_n_compared": (overlay or {}).get("n_compared"),
    }


def render_panel_decision_md(report: dict[str, Any]) -> str:
    m = report.get("metrics") or {}
    lines = [
        f"# Panel decision — {report.get('disease_id')}",
        "",
        f"Generated: {report.get('generated_utc')}",
        "",
        f"> {report.get('honesty')}",
        "",
        f"**Decision: `{report.get('decision')}`**",
        "",
        "## Reasons",
        "",
    ]
    for r in report.get("reasons") or []:
        lines.append(f"- {r}")
    lines += [
        "",
        "## Diligence metrics",
        "",
        f"- nested AUROC: {m.get('nested_auroc')}",
        f"- optimism gap: {m.get('optimism_gap')}",
        f"- n subjects: {m.get('n_subjects')}",
        f"- split sha256: `{(m.get('split_sha256') or '—')[:16]}…`",
        f"- claim grades: {m.get('claim_grade_counts')}",
        f"- outside-prior directional accuracy: {m.get('directional_accuracy_outside_prior')}",
        f"- confounder report available: {(m.get('confounders') or {}).get('available')}",
        "",
        "## Cited VOCs (ledger)",
        "",
    ]
    for vid, meta in (report.get("cited_vocs") or {}).items():
        lines.append(
            f"- `{vid}`: grade={meta.get('evidence_grade')} doi={meta.get('doi')}"
        )
    lines += [
        "",
        "## Next commands",
        "",
        "```bash",
        "voc eval-leaderboard",
        "voc export-claim-ledger --disease " + str(report.get("disease_id")),
        "voc eval-patient-diagnostic --study ST000883 --signature hybrid",
        "```",
        "",
    ]
    return "\n".join(lines)


def write_panel_decision(
    *,
    disease_id: str,
    out_dir: Path = Path("runs/panel_decision"),
    study_id: str | None = None,
) -> dict[str, Any]:
    report = decide_panel(disease_id=disease_id, study_id=study_id)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"PANEL_DECISION_{disease_id}"
    (out_dir / f"{stem}.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    (out_dir / f"{stem}.md").write_text(render_panel_decision_md(report))
    return report


__all__ = ["decide_panel", "write_panel_decision", "render_panel_decision_md"]
