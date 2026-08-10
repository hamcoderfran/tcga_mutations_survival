"""Mental-health literature panel eval with prior-masked (de-circularized) overlay.

For each MH disease with a measured literature panel:

1. **Raw overlay** — predict with full atlas prior, compare to panel
   (often circular when panel VOCs were copied into the prior).
2. **Panel-masked prior** — temporarily remove panel VOC keys from
   ``voc_log2fc_prior``, re-predict, compare to the same panel.

The masked score is the honest research metric for "does mechanism /
physiology still recover Magdeburg directions without peeking at the
panel VOC prior entries?"
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..data.literature_panels import (
    load_literature_panels,
    load_thin_evidence_notes,
    literature_panel_for_disease,
)
from ..knowledge.loader import clear_knowledge_cache
from ..viz.analysis_export import literature_overlay

MH_PANEL_DISEASES = (
    "schizophrenia",
    "major_depressive_disorder",
    "bipolar",
)
MH_THIN_DISEASES = (
    "anxiety",
    "ptsd",
    "adhd",
    "autism_spectrum_disorder",
)


def _predict_vec(engine: ExhaleBiomarkerEngine, disease_id: str) -> dict[str, float]:
    d = engine.kb.diseases.get(disease_id) or {}
    loc = d.get("default_site") or "brain"
    report = engine.predict(
        disease=disease_id,
        location=loc,
        top_n=80,
        explain=False,
        smoking_status="never",
        age_years=40,
        sex="male",
    )
    return {p.voc_id: float(p.log2_fold_change) for p in report.result.bundle.predictions}


def _with_masked_panel_prior(
    engine: ExhaleBiomarkerEngine,
    disease_id: str,
    panel_vocs: set[str],
):
    """Context-manager-like mask: yield after swapping prior, restore after."""

    class _Mask:
        def __enter__(self):
            original = engine.kb.diseases[disease_id]
            masked = dict(original)
            prior = dict(original.get("voc_log2fc_prior") or {})
            removed = {k: prior.pop(k) for k in list(prior) if k in panel_vocs}
            masked["voc_log2fc_prior"] = prior
            masked["_mh_panel_masked"] = True
            masked["_mh_removed_prior_vocs"] = removed
            self.original = original
            self.removed = removed
            engine.kb.diseases[disease_id] = masked
            return self

        def __exit__(self, *exc):
            engine.kb.diseases[disease_id] = self.original
            return False

    return _Mask()


def evaluate_mental_health_literature(
    *,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    """Run raw + panel-masked literature overlays for MH panel diseases."""
    clear_knowledge_cache()
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    panels = {p["disease_id"]: p for p in load_literature_panels()}
    thin = load_thin_evidence_notes()

    cases: list[dict[str, Any]] = []
    for did in MH_PANEL_DISEASES:
        panel = panels.get(did) or literature_panel_for_disease(did)
        if not panel:
            cases.append({"disease_id": did, "skipped": True, "reason": "no_panel"})
            continue
        measured = dict(panel.get("measured_log2fc") or {})
        if not measured:
            cases.append({"disease_id": did, "skipped": True, "reason": "empty_panel"})
            continue

        # raw
        pred_raw = _predict_vec(engine, did)
        prior = dict((engine.kb.diseases.get(did) or {}).get("voc_log2fc_prior") or {})
        ov_raw = literature_overlay(did, pred_raw, prior_log2fc=prior)

        # panel-masked prior
        with _with_masked_panel_prior(engine, did, set(measured)) as mask:
            pred_masked = _predict_vec(engine, did)
            ov_masked = literature_overlay(
                did,
                pred_masked,
                prior_log2fc=dict(
                    (engine.kb.diseases.get(did) or {}).get("voc_log2fc_prior") or {}
                ),
            )

        quant = [
            v
            for v, ev in (panel.get("voc_evidence") or {}).items()
            if (ev or {}).get("evidence") == "quantified"
        ]
        cases.append(
            {
                "disease_id": did,
                "skipped": False,
                "n_panel_vocs": len(measured),
                "n_quantified_vocs": len(quant),
                "quantified_vocs": quant,
                "n_prior_vocs_removed": len(mask.removed),
                "removed_prior_vocs": sorted(mask.removed),
                "raw": {
                    "n_agree": ov_raw.get("n_agree"),
                    "n_compared": ov_raw.get("n_compared"),
                    "directional_accuracy": ov_raw.get("directional_accuracy"),
                    "n_panel_vocs_also_in_prior": ov_raw.get("n_panel_vocs_also_in_prior"),
                    "directional_accuracy_outside_prior": ov_raw.get(
                        "directional_accuracy_outside_prior"
                    ),
                },
                "panel_masked_prior": {
                    "n_agree": ov_masked.get("n_agree"),
                    "n_compared": ov_masked.get("n_compared"),
                    "directional_accuracy": ov_masked.get("directional_accuracy"),
                    "n_panel_vocs_also_in_prior": ov_masked.get("n_panel_vocs_also_in_prior"),
                    "directional_accuracy_outside_prior": ov_masked.get(
                        "directional_accuracy_outside_prior"
                    ),
                },
                "refs": [
                    {"doi": r.get("doi"), "title": r.get("title"), "year": r.get("year")}
                    for r in (panel.get("refs") or [])
                    if r.get("doi")
                ],
            }
        )

    thin_rows = []
    for did in MH_THIN_DISEASES:
        note = thin.get(did) or {}
        d = engine.kb.diseases.get(did) or {}
        thin_rows.append(
            {
                "disease_id": did,
                "has_measured_panel": False,
                "atlas_source": d.get("atlas_source"),
                "n_prior_vocs": len(d.get("voc_log2fc_prior") or {}),
                "note": note.get("note"),
                "refs": note.get("refs") or [],
            }
        )

    scored = [c for c in cases if not c.get("skipped")]
    raw_accs = [
        c["raw"]["directional_accuracy"]
        for c in scored
        if c["raw"].get("directional_accuracy") is not None
    ]
    masked_accs = [
        c["panel_masked_prior"]["directional_accuracy"]
        for c in scored
        if c["panel_masked_prior"].get("directional_accuracy") is not None
    ]

    report = {
        "schema_version": "MentalHealthLitEval-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Raw directional accuracy can be curation-circular when panel VOCs "
            "are also in atlas priors. panel_masked_prior removes those prior "
            "entries before prediction — that score is the de-circularized metric. "
            "Neither score is a clinical AUROC; no MH patient intensity cohort is bundled."
        ),
        "n_panel_diseases": len(scored),
        "mean_raw_directional_accuracy": (sum(raw_accs) / len(raw_accs)) if raw_accs else None,
        "mean_panel_masked_directional_accuracy": (
            (sum(masked_accs) / len(masked_accs)) if masked_accs else None
        ),
        "cases": cases,
        "thin_evidence_conditions": thin_rows,
    }

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "mental_health_lit_eval.json").write_text(json.dumps(report, indent=2) + "\n")
        lines = [
            "# Mental-health literature eval (de-circularized)",
            "",
            f"Generated: {report['generated_utc']}",
            "",
            f"> {report['honesty']}",
            "",
            f"- Mean **raw** directional accuracy: "
            f"{report['mean_raw_directional_accuracy']}",
            f"- Mean **panel-masked prior** directional accuracy: "
            f"{report['mean_panel_masked_directional_accuracy']}",
            "",
            "## Panel diseases",
            "",
        ]
        for c in scored:
            lines.append(
                f"- **{c['disease_id']}**: raw="
                f"{c['raw']['directional_accuracy']} "
                f"({c['raw']['n_agree']}/{c['raw']['n_compared']}); "
                f"masked="
                f"{c['panel_masked_prior']['directional_accuracy']} "
                f"({c['panel_masked_prior']['n_agree']}/{c['panel_masked_prior']['n_compared']}); "
                f"removed_prior={c['n_prior_vocs_removed']}"
            )
        lines += ["", "## Thin-evidence conditions", ""]
        for t in thin_rows:
            lines.append(
                f"- **{t['disease_id']}**: no measured panel · "
                f"`atlas_source={t.get('atlas_source')}`"
            )
        lines.append("")
        (out_dir / "mental_health_lit_eval.md").write_text("\n".join(lines))

    return report


__all__ = ["evaluate_mental_health_literature", "MH_PANEL_DISEASES", "MH_THIN_DISEASES"]
