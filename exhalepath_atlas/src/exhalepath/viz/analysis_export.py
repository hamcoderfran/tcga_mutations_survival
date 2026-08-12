"""Shared research-facing analysis exports (literature overlay, Methods, Prism CSV).

Builds on existing curated literature panels + fused/biomarker predictions so
reports are useful for papers without inventing new databases.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Optional


def _panel_for_disease(disease_id: str) -> dict[str, Any] | None:
    from ..data.literature_panels import literature_panel_for_disease

    panel = literature_panel_for_disease(disease_id)
    if panel:
        return panel
    # also try knowledge literature benchmarks
    alt = Path(__file__).resolve().parents[1] / "data" / "knowledge" / "literature_benchmarks.json"
    if alt.exists():
        doc = json.loads(alt.read_text())
        want = {
            disease_id.lower(),
            disease_id.lower().replace(" ", "_"),
        }
        for c in doc.get("cases") or []:
            key = (c.get("disease_id") or c.get("disease") or "").lower().replace(" ", "_")
            if key in want:
                return c
    return None


def literature_overlay(
    disease_id: str,
    predicted_log2fc: dict[str, float],
    *,
    prior_log2fc: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Compare predicted directions to published panel (if available).

    When ``prior_log2fc`` is supplied (or auto-loaded), also reports an
    outside-prior concordance that excludes VOCs already present in the atlas
    disease prior — those agreements are often curation-circular.
    """
    panel = _panel_for_disease(disease_id)
    if not panel:
        return {
            "available": False,
            "disease_id": disease_id,
            "rows": [],
            "refs": [],
            "n_agree": 0,
            "n_compared": 0,
        }

    if prior_log2fc is None:
        try:
            from ..knowledge.loader import default_knowledge

            prior_log2fc = dict(
                (default_knowledge().diseases.get(disease_id) or {}).get("voc_log2fc_prior")
                or {}
            )
        except Exception:
            prior_log2fc = {}

    measured = dict(panel.get("measured_log2fc") or {})
    # also support elevated/suppressed lists
    for v in panel.get("expect_elevated") or panel.get("elevated") or []:
        measured.setdefault(str(v), 1.0)
    for v in panel.get("expect_suppressed") or panel.get("suppressed") or []:
        measured.setdefault(str(v), -1.0)

    voc_evidence = dict(panel.get("voc_evidence") or {})
    rows = []
    agree = 0
    agree_out = 0
    compared_out = 0
    for voc, lit_fc in measured.items():
        ev = voc_evidence.get(voc) or {}
        grade = ev.get("evidence") or panel.get("evidence_grade") or "mixed"
        doi = ev.get("source_doi")
        if not doi and (panel.get("refs") or []):
            doi = (panel.get("refs") or [{}])[0].get("doi")
        in_prior = voc in (prior_log2fc or {})
        pred = predicted_log2fc.get(voc)
        if pred is None:
            rows.append(
                {
                    "voc_id": voc,
                    "literature_log2fc": float(lit_fc),
                    "predicted_log2fc": None,
                    "agree": None,
                    "status": "not_in_prediction",
                    "evidence_grade": grade,
                    "doi": doi,
                    "in_atlas_prior": in_prior,
                    "circularity_risk": True if in_prior else grade in {"directional_only", "mixed", "atlas_prior"},
                }
            )
            continue
        lit_sign = 1 if float(lit_fc) >= 0 else -1
        pred_sign = 1 if float(pred) >= 0 else -1
        ok = lit_sign == pred_sign
        agree += int(ok)
        if not in_prior:
            compared_out += 1
            agree_out += int(ok)
        rows.append(
            {
                "voc_id": voc,
                "literature_log2fc": float(lit_fc),
                "predicted_log2fc": float(pred),
                "agree": ok,
                "status": "agree" if ok else "disagree",
                "evidence_grade": grade,
                "doi": doi,
                "in_atlas_prior": in_prior,
                "circularity_risk": True if in_prior else grade in {"directional_only", "mixed", "atlas_prior"},
            }
        )
    compared = sum(1 for r in rows if r["agree"] is not None)
    refs = []
    for r in panel.get("refs") or []:
        refs.append(
            {
                "doi": r.get("doi"),
                "title": r.get("title"),
                "year": r.get("year"),
                "method": r.get("method"),
            }
        )
    return {
        "available": True,
        "disease_id": panel.get("disease_id") or disease_id,
        "disease_name": panel.get("disease_name") or panel.get("disease"),
        "rows": rows,
        "refs": refs,
        "n_agree": agree,
        "n_compared": compared,
        "directional_accuracy": (agree / compared) if compared else None,
        "n_agree_outside_prior": agree_out,
        "n_compared_outside_prior": compared_out,
        "directional_accuracy_outside_prior": (agree_out / compared_out) if compared_out else None,
        "n_panel_vocs_also_in_prior": sum(1 for r in rows if r.get("in_atlas_prior")),
        "note": (
            "directional_accuracy includes VOCs that also sit in atlas priors (often "
            "curation-circular). Prefer directional_accuracy_outside_prior when non-null; "
            "if that is null, every panel VOC was already in the prior."
        ),
    }


def write_literature_overlay_csv(overlay: dict[str, Any], path: Path) -> Path:
    path = Path(path)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "voc_id",
                "literature_log2fc",
                "predicted_log2fc",
                "agree",
                "status",
                "evidence_grade",
                "doi",
                "in_atlas_prior",
                "circularity_risk",
            ],
        )
        w.writeheader()
        for row in overlay.get("rows") or []:
            w.writerow(row)
    return path


def write_graphpad_long_csv(
    rows: Iterable[dict[str, Any]],
    path: Path,
) -> Path:
    """Prism/GraphPad-friendly long table: group, voc, metric, estimate, ci_low, ci_high, source."""
    path = Path(path)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "group",
                "voc",
                "metric",
                "estimate",
                "ci_low",
                "ci_high",
                "source",
                "confidence",
            ],
        )
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k) for k in w.fieldnames})
    return path


def methods_markdown(
    *,
    tool: str,
    disease: str,
    location: Optional[str],
    model_version: str,
    extra_bullets: list[str] | None = None,
) -> str:
    bullets = [
        f"Predictions generated with **{tool}** (`{model_version}`).",
        f"Query disease: **{disease}**"
        + (f"; location hint: **{location}**." if location else "."),
        "Mechanism stack combines physiology (Farhi / cell-state), pathway priors, "
        "optional ML calibrator residuals, and multi-model fusion when using `voc stack`.",
        "Exhaled VOC outputs are hypothesized Δppb / log2 fold-changes vs a healthy reference — "
        "not measured clinical concentrations.",
        "Research / hypothesis-generation only; not a medical device.",
    ]
    bullets.extend(extra_bullets or [])
    lines = [
        "# Methods (paste-ready draft)",
        "",
        "Auto-generated from the run configuration. Edit before publication.",
        "",
    ]
    for b in bullets:
        lines.append(f"- {b}")
    lines += [
        "",
        "## Suggested citation framing",
        "",
        "We used ExhalePath Atlas / voc-breath, an open mechanism-aware exhaled-VOC "
        "prediction stack, to generate ranked breath VOC hypotheses and uncertainty "
        "estimates for the queried disease context.",
        "",
    ]
    return "\n".join(lines)


def next_experiments(
    *,
    top_vocs: list[str],
    overlay: dict[str, Any] | None = None,
    mean_epistemic: float | None = None,
    n_models_ok: int | None = None,
    n_models_total: int | None = None,
) -> list[str]:
    tips: list[str] = []
    if top_vocs:
        tips.append(
            "Validate top hypothesized VOCs with targeted GC-MS/PTR "
            f"({', '.join(top_vocs[:5])}) against matched healthy controls."
        )
    if overlay and overlay.get("available"):
        disag = [r["voc_id"] for r in overlay.get("rows") or [] if r.get("agree") is False]
        missing = [r["voc_id"] for r in overlay.get("rows") or [] if r.get("status") == "not_in_prediction"]
        if disag:
            tips.append(
                "Resolve literature disagreements with orthogonal assay / blank-aware "
                f"quantitation: {', '.join(disag[:5])}."
            )
        if missing:
            tips.append(
                "Expand peak→atlas VOC mapping for literature compounds not in the "
                f"prediction set: {', '.join(missing[:5])}."
            )
        acc = overlay.get("directional_accuracy")
        if acc is not None and acc < 0.7:
            tips.append(
                "Directional agreement with published panels is modest — prioritize "
                "patient-level holdout (voc eval-patient-diagnostic) over panel reuse."
            )
    if mean_epistemic is not None and mean_epistemic > 0.4:
        tips.append(
            "High cross-model epistemic disagreement — collect smoking/age/site metadata "
            "and re-score with confounder-aware splits."
        )
    if n_models_ok is not None and n_models_total is not None and n_models_ok < n_models_total:
        tips.append(
            f"Only {n_models_ok}/{n_models_total} stack heads fired — supply genes, "
            "comorbidities, or description text to activate more channels."
        )
    tips.append(
        "Lock a patient-level split (`voc lock-split`) before claiming diagnostic "
        "performance; report nested AUROC and optimism gap."
    )
    return tips[:6]
