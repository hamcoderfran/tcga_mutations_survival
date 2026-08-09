"""Evidence-graded VOC↔disease claim ledger (citeable research graph).

Unifies priority literature panels, atlas ``voc_log2fc_prior``, and (when
present) mechanism-pack WHY notes into one ledger with grades:

- quantified — numeric FC from a cited human breath study
- directional_only — up/down without a published number
- atlas_prior — mechanism/template prior (may be circular with overlays)
- mixed — panel mixes quantified + directional
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from ..config import DATA_DIR, PACKAGE_ROOT


def _load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def _priority_panels() -> list[dict[str, Any]]:
    path = DATA_DIR / "real_breath" / "literature_panels" / "priority10_voc_panels.json"
    if not path.exists():
        alt = PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "priority10_voc_panels.json"
        path = alt if alt.exists() else path
    if not path.exists():
        return []
    return list((_load_json(path).get("panels") or []))


def _disease_priors() -> list[dict[str, Any]]:
    path = DATA_DIR / "knowledge" / "disease_voc_priors.json"
    if not path.exists():
        alt = PACKAGE_ROOT / "data" / "knowledge" / "disease_voc_priors.json"
        path = alt if alt.exists() else path
    if not path.exists():
        return []
    return list((_load_json(path).get("diseases") or []))


def _mechanism_pack_vocs(disease_id: str) -> dict[str, Any]:
    candidates = [
        DATA_DIR / "knowledge" / "disease_mechanism_packs.json",
        PACKAGE_ROOT / "data" / "knowledge" / "disease_mechanism_packs.json",
    ]
    for path in candidates:
        if not path.exists():
            continue
        payload = _load_json(path)
        packs = payload.get("packs") or payload.get("diseases") or payload
        if isinstance(packs, dict):
            return packs.get(disease_id) or {}
        if isinstance(packs, list):
            for p in packs:
                if p.get("disease_id") == disease_id:
                    return p
    return {}


def build_claim_ledger(
    *,
    disease_ids: list[str] | None = None,
    include_atlas_priors: bool = True,
) -> dict[str, Any]:
    """Build the unified claim ledger."""
    panels = {p["disease_id"]: p for p in _priority_panels() if p.get("disease_id")}
    priors = {d["disease_id"]: d for d in _disease_priors() if d.get("disease_id")}
    ids = sorted(set(disease_ids or (set(panels) | set(priors))))

    claims: list[dict[str, Any]] = []
    for did in ids:
        panel = panels.get(did) or {}
        prior = priors.get(did) or {}
        voc_ev = panel.get("voc_evidence") or {}
        measured = panel.get("measured_log2fc") or {}
        refs = panel.get("refs") or []
        primary_doi = None
        if refs:
            primary_doi = refs[0].get("doi")

        # Literature panel claims
        for voc, fc in measured.items():
            ev = voc_ev.get(voc) or {}
            grade = ev.get("evidence") or panel.get("evidence_grade") or "mixed"
            claims.append(
                {
                    "claim_id": f"lit:{did}:{voc}",
                    "disease_id": did,
                    "voc_id": voc,
                    "log2fc": float(fc),
                    "direction": "increased"
                    if float(fc) > 0
                    else ("decreased" if float(fc) < 0 else "unchanged"),
                    "evidence_grade": grade,
                    "source": "literature_panel",
                    "doi": ev.get("source_doi") or primary_doi,
                    "circularity_risk": grade in {"directional_only", "mixed"},
                    "notes": ev.get("calculation"),
                }
            )

        # Atlas priors for VOCs not already covered
        if include_atlas_priors:
            covered = set(measured)
            for voc, fc in (prior.get("voc_log2fc_prior") or {}).items():
                if voc in covered:
                    continue
                claims.append(
                    {
                        "claim_id": f"prior:{did}:{voc}",
                        "disease_id": did,
                        "voc_id": voc,
                        "log2fc": float(fc),
                        "direction": "increased"
                        if float(fc) > 0
                        else ("decreased" if float(fc) < 0 else "unchanged"),
                        "evidence_grade": "atlas_prior",
                        "source": "disease_voc_priors",
                        "doi": None,
                        "circularity_risk": True,
                        "notes": "Mechanism/template prior — not an independent clinical measurement",
                    }
                )

        pack = _mechanism_pack_vocs(did)
        why = pack.get("voc_why") or pack.get("why") or {}
        if isinstance(why, dict):
            for voc, note in why.items():
                claims.append(
                    {
                        "claim_id": f"mech:{did}:{voc}",
                        "disease_id": did,
                        "voc_id": voc,
                        "log2fc": None,
                        "direction": None,
                        "evidence_grade": "mechanism_narrative",
                        "source": "disease_mechanism_packs",
                        "doi": None,
                        "circularity_risk": False,
                        "notes": str(note)[:500],
                    }
                )

    # summary
    by_grade: dict[str, int] = {}
    for c in claims:
        g = str(c.get("evidence_grade") or "unknown")
        by_grade[g] = by_grade.get(g, 0) + 1

    return {
        "schema_version": "ClaimLedger-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_diseases": len(ids),
        "n_claims": len(claims),
        "by_evidence_grade": by_grade,
        "diseases": ids,
        "claims": claims,
        "honesty": (
            "Literature quantified claims are the strongest open evidence here. "
            "Atlas priors are hypothesis-generating and may overlap literature overlays. "
            "Not a clinical diagnostic knowledge base."
        ),
    }


def export_claim_ledger(out_path: Path, **kwargs: Any) -> Path:
    ledger = build_claim_ledger(**kwargs)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(ledger, indent=2))
    md = out_path.with_suffix(".md")
    lines = [
        "# VOC↔disease claim ledger",
        "",
        f"Generated: {ledger['generated_utc']}",
        f"Claims: **{ledger['n_claims']}** across **{ledger['n_diseases']}** diseases",
        "",
        "## By evidence grade",
        "",
    ]
    for g, n in sorted((ledger.get("by_evidence_grade") or {}).items(), key=lambda x: -x[1]):
        lines.append(f"- `{g}`: {n}")
    lines += ["", f"> {ledger.get('honesty')}", ""]
    # top quantified
    quant = [c for c in ledger["claims"] if c.get("evidence_grade") == "quantified"][:15]
    if quant:
        lines += ["## Sample quantified claims", ""]
        for c in quant:
            lines.append(
                f"- `{c['disease_id']}` / `{c['voc_id']}` log2fc={c.get('log2fc')} "
                f"doi:{c.get('doi')}"
            )
        lines.append("")
    md.write_text("\n".join(lines))
    return out_path


__all__ = ["build_claim_ledger", "export_claim_ledger"]
