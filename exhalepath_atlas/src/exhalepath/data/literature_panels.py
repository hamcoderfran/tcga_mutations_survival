"""Load curated literature VOC panels from ``data/real_breath/literature_panels``.

Merges every ``*_voc_panels.json`` (priority10, mental_health, …) so claim ledger,
overlays, and lit_compare share one discovery path.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..config import DATA_DIR, PACKAGE_ROOT


def literature_panel_dirs() -> list[Path]:
    return [
        DATA_DIR / "real_breath" / "literature_panels",
        PACKAGE_ROOT / "data" / "real_breath" / "literature_panels",
        Path("data/real_breath/literature_panels"),
    ]


def load_literature_panels() -> list[dict[str, Any]]:
    """Return merged panel objects keyed later by disease_id (last file wins on clash)."""
    by_id: dict[str, dict[str, Any]] = {}
    seen_dirs: set[Path] = set()
    for d in literature_panel_dirs():
        try:
            resolved = d.resolve()
        except OSError:
            continue
        if not d.is_dir() or resolved in seen_dirs:
            continue
        seen_dirs.add(resolved)
        for path in sorted(d.glob("*_voc_panels.json")):
            try:
                doc = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError):
                continue
            for panel in doc.get("panels") or []:
                did = panel.get("disease_id")
                if did:
                    by_id[str(did)] = panel
    return list(by_id.values())


def literature_panel_for_disease(disease_id: str) -> dict[str, Any] | None:
    did = disease_id.lower().replace(" ", "_")
    for panel in load_literature_panels():
        pid = str(panel.get("disease_id") or "").lower()
        if pid == did or did in pid or pid in did:
            return panel
    return None


def load_thin_evidence_notes() -> dict[str, Any]:
    """Optional honesty notes for conditions without measured breath panels."""
    for d in literature_panel_dirs():
        path = d / "mental_health_voc_panels.json"
        if not path.exists():
            continue
        try:
            doc = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        notes = doc.get("thin_evidence_conditions") or {}
        if notes:
            return dict(notes)
    return {}


__all__ = [
    "literature_panel_dirs",
    "load_literature_panels",
    "literature_panel_for_disease",
    "load_thin_evidence_notes",
]
