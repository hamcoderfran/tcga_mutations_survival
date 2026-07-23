from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..config import KNOWLEDGE_DIR


def _norm(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


class WholeBodyMap:
    """Entire-human-body anatomic location index for affected-cell placement."""

    def __init__(self, path: Path | None = None):
        p = Path(path or KNOWLEDGE_DIR / "whole_body_tissues.json")
        doc = json.loads(p.read_text())
        self.tissues = {t["tissue_id"]: t for t in doc["tissues"]}
        self._alias = {}
        for tid, t in self.tissues.items():
            keys = [tid, t.get("name", "")] + list(t.get("aliases") or [])
            for k in keys:
                if k:
                    self._alias[_norm(k)] = tid

    def resolve(self, location: str | None) -> dict[str, Any]:
        if not location or not str(location).strip():
            return {
                "tissue_id": "systemic",
                "name": "Systemic / whole body",
                "matched": False,
                "query": location,
            }
        key = _norm(str(location))
        if key in self._alias:
            tid = self._alias[key]
            t = self.tissues[tid]
            return {**t, "matched": True, "query": location}
        # fuzzy containment
        candidates = []
        for alias, tid in self._alias.items():
            if key in alias or alias in key:
                candidates.append((len(alias), tid))
        if candidates:
            candidates.sort(reverse=True)
            tid = candidates[0][1]
            t = self.tissues[tid]
            return {**t, "matched": True, "query": location}
        # free-text location still usable as primary_site string
        return {
            "tissue_id": key.replace(" ", "_"),
            "name": str(location),
            "aliases": [str(location)],
            "n_census_healthy_cells": 0,
            "compartment": "custom",
            "matched": False,
            "query": location,
        }

    def list_locations(self) -> list[dict[str, Any]]:
        return sorted(self.tissues.values(), key=lambda t: (-(t.get("n_census_healthy_cells") or 0), t["name"]))


@lru_cache(maxsize=1)
def default_body_map() -> WholeBodyMap:
    return WholeBodyMap()


def resolve_location(location: str | None) -> dict[str, Any]:
    return default_body_map().resolve(location)
