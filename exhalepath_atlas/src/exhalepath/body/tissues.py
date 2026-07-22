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
        # Punctuation-only / garbage queries normalize to "" and must not
        # fuzzy-match every alias via empty-string containment ("" in alias).
        if not key:
            return {
                "tissue_id": f"custom_unmatched_{abs(hash(str(location))) % 10_000_000}",
                "name": str(location).strip(),
                "aliases": [str(location).strip()],
                "n_census_healthy_cells": 0,
                "compartment": "custom",
                "matched": False,
                "query": location,
            }
        # Underscore / hyphen variants of the same anatomic token
        key_compact = key.replace(" ", "")
        if key in self._alias:
            tid = self._alias[key]
            t = self.tissues[tid]
            return {**t, "matched": True, "query": location}
        for alias, tid in self._alias.items():
            if alias.replace(" ", "") == key_compact:
                t = self.tissues[tid]
                return {**t, "matched": True, "query": location}
        # Fuzzy: prefer longest whole-token alias containment so
        # "left breast upper outer" → breast, while "head neck" prefers
        # head_neck over the shorter "head" alias.
        key_tokens = [t for t in key.split() if t]
        key_token_set = set(key_tokens)
        candidates: list[tuple[float, int, str]] = []
        for alias, tid in self._alias.items():
            if not alias:
                continue
            alias_tokens = [t for t in alias.split() if t]
            alias_token_set = set(alias_tokens)
            if not alias_token_set:
                continue
            subset = alias_token_set <= key_token_set
            supersets = key_token_set <= alias_token_set
            contained = key in alias or alias in key
            if not (subset or supersets or contained):
                continue
            # Reject weak short substring-only hits that are not token-aligned
            # (keeps "head" from beating "head neck", but allows long clinical
            # stems like "peritoneum" inside "retroperitoneal …").
            if contained and not (subset or supersets):
                coverage = min(len(alias), len(key)) / max(len(alias), len(key))
                long_stem = len(alias) >= 8 and alias in key
                if coverage < 0.55 and not long_stem:
                    continue
            n_overlap = len(alias_token_set & key_token_set)
            score = (
                n_overlap * 100.0
                + len(alias)
                + (50.0 if subset else 0.0)
                + (30.0 if (contained and len(alias) >= 8) else 0.0)
            )
            candidates.append((score, len(alias), tid))
        if candidates:
            candidates.sort(reverse=True)
            tid = candidates[0][2]
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
