from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..config import KNOWLEDGE_DIR

# Aliases shorter than this must match exactly (prevents "cancer" → first *cancer* hit)
_MIN_FUZZY_ALIAS_LEN = 4


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as f:
        return json.load(f)


class KnowledgeBase:
    """Curated VOC / pathway / disease / physiology priors bundled with ExhalePath."""

    def __init__(self, knowledge_dir: Path | None = None):
        self.root = Path(knowledge_dir or KNOWLEDGE_DIR)
        self.voc_catalog = _load_json(self.root / "voc_catalog.json")
        self.pathway_map = _load_json(self.root / "pathway_voc_map.json")
        self.disease_priors = _load_json(self.root / "disease_voc_priors.json")
        self.pathway_chains_doc = _load_json(self.root / "voc_pathway_chains.json")
        self.cell_state_doc = _load_json(self.root / "cell_state_atlas.json")
        self.physio_constants = _load_json(self.root / "physio_constants.json")

        self.vocs = {v["voc_id"]: v for v in self.voc_catalog["vocs"]}
        self.pathways = {p["pathway_id"]: p for p in self.pathway_map["pathways"]}
        self.diseases = {d["disease_id"]: d for d in self.disease_priors["diseases"]}
        self.pathway_chains = list(self.pathway_chains_doc.get("chains") or [])
        self.cell_states = list(self.cell_state_doc.get("cell_states") or [])
        self._alias_index = self._build_alias_index()

    def _build_alias_index(self) -> dict[str, str]:
        idx: dict[str, str] = {}
        for did, d in self.diseases.items():
            keys = [did, d["name"], d.get("mondo_id", "")] + list(d.get("aliases", []))
            for k in keys:
                if not k:
                    continue
                idx[self._norm(k)] = did
        return idx

    @staticmethod
    def _norm(s: str) -> str:
        s = s.strip().lower()
        s = re.sub(r"^tcga-", "", s)
        s = re.sub(r"[^a-z0-9]+", " ", s)
        return re.sub(r"\s+", " ", s).strip()

    @staticmethod
    def _token_set(s: str) -> set[str]:
        stop = {"cancer", "carcinoma", "disease", "syndrome", "the", "of", "and"}
        return {t for t in s.split() if t and t not in stop}

    def resolve_disease(self, query: str) -> dict[str, Any]:
        if query is None or not str(query).strip():
            return self._unresolved("unknown")

        key = self._norm(str(query))
        if not key:
            return self._unresolved(str(query))

        if key in self._alias_index:
            return dict(self.diseases[self._alias_index[key]])

        # Prefer longest alias containment / token overlap over first substring hit
        q_tokens = self._token_set(key)
        candidates: list[tuple[float, int, str]] = []
        for alias, did in self._alias_index.items():
            if len(alias) < _MIN_FUZZY_ALIAS_LEN:
                continue
            score = 0.0
            if key == alias:
                score = 100.0
            elif key in alias:
                score = 40.0 * (len(key) / max(len(alias), 1))
            elif alias in key:
                score = 55.0 * (len(alias) / max(len(key), 1))
            else:
                a_tokens = self._token_set(alias)
                if not a_tokens or not q_tokens:
                    continue
                overlap = len(a_tokens & q_tokens) / len(a_tokens | q_tokens)
                if overlap < 0.5:
                    continue
                score = 50.0 * overlap
            if score > 0:
                candidates.append((score, len(alias), did))

        if candidates:
            candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
            best_score, _, did = candidates[0]
            if best_score >= 35.0:
                return dict(self.diseases[did])

        return self._unresolved(str(query))

    def _unresolved(self, query: str) -> dict[str, Any]:
        return {
            "disease_id": f"custom::{self._norm(query).replace(' ', '_') or 'unknown'}",
            "name": query,
            "aliases": [query],
            "mondo_id": None,
            "category": "unspecified",
            "default_site": None,
            "pathway_bias": {},
            "voc_log2fc_prior": {},
            "_unresolved": True,
        }

    def pathway_gene_universe(self) -> dict[str, set[str]]:
        return {
            pid: {g.upper() for g in p.get("seed_genes", [])}
            for pid, p in self.pathways.items()
        }

    def all_seed_genes(self) -> set[str]:
        genes: set[str] = set()
        for geneset in self.pathway_gene_universe().values():
            genes |= geneset
        return genes

    def chains_for_voc(self, voc_id: str) -> list[dict[str, Any]]:
        return [c for c in self.pathway_chains if c.get("voc_id") == voc_id]


@lru_cache(maxsize=1)
def default_knowledge() -> KnowledgeBase:
    return KnowledgeBase()


def clear_knowledge_cache() -> None:
    """Drop cached KnowledgeBase (use after swapping knowledge JSON in tests)."""
    default_knowledge.cache_clear()
