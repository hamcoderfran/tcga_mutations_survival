from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..config import KNOWLEDGE_DIR


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as f:
        return json.load(f)


class KnowledgeBase:
    """Curated VOC / pathway / disease priors bundled with ExhalePath."""

    def __init__(self, knowledge_dir: Path | None = None):
        self.root = Path(knowledge_dir or KNOWLEDGE_DIR)
        self.voc_catalog = _load_json(self.root / "voc_catalog.json")
        self.pathway_map = _load_json(self.root / "pathway_voc_map.json")
        self.disease_priors = _load_json(self.root / "disease_voc_priors.json")

        self.vocs = {v["voc_id"]: v for v in self.voc_catalog["vocs"]}
        self.pathways = {p["pathway_id"]: p for p in self.pathway_map["pathways"]}
        self.diseases = {d["disease_id"]: d for d in self.disease_priors["diseases"]}
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

    def resolve_disease(self, query: str) -> dict[str, Any]:
        key = self._norm(query)
        if key in self._alias_index:
            return self.diseases[self._alias_index[key]]

        # Fuzzy substring match on aliases/names
        for alias, did in self._alias_index.items():
            if key in alias or alias in key:
                return self.diseases[did]

        # Generic fallback: treat as unknown disease with neutral priors
        return {
            "disease_id": f"custom::{self._norm(query).replace(' ', '_')}",
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


@lru_cache(maxsize=1)
def default_knowledge() -> KnowledgeBase:
    return KnowledgeBase()
