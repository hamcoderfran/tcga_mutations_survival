from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from ..config import CACHE_DIR, REACTOME_CONTENT
from .http import request_json


class ReactomeClient:
    """Expand seed pathways with Reactome Content Service gene sets when online."""

    def __init__(self, cache_dir: Path | None = None):
        self.cache_dir = Path(cache_dir or CACHE_DIR / "reactome")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def participants(self, reactome_id: str) -> list[str]:
        cache = self.cache_dir / f"{reactome_id}_genes.txt"
        if cache.exists():
            return [g.strip() for g in cache.read_text().splitlines() if g.strip()]

        url = f"{REACTOME_CONTENT}/data/participants/{reactome_id}/participatingPhysicalEntities"
        try:
            data = request_json("GET", url, timeout=60)
        except RuntimeError:
            return []

        genes: set[str] = set()
        if isinstance(data, list):
            for ent in data:
                self._collect_genes(ent, genes)
        cache.write_text("\n".join(sorted(genes)))
        return sorted(genes)

    def _collect_genes(self, ent: dict[str, Any], genes: set[str]) -> None:
        display = (ent.get("displayName") or "").upper()
        # UniProt-style gene symbols often appear as "GENE gene"
        for token in display.replace(",", " ").split():
            if token.isalpha() and 2 <= len(token) <= 15 and token == token.upper():
                genes.add(token)
        for child_key in ("hasMember", "hasComponent", "repeatedUnit"):
            for child in ent.get(child_key) or []:
                if isinstance(child, dict):
                    self._collect_genes(child, genes)

    def expand_pathway_map(self, pathway_map: dict[str, dict]) -> pd.DataFrame:
        rows = []
        for pid, p in pathway_map.items():
            seed = {g.upper() for g in p.get("seed_genes", [])}
            expanded = set(seed)
            for rid in p.get("reactome_ids", []):
                expanded |= {g.upper() for g in self.participants(rid)}
            for g in sorted(expanded):
                rows.append(
                    {
                        "pathway_id": pid,
                        "pathway_name": p.get("name"),
                        "gene_symbol": g,
                        "is_seed": g in seed,
                    }
                )
        return pd.DataFrame(rows)
