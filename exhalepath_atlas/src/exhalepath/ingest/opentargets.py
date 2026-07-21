from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from ..config import CACHE_DIR, OPENTARGETS_API
from .http import request_json

_DISEASE_SEARCH = """
query searchDisease($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
    hits {
      id
      name
      entity
    }
  }
}
"""

_ASSOC = """
query diseaseAssoc($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: 0, size: 100}) {
      count
      rows {
        score
        target {
          approvedSymbol
          id
        }
      }
    }
  }
}
"""


class OpenTargetsClient:
    """Disease→target associations for non-cancer / pan-disease queries."""

    def __init__(self, cache_dir: Path | None = None):
        self.cache_dir = Path(cache_dir or CACHE_DIR / "opentargets")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _gql(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        return request_json(
            "POST",
            OPENTARGETS_API,
            json_body={"query": query, "variables": variables},
            timeout=60,
        )

    def disease_targets(self, disease_query: str, min_score: float = 0.2) -> pd.DataFrame:
        safe = "".join(c if c.isalnum() else "_" for c in disease_query.lower())[:80]
        cache = self.cache_dir / f"{safe}_targets.csv"
        if cache.exists():
            return pd.read_csv(cache)

        try:
            search = self._gql(_DISEASE_SEARCH, {"q": disease_query})
            hits = (((search.get("data") or {}).get("search") or {}).get("hits")) or []
            if not hits:
                return pd.DataFrame(columns=["disease_id", "disease_name", "gene_symbol", "score"])
            efo = hits[0]["id"]
            assoc = self._gql(_ASSOC, {"efoId": efo})
            disease = ((assoc.get("data") or {}).get("disease")) or {}
            rows = []
            for r in ((disease.get("associatedTargets") or {}).get("rows")) or []:
                score = float(r.get("score") or 0)
                if score < min_score:
                    continue
                sym = ((r.get("target") or {}).get("approvedSymbol")) or ""
                if not sym:
                    continue
                rows.append(
                    {
                        "disease_id": disease.get("id"),
                        "disease_name": disease.get("name"),
                        "gene_symbol": sym.upper(),
                        "score": score,
                    }
                )
            df = pd.DataFrame(rows)
            df.to_csv(cache, index=False)
            return df
        except RuntimeError:
            return pd.DataFrame(columns=["disease_id", "disease_name", "gene_symbol", "score"])
