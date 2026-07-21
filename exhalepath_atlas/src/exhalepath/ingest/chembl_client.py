"""ChEMBL web-services client (EBI).

ChEMBL 36+ holds ~2.9M compounds and ~24M bioactivities. ExhalePath uses a
pathway-gene slice of that graph — not as exhaled-ppb labels, but as
chemogenomic supervision for enzyme/pathway modulation and VOC physicochemical
priors.
"""

from __future__ import annotations

from typing import Any, Iterator

from ..config import CHEMBL_API, USER_AGENT
from .http import request_json

# Prefer binding / functional potency endpoints with pChEMBL
DEFAULT_STANDARD_TYPES = ("IC50", "Ki", "Kd", "EC50", "AC50", "Potency")


class ChEMBLClient:
    def __init__(self, base_url: str | None = None):
        self.base = (base_url or CHEMBL_API).rstrip("/")

    def status(self) -> dict[str, Any]:
        return request_json("GET", f"{self.base}/status.json")

    def find_targets_for_gene(self, gene: str, *, limit: int = 20) -> list[dict[str, Any]]:
        gene = gene.strip().upper()
        if not gene:
            return []
        data = request_json(
            "GET",
            f"{self.base}/target.json",
            params={"target_synonym__icontains": gene, "limit": limit},
        )
        targets = list(data.get("targets") or [])
        # Prefer human single-protein targets whose component synonyms/accessions match
        scored: list[tuple[int, dict[str, Any]]] = []
        for t in targets:
            score = 0
            if t.get("target_type") == "SINGLE PROTEIN":
                score += 5
            org = (t.get("organism") or "").lower()
            if "homo sapiens" in org or org == "homo sapiens":
                score += 5
            pref = (t.get("pref_name") or "").upper()
            if gene in pref.replace(" ", ""):
                score += 2
            for syn in t.get("target_components") or []:
                for s in syn.get("target_component_synonyms") or []:
                    if str(s.get("component_synonym", "")).upper() == gene:
                        score += 8
            scored.append((score, t))
        scored.sort(key=lambda x: x[0], reverse=True)
        return [t for _, t in scored]

    def iter_activities(
        self,
        *,
        target_chembl_id: str,
        standard_types: tuple[str, ...] = DEFAULT_STANDARD_TYPES,
        page_size: int = 1000,
        max_rows: int | None = None,
        require_pchembl: bool = True,
    ) -> Iterator[dict[str, Any]]:
        """Yield activity records for a target (paginated)."""
        fetched = 0
        for st in standard_types:
            offset = 0
            while True:
                if max_rows is not None and fetched >= max_rows:
                    return
                limit = page_size
                if max_rows is not None:
                    limit = min(page_size, max_rows - fetched)
                params: dict[str, Any] = {
                    "target_chembl_id": target_chembl_id,
                    "standard_type": st,
                    "limit": limit,
                    "offset": offset,
                }
                if require_pchembl:
                    params["pchembl_value__isnull"] = "False"
                data = request_json(
                    "GET", f"{self.base}/activity.json", params=params, timeout=120
                )
                acts = list(data.get("activities") or [])
                if not acts:
                    break
                for a in acts:
                    yield a
                    fetched += 1
                    if max_rows is not None and fetched >= max_rows:
                        return
                meta = data.get("page_meta") or {}
                total = int(meta.get("total_count") or 0)
                offset += len(acts)
                if offset >= total or meta.get("next") is None:
                    break

    def molecule(self, chembl_id: str) -> dict[str, Any]:
        return request_json("GET", f"{self.base}/molecule/{chembl_id}.json")

    def search_molecule(self, query: str, *, limit: int = 5) -> list[dict[str, Any]]:
        data = request_json(
            "GET",
            f"{self.base}/molecule/search.json",
            params={"q": query, "limit": limit},
        )
        return list(data.get("molecules") or [])
