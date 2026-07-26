"""Ontology / alias nearest-neighbor transfer for unresolved diseases.

Transfers category, site, pathway_bias, driver genes, and cell-state donor IDs.
Never copies voc_log2fc_prior unless explicitly allowed.
"""

from __future__ import annotations

import re
from typing import Any

from .category_templates import merge_pathway_bias, template_for_category


def _norm(s: str) -> str:
    s = str(s or "").strip().lower()
    s = re.sub(r"^tcga-", "", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _tokens(s: str) -> set[str]:
    stop = {
        "cancer",
        "carcinoma",
        "disease",
        "syndrome",
        "disorder",
        "the",
        "of",
        "and",
        "type",
        "primary",
        "chronic",
        "acute",
        "benign",
        "essential",
    }
    return {t for t in _norm(s).split() if t and t not in stop and len(t) > 2}


def nearest_atlas_diseases(
    query: str,
    diseases: dict[str, dict[str, Any]],
    *,
    top_k: int = 3,
    min_score: float = 0.28,
) -> list[tuple[float, dict[str, Any]]]:
    """Return ranked (score, disease_dict) neighbors by token / alias / MONDO overlap."""
    q = _norm(query)
    q_tokens = _tokens(query)
    if not q and not q_tokens:
        return []

    scored: list[tuple[float, dict[str, Any]]] = []
    for d in diseases.values():
        keys = [d.get("disease_id", ""), d.get("name", ""), d.get("mondo_id") or ""]
        keys.extend(list(d.get("aliases") or []))
        best = 0.0
        for raw in keys:
            if not raw:
                continue
            a = _norm(str(raw))
            if not a:
                continue
            if q == a:
                best = max(best, 1.0)
                continue
            # Avoid gene-symbol / short-token traps (BRCA1 ⊂ BRCA alias, etc.)
            if len(q) <= 6 and q.isalpha() and any(ch.isdigit() for ch in query):
                continue
            if q in a or a in q:
                # Require substantial coverage so "cancer" does not NN-match every *cancer*
                cov = min(len(q), len(a)) / max(len(q), len(a), 1)
                if cov >= 0.55:
                    best = max(best, 0.55 * cov)
            a_tokens = _tokens(a)
            if q_tokens and a_tokens:
                overlap = len(q_tokens & a_tokens) / len(q_tokens | a_tokens)
                # Prefer neighbors that cover most query tokens
                coverage = len(q_tokens & a_tokens) / max(len(q_tokens), 1)
                best = max(best, 0.65 * overlap + 0.35 * coverage)
            # MONDO id exact
            if str(raw).upper().startswith("MONDO:") and q.replace(" ", "") == a.replace(" ", ""):
                best = max(best, 1.0)
        # Category token boost when query mentions neighbor category
        cat = (d.get("category") or "").lower()
        if cat and cat in q:
            best = max(best, 0.3)
        if best >= min_score:
            scored.append((best, d))

    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[:top_k]


def transfer_from_neighbor(
    neighbor: dict[str, Any],
    *,
    score: float,
    copy_voc_priors: bool = False,
    pathway_scale: float = 0.85,
) -> dict[str, Any]:
    """Build a mechanistic transfer payload from an atlas neighbor."""
    bias = dict(neighbor.get("pathway_bias") or {})
    if not bias:
        bias = template_for_category(neighbor.get("category"))
    # Soften transferred bias so neighbors don't fully impersonate curated diseases
    scaled = {
        pid: 1.0 + pathway_scale * score * (float(v) - 1.0)
        for pid, v in bias.items()
        if float(v) > 1.0
    }
    out: dict[str, Any] = {
        "category": neighbor.get("category") or "default",
        "default_site": neighbor.get("default_site"),
        "pathway_bias": scaled,
        "driver_genes": list(neighbor.get("driver_genes") or [])[:12],
        "mondo_id": neighbor.get("mondo_id"),
        "neighbor_disease_id": neighbor.get("disease_id"),
        "neighbor_name": neighbor.get("name"),
        "neighbor_score": float(score),
        "cell_state_donor_ids": [neighbor.get("disease_id")]
        if neighbor.get("disease_id")
        else [],
        "voc_log2fc_prior": dict(neighbor.get("voc_log2fc_prior") or {})
        if copy_voc_priors
        else {},
    }
    return out


def ontology_transfer(
    query: str,
    diseases: dict[str, dict[str, Any]],
    *,
    copy_voc_priors: bool = False,
) -> dict[str, Any] | None:
    neighbors = nearest_atlas_diseases(query, diseases, top_k=2, min_score=0.34)
    if not neighbors:
        return None
    score, neigh = neighbors[0]
    # Require a moderately strong match before transferring
    if score < 0.34:
        return None
    payload = transfer_from_neighbor(
        neigh, score=score, copy_voc_priors=copy_voc_priors
    )
    # Merge soft pathway hints from #2 if same category
    if len(neighbors) > 1 and neighbors[1][1].get("category") == neigh.get("category"):
        extra = transfer_from_neighbor(
            neighbors[1][1], score=neighbors[1][0] * 0.7, copy_voc_priors=False
        )
        payload["pathway_bias"] = merge_pathway_bias(
            payload.get("pathway_bias"), extra.get("pathway_bias")
        )
        donors = list(payload.get("cell_state_donor_ids") or [])
        donors.extend(extra.get("cell_state_donor_ids") or [])
        payload["cell_state_donor_ids"] = list(dict.fromkeys(donors))
    payload["neighbors"] = [
        {
            "disease_id": d.get("disease_id"),
            "name": d.get("name"),
            "score": float(s),
            "category": d.get("category"),
        }
        for s, d in neighbors
    ]
    return payload
