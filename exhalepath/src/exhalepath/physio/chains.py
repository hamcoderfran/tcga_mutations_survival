from __future__ import annotations

import math
from typing import Any

from ..schemas import PathwayScore


def score_pathway_chain(
    chain: dict[str, Any],
    pathway_scores: dict[str, float],
    mutated_genes: set[str],
) -> dict[str, Any]:
    """
    Score a multi-step VOC biosynthetic chain.

    Flux ≈ geometric-mean of step supports (pathway activity + enzyme hits),
    with detox steps subtracting capacity for net VOC (e.g. urea cycle).
    """
    step_scores: list[float] = []
    enzyme_hits: list[str] = []
    pathway_drivers: list[str] = []
    detox_capacity = 0.0

    for step in chain.get("steps") or []:
        pids = list(step.get("pathway_ids") or [])
        pw = max((pathway_scores.get(p, 0.0) for p in pids), default=0.0)
        enzymes = [e.upper() for e in step.get("enzymes") or []]
        e_hits = [e for e in enzymes if e in mutated_genes]
        enzyme_hits.extend(e_hits)
        enzyme_frac = (len(e_hits) / len(enzymes)) if enzymes else (0.35 if step.get("spontaneous") else 0.0)
        support = max(pw, 0.0) * 0.7 + enzyme_frac * 0.5 + (0.15 if step.get("spontaneous") else 0.0)
        if step.get("detox"):
            detox_capacity = max(detox_capacity, support)
            # Detox reduces net VOC; represent as inverse pressure later
            support = max(0.05, 0.4 - 0.3 * support)
        step_scores.append(max(support, 1e-6))
        for p in pids:
            if pathway_scores.get(p, 0.0) > 0:
                pathway_drivers.append(p)

    # Bottleneck-aware flux: geometric mean
    if not step_scores:
        flux = 0.0
    else:
        flux = math.exp(sum(math.log(s) for s in step_scores) / len(step_scores))
        flux *= float(chain.get("stoichiometry") or 1.0)

    return {
        "chain_id": chain["chain_id"],
        "voc_id": chain["voc_id"],
        "name": chain.get("name"),
        "flux": float(flux),
        "step_scores": step_scores,
        "enzyme_hits": sorted(set(enzyme_hits)),
        "pathway_drivers": sorted(set(pathway_drivers)),
        "detox_capacity": float(detox_capacity),
        "primary_tissues": list(chain.get("primary_tissues") or []),
    }


def score_all_chains(
    chains: list[dict[str, Any]],
    pathway_score_list: list[PathwayScore],
    mutated_genes: set[str],
) -> list[dict[str, Any]]:
    pw = {p.pathway_id: p.score for p in pathway_score_list}
    return [score_pathway_chain(c, pw, mutated_genes) for c in chains]
