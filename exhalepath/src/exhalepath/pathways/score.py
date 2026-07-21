from __future__ import annotations

from typing import Iterable

from ..knowledge.loader import KnowledgeBase
from ..schemas import PathwayScore, TumorContext


def _stage_multiplier(tumor: TumorContext | None) -> float:
    if tumor is None or not tumor.stage:
        return 1.0
    s = tumor.stage.upper().replace("STAGE", "").strip()
    mapping = {"0": 0.85, "I": 1.0, "IA": 1.0, "IB": 1.05, "II": 1.15, "IIA": 1.15,
               "IIB": 1.2, "III": 1.35, "IIIA": 1.35, "IIIB": 1.4, "IV": 1.55}
    for key, val in mapping.items():
        if s.startswith(key):
            return val
    # Roman / numeric fallback
    for roman, val in [("IV", 1.55), ("III", 1.35), ("II", 1.15), ("I", 1.0)]:
        if roman in s:
            return val
    return 1.0


def _site_multiplier(pathway_id: str, tumor: TumorContext | None, default_site: str | None) -> float:
    site = (tumor.primary_site if tumor and tumor.primary_site else default_site or "").lower()
    if not site:
        return 1.0
    boosts = {
        "liver": {"methionine_transsulfuration": 1.35, "urea_cycle": 1.3, "cytochrome_p450_detox": 1.25},
        "lung": {"lipid_peroxidation": 1.25, "cytochrome_p450_detox": 1.2, "Kras_mapk_proliferation": 1.15},
        "pancreas": {"Kras_mapk_proliferation": 1.3, "glycolysis_warburg": 1.2},
        "colon": {"glycolysis_warburg": 1.15, "one_carbon_folate": 1.1},
        "breast": {"pi3k_akt_mtor": 1.2, "lipid_peroxidation": 1.1},
        "kidney": {"urea_cycle": 1.2},
        "ovary": {"lipid_peroxidation": 1.15, "glycolysis_warburg": 1.1},
    }
    for key, pathway_boost in boosts.items():
        if key in site:
            return float(pathway_boost.get(pathway_id, 1.0))
    return 1.0


def score_pathways(
    *,
    kb: KnowledgeBase,
    disease: dict,
    mutated_genes: Iterable[str],
    tumor: TumorContext | None = None,
    pathway_overrides: dict[str, float] | None = None,
    associated_genes: dict[str, float] | None = None,
) -> list[PathwayScore]:
    """
    Score metabolic pathway dysregulation from mutated / associated genes.

    Score ≈ disease_bias * site_mult * stage_mult * (mutation hits + soft OT associations)
    """
    mutated = {g.upper() for g in mutated_genes}
    associated_genes = associated_genes or {}
    overrides = pathway_overrides or {}
    stage_m = _stage_multiplier(tumor)
    burden = 1.0
    if tumor and tumor.tumor_burden_proxy is not None:
        burden = 0.85 + 0.5 * tumor.tumor_burden_proxy
    if tumor and tumor.metastatic:
        burden *= 1.15

    scores: list[PathwayScore] = []
    for pid, p in kb.pathways.items():
        genes = {g.upper() for g in p.get("seed_genes", [])}
        if not genes:
            continue
        hits = sorted(genes & mutated)
        hit_frac = len(hits) / len(genes)
        soft = 0.0
        soft_hits = []
        for g, w in associated_genes.items():
            if g.upper() in genes:
                soft += float(w)
                soft_hits.append(g.upper())
        soft = min(soft / max(len(genes), 1), 1.0)

        base = hit_frac + 0.5 * soft
        if pid in overrides:
            base = float(overrides[pid])

        bias = float(disease.get("pathway_bias", {}).get(pid, 1.0))
        site_m = _site_multiplier(pid, tumor, disease.get("default_site"))
        score = base * bias * stage_m * site_m * burden

        # Mild histology modulation
        hist = (tumor.histology if tumor and tumor.histology else "").lower()
        if "squamous" in hist and pid == "lipid_peroxidation":
            score *= 1.1
        if "adenocarcinoma" in hist and pid == "glycolysis_warburg":
            score *= 1.05

        scores.append(
            PathwayScore(
                pathway_id=pid,
                name=p["name"],
                score=float(score),
                hit_genes=sorted(set(hits) | set(soft_hits)),
                disease_bias=bias,
            )
        )

    scores.sort(key=lambda x: x.score, reverse=True)
    return scores
