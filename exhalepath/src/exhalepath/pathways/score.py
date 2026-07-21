from __future__ import annotations

from typing import Iterable

from ..knowledge.loader import KnowledgeBase
from ..schemas import PathwayScore, TumorContext
from ..utils import stage_multiplier


def _site_multiplier(pathway_id: str, tumor: TumorContext | None, default_site: str | None) -> float:
    site = (tumor.primary_site if tumor and tumor.primary_site else default_site or "").lower()
    if not site:
        return 1.0
    boosts = {
        "liver": {
            "methionine_transsulfuration": 1.35,
            "urea_cycle": 1.3,
            "cytochrome_p450_detox": 1.25,
        },
        "lung": {
            "lipid_peroxidation": 1.25,
            "cytochrome_p450_detox": 1.2,
            "Kras_mapk_proliferation": 1.15,
        },
        "pancreas": {"Kras_mapk_proliferation": 1.3, "glycolysis_warburg": 1.2},
        "colon": {
            "glycolysis_warburg": 1.15,
            "one_carbon_folate": 1.1,
            "gut_microbiome_fermentation": 1.25,
            "microbial_proteolysis_putrefaction": 1.2,
        },
        "gut": {
            "gut_microbiome_fermentation": 1.4,
            "microbial_proteolysis_putrefaction": 1.35,
        },
        "intestin": {  # intestine / small_intestine
            "gut_microbiome_fermentation": 1.35,
            "microbial_proteolysis_putrefaction": 1.3,
        },
        "stomach": {"urea_cycle": 1.25, "gut_microbiome_fermentation": 1.15},
        "brain": {
            "neuroinflammation": 1.35,
            "neurotransmitter_metabolism": 1.3,
            "brain_energy_metabolism": 1.3,
            "lipid_peroxidation": 1.15,
        },
        "breast": {"pi3k_akt_mtor": 1.2, "lipid_peroxidation": 1.1},
        "kidney": {"urea_cycle": 1.2},
        "ovary": {"lipid_peroxidation": 1.15, "glycolysis_warburg": 1.1},
        "heart": {
            "lipid_peroxidation": 1.2,
            "fatty_acid_oxidation": 1.25,
            "mevalonate_cholesterol": 1.15,
        },
        "prostate": {"lipid_peroxidation": 1.1, "glycolysis_warburg": 1.1},
        "skin": {"lipid_peroxidation": 1.15, "apoptosis_necrosis": 1.1},
        "adipose": {"fatty_acid_oxidation": 1.35, "ketone_body_metabolism": 1.2},
        "blood": {"glycolysis_warburg": 1.15, "apoptosis_necrosis": 1.15},
        "bladder": {"lipid_peroxidation": 1.1, "urea_cycle": 1.1},
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
    mutated = {g.upper() for g in mutated_genes if g}
    associated_genes = associated_genes or {}
    overrides = pathway_overrides or {}
    stage_m = stage_multiplier(tumor.stage if tumor else None)
    burden = 1.0
    if tumor and tumor.tumor_burden_proxy is not None:
        burden = 0.85 + 0.5 * tumor.tumor_burden_proxy
    if tumor and tumor.metastatic:
        burden *= 1.15

    # When no molecular evidence is available, apply a mild disease-level baseline
    # so curated disease priors (e.g. T2D acetone) are not wiped out by zero pathway scores.
    has_molecular = bool(mutated) or bool(associated_genes) or bool(overrides)

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
            gu = str(g).upper()
            if gu in genes:
                soft += float(w)
                soft_hits.append(gu)
        soft = min(soft / max(len(genes), 1), 1.0)

        base = hit_frac + 0.5 * soft
        if pid in overrides:
            base = float(overrides[pid])

        bias = float(disease.get("pathway_bias", {}).get(pid, 1.0))
        # Disease-atlas pathway_bias always contributes a floor so VOC panel members
        # linked to biased pathways activate even when supplied genes hit other sets.
        if bias > 1.0 and pid not in overrides:
            prior_floor = 0.32 * (bias - 1.0)
            if not has_molecular:
                base = max(base, prior_floor)
            else:
                base = base + prior_floor

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
