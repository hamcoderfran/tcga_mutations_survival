"""Zero-shot mechanism resolver for custom:: / unresolved diseases."""

from __future__ import annotations

from typing import Any

import re

from .category_templates import (
    CATEGORY_DEFAULT_SITE,
    CATEGORY_DRIVER_GENES,
    GENERIC_DISEASE_TOKENS,
    infer_site_from_text,
    match_category_keywords,
    match_token_cues,
    merge_pathway_bias,
    template_for_category,
)
from .disease_gene_index import lookup_disease_genes, lookup_zero_shot_record
from .ontology_nn import ontology_transfer


def _is_gene_like(raw: str) -> bool:
    s = str(raw).strip()
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9-]{1,14}", s):
        return False
    return bool(re.search(r"\d", s)) or s.upper() in {
        "BRCA",
        "KRAS",
        "EGFR",
        "TP53",
        "PTEN",
        "MYC",
        "ALK",
        "MET",
    }


def _is_generic_only(query: str, description: str | None) -> bool:
    """True when the query is only a generic token and no phenotype text was given."""
    if description and str(description).strip():
        return False
    toks = re.sub(r"[^a-z0-9]+", " ", str(query).lower()).split()
    return bool(toks) and all(t in GENERIC_DISEASE_TOKENS for t in toks)


def enrich_unresolved_disease(
    base: dict[str, Any],
    *,
    diseases: dict[str, dict[str, Any]],
    datasources: dict[str, Any] | None = None,
    location_hint: str | None = None,
    description: str | None = None,
    copy_voc_priors: bool = False,
) -> dict[str, Any]:
    """
    Fill category / site / pathway_bias / driver_genes for an unresolved disease.

    Always keeps voc_log2fc_prior empty unless copy_voc_priors=True.
    Sets _zero_shot metadata and mechanism confidence.
    """
    out = dict(base)
    query = str(out.get("name") or "")
    # Gene symbols must not ontology-NN into a cancer (BRCA1 ↛ breast VOC panel)
    if _is_gene_like(query):
        out["voc_log2fc_prior"] = {}
        out["_unresolved"] = True
        out["_zero_shot"] = True
        out["_zero_shot_evidence"] = [{"source": "gene_symbol_guard", "query": query}]
        out["_mechanism_confidence"] = 0.0
        out["_zero_shot_mode"] = "near_healthy_fallback"
        out["_associated_gene_scores"] = {query.upper(): 0.9}
        out["driver_genes"] = [query.upper()]
        out["_cell_state_donor_ids"] = []
        return out

    # Bare "cancer"/"disease" alone stays near-healthy; need a more specific phrase
    if _is_generic_only(query, description):
        out["voc_log2fc_prior"] = {}
        out["_unresolved"] = True
        out["_zero_shot"] = True
        out["_zero_shot_evidence"] = [{"source": "generic_token_guard", "query": query}]
        out["_mechanism_confidence"] = 0.0
        out["_zero_shot_mode"] = "near_healthy_fallback"
        out["_associated_gene_scores"] = {}
        out["driver_genes"] = []
        out["_cell_state_donor_ids"] = []
        return out

    text = " ".join(x for x in [query, description or "", location_hint or ""] if x)
    evidence: list[dict[str, Any]] = []
    pathway_bias: dict[str, float] = {}
    genes: list[str] = []
    gene_scores: dict[str, float] = {}
    category: str | None = None
    site: str | None = None
    mondo_id = out.get("mondo_id")
    donors: list[str] = []
    confidence = 0.0

    # 1) Curated rare-disease catalog (strongest offline signal)
    zs = lookup_zero_shot_record(query)
    if zs:
        category = zs.get("category") or category
        site = zs.get("default_site") or site
        pathway_bias = merge_pathway_bias(pathway_bias, zs.get("pathway_bias"))
        if not pathway_bias and category:
            pathway_bias = merge_pathway_bias(pathway_bias, template_for_category(category))
        for g in zs.get("genes") or []:
            if isinstance(g, dict):
                gu = str(g.get("gene") or "").upper()
                gene_scores[gu] = max(gene_scores.get(gu, 0.0), float(g.get("score") or 0.75))
            else:
                gu = str(g).upper()
                gene_scores[gu] = max(gene_scores.get(gu, 0.0), 0.75)
            if gu:
                genes.append(gu)
        mondo_id = zs.get("mondo_id") or mondo_id
        if zs.get("cell_state_donor_ids"):
            donors.extend(list(zs["cell_state_donor_ids"]))
        elif zs.get("atlas_neighbor"):
            donors.append(str(zs["atlas_neighbor"]))
        confidence = max(confidence, float(zs.get("confidence") or 0.72))
        evidence.append(
            {
                "source": "zero_shot_disease_genes",
                "disease_id": zs.get("disease_id"),
                "strength": confidence,
            }
        )

    # 2) Mechanism token cues
    for cue in match_token_cues(text):
        category = category or cue["category"]
        pathway_bias = merge_pathway_bias(pathway_bias, cue.get("pathway_bias"))
        site = site or cue.get("default_site")
        genes.extend(list(cue.get("driver_genes") or []))
        for g in cue.get("driver_genes") or []:
            gene_scores[str(g).upper()] = max(gene_scores.get(str(g).upper(), 0.0), 0.65)
        strength = float(cue.get("strength") or 0.55)
        confidence = max(confidence, strength)
        evidence.append(
            {
                "source": "token_cue",
                "matched": cue.get("matched"),
                "category": cue.get("category"),
                "strength": strength,
            }
        )

    # 3) Category keywords → pathway template (no VOC priors)
    cat_hits = match_category_keywords(text)
    if cat_hits:
        # Prefer first non-default clinical category
        pick = cat_hits[0]
        if category is None:
            category = pick
        pathway_bias = merge_pathway_bias(pathway_bias, template_for_category(pick))
        for g in CATEGORY_DRIVER_GENES.get(pick, [])[:6]:
            genes.append(g)
            gene_scores[g.upper()] = max(gene_scores.get(g.upper(), 0.0), 0.5)
        strength = 0.38 if len(cat_hits) == 1 else 0.45
        confidence = max(confidence, strength)
        evidence.append(
            {
                "source": "category_keyword",
                "categories": cat_hits,
                "strength": strength,
            }
        )

    # 4) Ontology / atlas nearest neighbor (mechanistic transfer only)
    ont = ontology_transfer(query, diseases, copy_voc_priors=copy_voc_priors)
    if ont:
        if category is None or category in {"unspecified", "default"}:
            category = ont.get("category") or category
        site = site or ont.get("default_site")
        pathway_bias = merge_pathway_bias(pathway_bias, ont.get("pathway_bias"))
        genes.extend(list(ont.get("driver_genes") or []))
        for g in ont.get("driver_genes") or []:
            gene_scores[str(g).upper()] = max(
                gene_scores.get(str(g).upper(), 0.0), 0.55 * float(ont.get("neighbor_score") or 0.5)
            )
        donors.extend(list(ont.get("cell_state_donor_ids") or []))
        mondo_id = mondo_id or ont.get("mondo_id")
        strength = 0.4 + 0.45 * float(ont.get("neighbor_score") or 0.0)
        confidence = max(confidence, min(0.85, strength))
        evidence.append(
            {
                "source": "ontology_nn",
                "neighbor_disease_id": ont.get("neighbor_disease_id"),
                "neighbor_name": ont.get("neighbor_name"),
                "neighbor_score": ont.get("neighbor_score"),
                "neighbors": ont.get("neighbors"),
                "strength": strength,
                "copied_voc_priors": bool(copy_voc_priors),
            }
        )

    # 5) Free-text gene index (OT + zero-shot catalog)
    for g, s in lookup_disease_genes(
        query, kb_datasources=datasources, atlas_disease_id=None
    ).items():
        gene_scores[g] = max(gene_scores.get(g, 0.0), float(s))
        genes.append(g)
        confidence = max(confidence, min(0.8, 0.5 + 0.3 * float(s)))
    if gene_scores:
        evidence.append(
            {
                "source": "disease_gene_index",
                "n_genes": len(gene_scores),
                "top_genes": sorted(gene_scores, key=gene_scores.get, reverse=True)[:8],
            }
        )

    # 6) Site from disease text / location hint / category default
    site = (
        site
        or infer_site_from_text(text)
        or (infer_site_from_text(location_hint) if location_hint else None)
        or CATEGORY_DEFAULT_SITE.get(category or "default")
    )

    # If we have any mechanism signal, ensure category template floor
    if confidence >= 0.34 and category:
        pathway_bias = merge_pathway_bias(pathway_bias, template_for_category(category))
        if not genes:
            genes.extend(CATEGORY_DRIVER_GENES.get(category, [])[:6])
            for g in CATEGORY_DRIVER_GENES.get(category, [])[:6]:
                gene_scores[g.upper()] = max(gene_scores.get(g.upper(), 0.0), 0.45)

    # Dedupe genes
    seen: set[str] = set()
    driver_genes: list[str] = []
    for g in genes:
        gu = str(g).upper().replace(" ", "")
        if not gu or gu in seen:
            continue
        seen.add(gu)
        driver_genes.append(gu)

    out["category"] = category or out.get("category") or "unspecified"
    out["default_site"] = site
    out["pathway_bias"] = pathway_bias if confidence >= 0.34 else {}
    out["driver_genes"] = driver_genes[:16] if confidence >= 0.34 else []
    out["mondo_id"] = mondo_id
    if not copy_voc_priors:
        out["voc_log2fc_prior"] = {}
    out["_unresolved"] = True
    out["_zero_shot"] = True
    out["_zero_shot_evidence"] = evidence
    out["_mechanism_confidence"] = float(confidence)
    out["_associated_gene_scores"] = {
        g: float(gene_scores[g]) for g in driver_genes[:16] if g in gene_scores
    }
    out["_cell_state_donor_ids"] = list(dict.fromkeys(d for d in donors if d))
    out["_zero_shot_mode"] = (
        "mechanism"
        if confidence >= 0.34 and (pathway_bias or driver_genes)
        else "near_healthy_fallback"
    )
    return out
