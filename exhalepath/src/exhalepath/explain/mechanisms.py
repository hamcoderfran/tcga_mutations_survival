"""
Mechanism explainability: WHY a VOC is up/down for a disease.

Links each predicted VOC shift to:
  • dysregulated metabolic pathways (+ seed / driver genes)
  • biosynthetic VOC pathway chains
  • affected cell populations (atlas states + CELLxGENE Census fractions / top cell types)
  • category-default genetic alterations when patient mutations are not supplied

Builds per-disease mechanism packs for the full ~100-disease atlas.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from ..biomarker import BiomarkerReport, ExhaleBiomarkerEngine
from ..config import DATA_DIR, KNOWLEDGE_DIR
from ..knowledge.loader import KnowledgeBase, clear_knowledge_cache
from ..physio.census_fractions import (
    census_cell_state_fractions_for_disease,
    resolve_us_disease_id,
)
from ..schemas import VOCPrediction

CENSUS_DIR = DATA_DIR / "census"

# Category → canonical driver genes used when query has no mutations
CATEGORY_DRIVER_GENES: dict[str, list[str]] = {
    "cancer": ["KRAS", "TP53", "MYC", "HIF1A", "LDHA", "HK2", "EGFR", "PIK3CA"],
    "neurological": ["APOE", "SNCA", "MAPT", "PSEN1", "SOD1", "GRIN1", "IDO1"],
    "neurodegenerative": ["APOE", "SNCA", "MAPT", "PSEN1", "SOD1"],
    "metabolic": ["HMGCS2", "CPT1A", "INSR", "PPARG", "GCK"],
    "microbiome": ["IDO1", "TDO2", "AHR", "NOS2", "FMO3"],
    "inflammatory": ["TNF", "IL6", "PTGS2", "NOS2", "CYBB", "NFE2L2"],
    "infectious": ["TLR4", "MYD88", "NOS2", "CYBB"],
    "cardiovascular": ["APOE", "PCSK9", "NOS3", "HMGCR"],
    "pulmonary": ["CYP1A1", "CYP2E1", "NOS2", "CYBB", "HIF1A", "GPX4"],
}


@dataclass
class VocMechanism:
    voc_id: str
    name: str
    direction: str
    delta_ppb: float
    fold_change: float
    predicted_ppb: float
    healthy_ppb: float
    why: str
    pathways: list[dict[str, Any]] = field(default_factory=list)
    chains: list[str] = field(default_factory=list)
    cell_states: list[dict[str, Any]] = field(default_factory=list)
    driver_genes: list[str] = field(default_factory=list)
    genetic_alterations: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "voc_id": self.voc_id,
            "name": self.name,
            "direction": self.direction,
            "delta_ppb": self.delta_ppb,
            "fold_change": self.fold_change,
            "predicted_ppb": self.predicted_ppb,
            "healthy_ppb": self.healthy_ppb,
            "why": self.why,
            "pathways": self.pathways,
            "chains": self.chains,
            "cell_states": self.cell_states,
            "driver_genes": self.driver_genes,
            "genetic_alterations": self.genetic_alterations,
        }


class MechanismExplainer:
    def __init__(self, *, knowledge: KnowledgeBase | None = None, use_opentargets: bool = False):
        self.kb = knowledge or KnowledgeBase()
        self.engine = ExhaleBiomarkerEngine(
            knowledge=self.kb, use_opentargets=use_opentargets
        )
        self._celltype_table = self._load_celltypes()

    def _load_celltypes(self) -> pd.DataFrame | None:
        path = CENSUS_DIR / "disease_tissue_celltype_composition.csv"
        if not path.exists():
            return None
        return pd.read_csv(path)

    @staticmethod
    def _census_labels_for_us(us_id: str, disease_name: str | None = None) -> list[str]:
        match_path = CENSUS_DIR / "top100_disease_census_match.csv"
        if not match_path.exists() or not us_id:
            return []
        m = pd.read_csv(match_path)
        sub = m[m["us_disease_id"] == us_id]
        if sub.empty:
            return []
        labels: list[str] = []
        for raw in sub["census_diseases"].astype(str):
            for part in raw.split("|"):
                p = part.strip()
                if not p or p.lower() == "nan":
                    continue
                # Skip comorbidity combo labels ("dementia || diabetes …")
                if "||" in p:
                    continue
                labels.append(p)
        # Prefer labels whose text matches the disease / us_id
        text = f"{us_id.replace('_', ' ')} {disease_name or ''}".lower()
        tokens = {
            t
            for t in text.split()
            if len(t) > 3 and t not in {"type", "disease", "cancer", "syndrome"}
        }
        if tokens:
            preferred = [
                l for l in labels if any(t in l.lower() for t in tokens)
            ]
            if preferred:
                return list(dict.fromkeys(preferred))
        return list(dict.fromkeys(labels))

    def default_genes_for_disease(self, disease: dict[str, Any]) -> list[str]:
        cat = (disease.get("category") or "default").lower()
        genes = list(CATEGORY_DRIVER_GENES.get(cat, []))
        # Add seed genes from strongly biased pathways
        for pid, bias in (disease.get("pathway_bias") or {}).items():
            if float(bias) < 1.15:
                continue
            pw = self.kb.pathways.get(pid) or {}
            genes.extend(list(pw.get("seed_genes") or [])[:6])
        # unique preserve order
        seen = set()
        out = []
        for g in genes:
            gu = str(g).upper()
            if gu not in seen:
                seen.add(gu)
                out.append(gu)
        return out[:16]

    def census_context(
        self, disease_id: str, *, disease_name: str | None = None, preferred_tissue: str | None = None
    ) -> dict[str, Any]:
        us_id = resolve_us_disease_id(disease_id)
        fracs = census_cell_state_fractions_for_disease(disease_id)
        top_cts: list[dict[str, Any]] = []
        n_cells = 0
        labels_used: list[str] = []
        if us_id and self._celltype_table is not None:
            labels = self._census_labels_for_us(us_id, disease_name=disease_name)
            labels_used = labels
            sub = self._celltype_table
            if labels:
                sub = sub[sub["disease"].isin(labels)]
            else:
                sub = sub.iloc[0:0]
            if not sub.empty and "n_cells" in sub.columns:
                tissue_col = "tissue_general" if "tissue_general" in sub.columns else "tissue"
                # Prefer anatomic site of the disease when available
                if preferred_tissue:
                    pt = preferred_tissue.lower()
                    tissue_hit = sub[
                        sub[tissue_col].astype(str).str.lower().str.contains(pt, na=False)
                    ]
                    if not tissue_hit.empty and tissue_hit["n_cells"].sum() >= 500:
                        sub = tissue_hit
                n_cells = int(sub["n_cells"].sum())
                g = (
                    sub.groupby([tissue_col, "cell_type"], dropna=False)["n_cells"]
                    .sum()
                    .reset_index()
                    .sort_values("n_cells", ascending=False)
                    .head(12)
                )
                total = float(g["n_cells"].sum()) or 1.0
                for r in g.itertuples(index=False):
                    top_cts.append(
                        {
                            "tissue": getattr(r, tissue_col),
                            "cell_type": r.cell_type,
                            "n_cells": int(r.n_cells),
                            "fraction": float(r.n_cells) / total,
                        }
                    )
        return {
            "us_disease_id": us_id,
            "census_labels_used": labels_used[:12],
            "n_census_cells": n_cells,
            "state_fractions": fracs,
            "top_cell_types": top_cts,
        }

    def explain_voc(
        self,
        *,
        voc: VOCPrediction,
        disease: dict[str, Any],
        pathway_scores: list[Any],
        cell_states: list[Any],
        genes: list[str],
        census: dict[str, Any],
    ) -> VocMechanism:
        score_map = {p.pathway_id: p for p in pathway_scores}
        # Pathways that emit this VOC
        pw_hits = []
        for pid, pdef in self.kb.pathways.items():
            coef = float((pdef.get("voc_effects") or {}).get(voc.voc_id, 0.0))
            if abs(coef) < 1e-9:
                continue
            ps = score_map.get(pid)
            score = float(ps.score) if ps else 0.0
            seeds = [g.upper() for g in (pdef.get("seed_genes") or [])]
            hit_genes = sorted(set(seeds) & set(g.upper() for g in genes))
            if ps and ps.hit_genes:
                hit_genes = sorted(set(hit_genes) | {h.upper() for h in ps.hit_genes})
            pw_hits.append(
                {
                    "pathway_id": pid,
                    "name": pdef.get("name"),
                    "score": score,
                    "emission_coef": coef,
                    "contribution": score * coef,
                    "seed_genes": seeds[:12],
                    "hit_genes": hit_genes,
                }
            )
        pw_hits.sort(key=lambda x: abs(x["contribution"]), reverse=True)
        pw_hits = pw_hits[:6]

        chains = []
        if voc.physiology and voc.physiology.contributing_chains:
            chains = list(voc.physiology.contributing_chains)
        else:
            chains = [
                c.get("chain_id")
                for c in self.kb.pathway_chains
                if c.get("voc_id") == voc.voc_id
            ][:4]

        # Cell states that produce those chains or are disease-modulated
        state_rows = []
        chain_set = set(chains)
        for st in cell_states:
            produces = set(st.produces_chains or [])
            relevant = bool(produces & chain_set) or st.disease_modulated or st.marker_hit_score > 0
            if not relevant and st.effective_source < 0.02:
                continue
            if not relevant and not (produces & chain_set):
                # keep top effective sources lightly
                if st.effective_source < 0.05:
                    continue
            census_frac = (census.get("state_fractions") or {}).get(st.state_id)
            state_rows.append(
                {
                    "state_id": st.state_id,
                    "name": st.name,
                    "tissue": st.tissue,
                    "density": st.density,
                    "activity": st.activity,
                    "effective_source": st.effective_source,
                    "marker_genes_hit": list(st.marker_genes_hit or []),
                    "produces_chains": list(st.produces_chains or []),
                    "census_fraction": census_frac,
                    "disease_modulated": st.disease_modulated,
                }
            )
        state_rows.sort(key=lambda x: x["effective_source"], reverse=True)
        state_rows = state_rows[:8]

        driver_genes = []
        for p in pw_hits:
            driver_genes.extend(p.get("hit_genes") or [])
            if not p.get("hit_genes"):
                driver_genes.extend((p.get("seed_genes") or [])[:3])
        for st in state_rows:
            driver_genes.extend(st.get("marker_genes_hit") or [])
        # unique
        seen = set()
        drivers = []
        for g in driver_genes:
            if g not in seen:
                seen.add(g)
                drivers.append(g)

        genetic = []
        gene_set = set(g.upper() for g in genes)
        for g in drivers:
            role = "query_or_category_driver" if g in gene_set else "pathway_seed"
            genetic.append(
                {
                    "gene": g,
                    "role": role,
                    "source": "mutated_or_category_default"
                    if g in gene_set
                    else "pathway_seed_gene",
                }
            )

        direction = (
            "up" if voc.delta_ppb > 0.05 else "down" if voc.delta_ppb < -0.05 else "near_baseline"
        )
        why = _compose_why(
            voc_name=voc.name,
            direction=direction,
            disease_name=disease.get("name") or disease.get("disease_id"),
            pathways=pw_hits,
            cell_states=state_rows,
            chains=chains,
            drivers=drivers,
            census=census,
        )
        return VocMechanism(
            voc_id=voc.voc_id,
            name=voc.name,
            direction=direction,
            delta_ppb=float(voc.delta_ppb),
            fold_change=float(voc.fold_change),
            predicted_ppb=float(voc.predicted_ppb),
            healthy_ppb=float(voc.healthy_ppb),
            why=why,
            pathways=pw_hits,
            chains=list(chains),
            cell_states=state_rows,
            driver_genes=drivers[:12],
            genetic_alterations=genetic[:16],
        )

    def build_disease_pack(
        self,
        disease_query: str,
        *,
        location: str | None = None,
        genes: list[str] | None = None,
        top_n: int = 20,
        mode: str = "hybrid",
    ) -> dict[str, Any]:
        disease = self.kb.resolve_disease(disease_query)
        use_genes = genes or self.default_genes_for_disease(disease)
        loc = location or disease.get("default_site") or "systemic"
        report = self.engine.predict(
            disease_query,
            location=loc,
            genes=use_genes,
            top_n=max(top_n, 50),
            mode=mode,
            explain=False,  # avoid recursive explanation
        )
        site = (loc or "").split()[0] if loc else None
        census = self.census_context(
            report.disease_id,
            disease_name=report.disease_name,
            preferred_tissue=site,
        )
        if not census.get("state_fractions"):
            us = resolve_us_disease_id(report.disease_id)
            if us:
                census = self.census_context(
                    us, disease_name=report.disease_name, preferred_tissue=site
                )

        pathway_scores = report.result.bundle.pathway_scores
        cell_states = report.result.bundle.cell_states
        mechanisms = []
        for voc in report.top_vocs[:top_n]:
            if abs(voc.delta_ppb) < 0.01 and abs(voc.log2_fold_change) < 0.05:
                continue
            mechanisms.append(
                self.explain_voc(
                    voc=voc,
                    disease=disease,
                    pathway_scores=pathway_scores,
                    cell_states=cell_states,
                    genes=use_genes,
                    census=census,
                ).to_dict()
            )

        pack = {
            "disease_id": report.disease_id,
            "disease_name": report.disease_name,
            "category": disease.get("category"),
            "location": report.location,
            "genes_used": use_genes,
            "genes_source": "user" if genes else "category_and_pathway_defaults",
            "census": census,
            "top_pathways": [
                {
                    "pathway_id": p.pathway_id,
                    "name": p.name,
                    "score": p.score,
                    "hit_genes": p.hit_genes,
                    "disease_bias": p.disease_bias,
                }
                for p in pathway_scores[:12]
            ],
            "cell_states": [
                {
                    "state_id": s.state_id,
                    "name": s.name,
                    "tissue": s.tissue,
                    "density": s.density,
                    "activity": s.activity,
                    "effective_source": s.effective_source,
                    "marker_genes_hit": s.marker_genes_hit,
                    "produces_chains": s.produces_chains,
                    "census_fraction": (census.get("state_fractions") or {}).get(s.state_id),
                }
                for s in cell_states[:12]
            ],
            "voc_mechanisms": mechanisms,
            "model_version": report.model_version,
            "n_vocs_explained": len(mechanisms),
        }
        return pack


def _compose_why(
    *,
    voc_name: str,
    direction: str,
    disease_name: str,
    pathways: list[dict[str, Any]],
    cell_states: list[dict[str, Any]],
    chains: list[str],
    drivers: list[str],
    census: dict[str, Any],
) -> str:
    dir_word = {"up": "elevated", "down": "reduced", "near_baseline": "near baseline"}[direction]
    bits = [f"{voc_name} is predicted {dir_word} in {disease_name}"]
    if pathways:
        top = pathways[0]
        genes = ", ".join((top.get("hit_genes") or top.get("seed_genes") or [])[:4]) or "pathway enzymes"
        bits.append(
            f"primarily via {top.get('name')} (score={top.get('score', 0):.2f}; genes {genes})"
        )
    if cell_states:
        st = cell_states[0]
        frac = st.get("census_fraction")
        frac_s = f", Census fraction≈{frac:.3f}" if isinstance(frac, float) else ""
        bits.append(
            f"sourced from {st.get('name')} cells in {st.get('tissue')} "
            f"(density×activity={st.get('effective_source', 0):.3f}{frac_s})"
        )
    if chains:
        bits.append(f"biosynthetic chain(s): {', '.join(chains[:3])}")
    if drivers:
        bits.append(f"key genetic/pathway nodes: {', '.join(drivers[:6])}")
    if census.get("top_cell_types"):
        ct = census["top_cell_types"][0]
        bits.append(
            f"Census enriched population example: {ct.get('cell_type')} in {ct.get('tissue')}"
        )
    return "; ".join(bits) + "."


def build_all_mechanism_packs(
    *,
    out_path: Path | None = None,
    top_n_vocs: int = 15,
    mode: str = "hybrid",
) -> Path:
    """Build mechanism packs for every disease in the atlas (~100)."""
    clear_knowledge_cache()
    explainer = MechanismExplainer(use_opentargets=False)
    packs = []
    for did in sorted(explainer.kb.diseases.keys()):
        packs.append(
            explainer.build_disease_pack(did, top_n=top_n_vocs, mode=mode)
        )
    doc = {
        "version": "1.0.0",
        "description": (
            "Per-disease exhaled VOC mechanism packs: pathways, Census cell populations, "
            "driver genes, and plain-language WHY statements for top VOC changes."
        ),
        "n_diseases": len(packs),
        "diseases": packs,
    }
    out_path = Path(out_path or KNOWLEDGE_DIR / "disease_mechanism_packs.json")
    out_path.write_text(json.dumps(doc, indent=2))
    # also write slim index
    index = [
        {
            "disease_id": p["disease_id"],
            "disease_name": p["disease_name"],
            "category": p.get("category"),
            "n_vocs_explained": p.get("n_vocs_explained"),
            "has_census": bool((p.get("census") or {}).get("state_fractions")),
            "n_census_cells": (p.get("census") or {}).get("n_census_cells"),
        }
        for p in packs
    ]
    idx_path = out_path.with_name("disease_mechanism_index.json")
    idx_path.write_text(json.dumps({"n_diseases": len(index), "diseases": index}, indent=2))
    return out_path
