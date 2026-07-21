"""
ExhalePath Biomarker Engine — whole-body exhaled VOC biomarker prediction.

Analogous in ambition to structure predictors for proteins: given a disease and
any anatomic location of affected cells, return the ranked exhaled VOC changes
(top 50) with predicted quantities in ppb.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from .body.tissues import WholeBodyMap, resolve_location
from .knowledge.loader import KnowledgeBase, clear_knowledge_cache, default_knowledge
from .model.predict import ExhalePathPredictor, PredictionResult
from .schemas import DiseaseQuery, TumorContext, VOCPrediction


@dataclass
class BiomarkerReport:
    """Top-N exhaled VOC biomarker prediction for a disease × location."""

    disease_query: str
    disease_id: str
    disease_name: str
    location_query: str
    location: dict[str, Any]
    top_vocs: list[VOCPrediction]
    n_vocs_modeled: int
    model_version: str
    notes: list[str]
    result: PredictionResult
    mechanisms: list[dict[str, Any]] | None = None
    census: dict[str, Any] | None = None

    @property
    def ranked(self) -> list[VOCPrediction]:
        return self.top_vocs

    @property
    def tissue(self) -> str:
        return str(self.location.get("tissue_id") or "")

    @property
    def resolved_location(self) -> str:
        return str(self.location.get("name") or self.location_query)

    def to_dataframe(self) -> pd.DataFrame:
        rows = []
        for i, p in enumerate(self.top_vocs, 1):
            rows.append(
                {
                    "rank": i,
                    "voc_id": p.voc_id,
                    "name": p.name,
                    "cas": p.cas,
                    "healthy_ppb": p.healthy_ppb,
                    "predicted_ppb": p.predicted_ppb,
                    "delta_ppb": p.delta_ppb,
                    "fold_change": p.fold_change,
                    "log2_fold_change": p.log2_fold_change,
                    "ci_low_ppb": p.ci_low_ppb,
                    "ci_high_ppb": p.ci_high_ppb,
                    "confidence": p.confidence,
                    "drivers": ",".join(p.top_pathway_drivers),
                    "tissue_production": None
                    if not p.physiology
                    else p.physiology.tissue_production,
                    "blood_delivery": None if not p.physiology else p.physiology.blood_delivery,
                    "alveolar_fraction": None
                    if not p.physiology
                    else p.physiology.alveolar_fraction,
                }
            )
        return pd.DataFrame(rows)

    def summary(self) -> str:
        lines = [
            f"ExhalePath Biomarker Report",
            f"  Disease:  {self.disease_name} ({self.disease_id})",
            f"  Location: {self.location.get('name')} [{self.location.get('tissue_id')}]",
            f"  Modeled:  {self.n_vocs_modeled} VOCs · showing top {len(self.top_vocs)}",
            f"  Model:    {self.model_version}",
            "",
            f"{'Rank':<5}{'VOC':<22}{'Healthy':>10}{'Predicted':>12}{'Δ ppb':>12}{'Fold':>8}",
        ]
        for i, p in enumerate(self.top_vocs, 1):
            lines.append(
                f"{i:<5}{p.name[:21]:<22}{p.healthy_ppb:10.2f}{p.predicted_ppb:12.2f}"
                f"{p.delta_ppb:12.2f}{p.fold_change:8.2f}x"
            )
        if self.mechanisms:
            lines.append("")
            lines.append("Why (top mechanisms):")
            for m in self.mechanisms[:5]:
                lines.append(f"  • {m.get('why')}")
        return "\n".join(lines)


class ExhaleBiomarkerEngine:
    """
    Whole-body exhaled biomarker predictor.

    Parameters
    ----------
    disease : any of ~100 atlas diseases (or free-text / Open Targets)
    location : any body tissue / organ (Census whole-body map + free text)
    affected_fraction : optional density of affected cells at that location (0–1)
    """

    def __init__(
        self,
        *,
        knowledge: KnowledgeBase | None = None,
        use_opentargets: bool = False,
        reload_knowledge: bool = False,
    ):
        if reload_knowledge:
            clear_knowledge_cache()
        self.kb = knowledge or default_knowledge()
        self.body = WholeBodyMap()
        self.predictor = ExhalePathPredictor(
            knowledge=self.kb, use_opentargets=use_opentargets
        )

    def list_diseases(self) -> list[dict[str, Any]]:
        return [
            {
                "disease_id": d["disease_id"],
                "name": d["name"],
                "category": d.get("category"),
                "default_site": d.get("default_site"),
                "aliases": d.get("aliases") or [],
            }
            for d in self.kb.diseases.values()
        ]

    def list_locations(self) -> list[dict[str, Any]]:
        return self.body.list_locations()

    def predict(
        self,
        disease: str,
        location: str | None = None,
        *,
        top_n: int = 50,
        stage: str | None = None,
        genes: list[str] | None = None,
        affected_fraction: float | None = None,
        affected_activity: float | None = None,
        mode: str = "hybrid",
        sex: str | None = None,
        age_years: float | None = None,
        smoking_status: str | None = None,
        metastatic: bool = False,
        explain: bool = True,
        comorbidities: list[str] | None = None,
        comorbidity_weight: float = 0.65,
    ) -> BiomarkerReport:
        loc = resolve_location(location) if location else resolve_location(
            self.kb.resolve_disease(disease).get("default_site") or "systemic"
        )
        disease_obj = self.kb.resolve_disease(disease)

        # Place affected cells at the requested location
        cell_activity = {}
        cell_fractions = {}
        if affected_fraction is not None or affected_activity is not None:
            # Bias dominant tissue-linked states
            tissue = (loc.get("tissue_id") or "").lower()
            state_hints = _states_for_tissue(tissue, disease_obj.get("category"))
            frac = float(affected_fraction) if affected_fraction is not None else 0.35
            act = float(affected_activity) if affected_activity is not None else 1.0
            for sid in state_hints:
                cell_fractions[sid] = max(frac, 0.05)
                cell_activity[sid] = act

        tumor = TumorContext(
            primary_site=loc.get("name") or loc.get("tissue_id"),
            stage=stage,
            metastatic=metastatic,
            tumor_burden_proxy=affected_fraction,
        )
        query = DiseaseQuery(
            disease=disease,
            tumor=tumor,
            mutated_genes=genes or [],
            cell_state_fractions=cell_fractions,
            cell_state_activity=cell_activity,
            mode=mode,  # type: ignore[arg-type]
            sex=sex if sex in {"female", "male", "other"} else None,
            age_years=age_years,
            smoking_status=smoking_status
            if smoking_status in {"never", "former", "current"}
            else None,
            comorbidities=list(comorbidities or []),
            comorbidity_weight=float(comorbidity_weight),
        )
        result = self.predictor.predict(query)
        # Rank by absolute delta ppb (quantity change), then |log2fc|
        ranked = sorted(
            result.bundle.predictions,
            key=lambda p: (abs(p.delta_ppb), abs(p.log2_fold_change)),
            reverse=True,
        )
        top = ranked[: max(1, min(top_n, len(ranked)))]
        notes = list(result.bundle.notes)
        notes.insert(
            0,
            "ExhalePath Biomarker Engine: whole-body location-aware exhaled VOC panel "
            f"({len(result.bundle.predictions)} VOCs modeled).",
        )
        if not loc.get("matched"):
            notes.append(
                f"Location '{location}' not in Census tissue index; used as free-text primary site."
            )

        mechanisms = None
        census = None
        if explain:
            from .explain.mechanisms import MechanismExplainer
            from .knowledge.comorbidity import (
                merge_disease_with_comorbidities,
                resolve_comorbid_diseases,
            )

            explainer = MechanismExplainer(knowledge=self.kb, use_opentargets=False)
            fused = disease_obj
            if comorbidities:
                fused = merge_disease_with_comorbidities(
                    disease_obj,
                    resolve_comorbid_diseases(self.kb, list(comorbidities)),
                    weight=float(comorbidity_weight),
                )
            use_genes = genes or explainer.default_genes_for_disease(fused)
            site = (loc.get("tissue_id") or loc.get("name") or "").split()[0]
            census = explainer.census_context(
                result.bundle.disease_id,
                disease_name=result.bundle.disease_name,
                preferred_tissue=site,
            )
            mechanisms = [
                explainer.explain_voc(
                    voc=p,
                    disease=fused,
                    pathway_scores=result.bundle.pathway_scores,
                    cell_states=result.bundle.cell_states,
                    genes=use_genes,
                    census=census,
                ).to_dict()
                for p in top
                if abs(p.delta_ppb) >= 0.01 or abs(p.log2_fold_change) >= 0.05
            ][: min(15, len(top))]
            if census.get("n_census_cells"):
                notes.append(
                    f"Census single-cell context: {census['n_census_cells']:,} cells "
                    f"(us_id={census.get('us_disease_id')})."
                )

        return BiomarkerReport(
            disease_query=disease,
            disease_id=result.bundle.disease_id,
            disease_name=result.bundle.disease_name,
            location_query=str(location or loc.get("name")),
            location=loc,
            top_vocs=top,
            n_vocs_modeled=len(result.bundle.predictions),
            model_version=result.bundle.model_version,
            notes=notes,
            result=result,
            mechanisms=mechanisms,
            census=census,
        )


def _states_for_tissue(tissue: str, category: str | None) -> list[str]:
    t = tissue.lower()
    cat = (category or "").lower()
    if any(x in t for x in ("brain", "nervous", "spinal", "cortex")):
        return ["neuron_stressed", "microglia_activated", "oxidative_stress_cell"]
    if any(x in t for x in ("liver",)):
        return ["hepatocyte_ketogenic", "hepatocyte_sulfur", "oxidative_stress_cell"]
    if any(x in t for x in ("gut", "intestin", "colon", "stomach")):
        return ["gut_fermentative_microbe", "gut_putrefactive_microbe"]
    if any(x in t for x in ("lung", "respirat", "nose")):
        return ["tumor_epithelial_warburg", "oxidative_stress_cell"] if cat == "cancer" else [
            "oxidative_stress_cell"
        ]
    if "adipose" in t or "fat" in t:
        return ["adipocyte_lipolytic"]
    if cat == "cancer":
        return ["tumor_epithelial_warburg", "oxidative_stress_cell"]
    if cat == "metabolic":
        return ["hepatocyte_ketogenic", "adipocyte_lipolytic"]
    if cat == "microbiome":
        return ["gut_fermentative_microbe", "gut_putrefactive_microbe"]
    return ["oxidative_stress_cell", "tumor_epithelial_warburg"]
