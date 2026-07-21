from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..knowledge.loader import KnowledgeBase
from ..schemas import (
    CellStateActivity,
    DiseaseQuery,
    PathwayScore,
    PhysiologyTrace,
)
from .cell_states import estimate_cell_states
from .chains import score_all_chains
from .transport import integrate_production


@dataclass
class PhysiologyResult:
    cell_states: list[CellStateActivity]
    chain_fluxes: list[dict[str, Any]]
    traces: dict[str, PhysiologyTrace]
    mode: str = "physiology"
    notes: list[str] = field(default_factory=list)

    def ppb_map(self) -> dict[str, float]:
        return {vid: t.predicted_ppb for vid, t in self.traces.items()}


class PhysiologyEngine:
    """
    End-to-end VOC physiology:
      pathway chains × affected cell states → production
      → blood transfer → alveolar release → exhaled ppb
    """

    def __init__(self, kb: KnowledgeBase):
        self.kb = kb

    def run(
        self,
        *,
        disease: dict[str, Any],
        query: DiseaseQuery,
        pathway_scores: list[PathwayScore],
        associated_genes: dict[str, float] | None = None,
    ) -> PhysiologyResult:
        mutated = {g.upper() for g in query.mutated_genes if g}
        cell_states = estimate_cell_states(
            atlas=self.kb.cell_states,
            disease=disease,
            mutated_genes=mutated,
            associated_genes=associated_genes,
            tumor=query.tumor,
            cell_state_fractions=query.cell_state_fractions,
            cell_state_activity=query.cell_state_activity,
        )
        chain_fluxes = score_all_chains(self.kb.pathway_chains, pathway_scores, mutated)

        # Optional ventilation / cardiac output overrides
        physio = self.kb.physio_constants
        if query.alveolar_ventilation_l_per_min or query.cardiac_output_l_per_min:
            physio = {
                **physio,
                "physiology": {
                    **physio.get("physiology", {}),
                    **(
                        {
                            "alveolar_ventilation_l_per_min": query.alveolar_ventilation_l_per_min
                        }
                        if query.alveolar_ventilation_l_per_min
                        else {}
                    ),
                    **(
                        {"cardiac_output_l_per_min": query.cardiac_output_l_per_min}
                        if query.cardiac_output_l_per_min
                        else {}
                    ),
                },
            }

        traces: dict[str, PhysiologyTrace] = {}
        chained_vocs = {c["voc_id"] for c in self.kb.pathway_chains}
        for voc_id, voc in self.kb.vocs.items():
            if voc_id not in chained_vocs:
                continue
            prior = float(disease.get("voc_log2fc_prior", {}).get(voc_id, 0.0))
            traces[voc_id] = integrate_production(
                voc_id=voc_id,
                chain_fluxes=chain_fluxes,
                cell_states=cell_states,
                physio=physio,
                healthy_ppb=float(voc["healthy_ppb_median"]),
                disease_prior_log2fc=prior,
            )

        notes = [
            "Physiology mode: cell-state density×activity × pathway-chain flux → "
            "blood transfer (perfusion, first-pass) → Farhi alveolar release → ppb.",
            "Cell-state densities are atlas priors / optional single-cell fraction overrides, "
            "not raw patient scRNA-seq counts unless you supply cell_state_fractions.",
        ]
        return PhysiologyResult(
            cell_states=cell_states,
            chain_fluxes=chain_fluxes,
            traces=traces,
            notes=notes,
        )

    @staticmethod
    def blend_with_legacy(
        *,
        physio_ppb: float | None,
        legacy_ppb: float,
        healthy_ppb: float,
        weight_physio: float = 0.65,
    ) -> float:
        if physio_ppb is None or not np.isfinite(physio_ppb):
            return legacy_ppb
        w = min(max(weight_physio, 0.0), 1.0)
        return float(w * physio_ppb + (1.0 - w) * legacy_ppb)
