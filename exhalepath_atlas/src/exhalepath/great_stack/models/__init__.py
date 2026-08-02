"""Register all stack model adapters."""

from __future__ import annotations

from typing import Any

from ..base import StackModel
from .core_exhalepath import (
    CalibratorModel,
    HybridExhalePathModel,
    PhysiologyExhalePathModel,
    ZeroShotModel,
)
from .cutting_edge import CounterfactualNullModel, MetaEnsembleModel, PhenotypeMondoModel
from .graph_genetics import OmniPathSignalingModel, OpenTargetsGenesModel, PrimeKGGraphModel
from .literature_clinical import (
    CellCensusModel,
    ChemblPharmModel,
    ComorbidityModel,
    LiteraturePriorModel,
    PathwayEnrichmentModel,
)
from .systems import AgoraMicrobiomeModel, HumanGEMFluxModel, OperaPhyschemModel, PBPKTransportModel
from .zero_shot_evidence import ZeroShotEvidenceModel


def build_all_models(ctx: dict[str, Any]) -> list[StackModel]:
    """Instantiate the full cutting-edge multi-model stack (20 heads)."""
    classes: list[type[StackModel]] = [
        HybridExhalePathModel,
        PhysiologyExhalePathModel,
        CalibratorModel,
        ZeroShotModel,
        ZeroShotEvidenceModel,
        LiteraturePriorModel,
        OpenTargetsGenesModel,
        HumanGEMFluxModel,
        OperaPhyschemModel,
        PrimeKGGraphModel,
        OmniPathSignalingModel,
        AgoraMicrobiomeModel,
        PBPKTransportModel,
        CellCensusModel,
        ComorbidityModel,
        ChemblPharmModel,
        PathwayEnrichmentModel,
        PhenotypeMondoModel,
        CounterfactualNullModel,
        MetaEnsembleModel,
    ]
    # apply registry weights if present
    reg = {m["id"]: m for m in (ctx.get("registry") or {}).get("models") or []}
    models: list[StackModel] = []
    for cls in classes:
        m = cls(ctx)
        meta = reg.get(m.model_id) or {}
        if "weight" in meta:
            m.default_weight = float(meta["weight"])
        models.append(m)
    return models


__all__ = ["build_all_models"]
