"""Orchestrate all stack models and fuse into one comprehensive prediction."""

from __future__ import annotations

from typing import Any, Optional

from ..knowledge.loader import KnowledgeBase, default_knowledge
from .data import model_registry
from .fusion import fuse_aspects, fuse_vocs, fusion_summary
from .models import build_all_models
from .types import StackQuery, StackResult


class GreatDiseaseStack:
    """16-model fused disease × VOC × mechanism predictor."""

    def __init__(self, kb: KnowledgeBase | None = None):
        self.kb = kb or default_knowledge()
        self.registry = model_registry()
        self.ctx: dict[str, Any] = {
            "kb": self.kb,
            "registry": self.registry,
            "engine": None,
            "biomarker_report": None,
        }

    def predict(self, query: StackQuery | dict[str, Any]) -> StackResult:
        if isinstance(query, dict):
            query = StackQuery(**query)
        models = build_all_models(self.ctx)
        outputs = []
        # run hybrid first so others can reuse biomarker_report
        order = {m.model_id: i for i, m in enumerate(models)}
        models.sort(
            key=lambda m: (
                0 if m.model_id == "exhalepath_hybrid" else 1,
                order.get(m.model_id, 99),
            )
        )
        for model in models:
            outputs.append(model.predict(query))

        resolved = self.ctx.get("disease_resolved") or {}
        if not resolved:
            d = self.kb.resolve_disease(
                query.disease, location_hint=query.location, description=query.description
            )
            resolved = {
                "disease_id": d.get("disease_id"),
                "disease_name": d.get("name"),
                "location": {"name": query.location or d.get("default_site")},
            }

        fused_vocs = fuse_vocs(outputs, voc_catalog=self.kb.vocs, top_n=query.top_n)
        aspects = fuse_aspects(outputs)
        weights = {m.model_id: m.weight for m in outputs}
        summary = fusion_summary(query, outputs, fused_vocs)
        notes = [
            "Great Disease Prediction Stack — multi-model fusion across VOC + disease biology.",
            "Research / hypothesis-generation only — not a medical device.",
        ]
        for mo in outputs:
            notes.extend(mo.notes[:1])

        return StackResult(
            query=query,
            disease_id=str(resolved.get("disease_id") or query.disease),
            disease_name=str(resolved.get("disease_name") or query.disease),
            location=dict(resolved.get("location") or {"name": query.location}),
            model_outputs=outputs,
            fused_vocs=fused_vocs,
            fused_aspects=aspects,
            fusion_weights=weights,
            summary=summary,
            notes=notes,
            biomarker_report=self.ctx.get("biomarker_report"),
        )


def run_great_stack(
    disease: str,
    *,
    location: Optional[str] = None,
    comorbidities: Optional[list[str]] = None,
    age: Optional[float] = None,
    sex: Optional[str] = None,
    genes: Optional[list[str]] = None,
    description: Optional[str] = None,
    smoking: Optional[str] = None,
    top_n: int = 20,
) -> StackResult:
    stack = GreatDiseaseStack()
    return stack.predict(
        StackQuery(
            disease=disease,
            location=location,
            comorbidities=list(comorbidities or []),
            age=age,
            sex=sex,
            genes=list(genes or []),
            description=description,
            smoking=smoking,
            top_n=top_n,
        )
    )
