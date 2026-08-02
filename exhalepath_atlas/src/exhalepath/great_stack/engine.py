"""Orchestrate all stack models and fuse into one comprehensive prediction."""

from __future__ import annotations

from typing import Any, Optional

from ..knowledge.loader import KnowledgeBase, default_knowledge
from .data import model_registry
from .fusion import (
    calibrated_weight_overrides,
    detect_zero_shot,
    effective_weight,
    fuse_aspects,
    fuse_vocs,
    fusion_summary,
)
from .models import build_all_models
from .types import StackQuery, StackResult


class GreatDiseaseStack:
    """20-model cutting-edge fused disease × VOC × mechanism predictor."""

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
        # fresh per-query context fields
        self.ctx["biomarker_report"] = None
        self.ctx["disease_resolved"] = None

        models = build_all_models(self.ctx)
        outputs = []
        # hybrid first (ctx), meta last (needs hybrid report)
        priority = {
            "exhalepath_hybrid": 0,
            "exhalepath_physiology": 1,
            "exhalepath_calibrator": 2,
            "meta_ensemble": 90,
        }
        order = {m.model_id: i for i, m in enumerate(models)}
        models.sort(
            key=lambda m: (
                priority.get(m.model_id, 10),
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

        zero_shot = detect_zero_shot(outputs, query)
        fused_vocs = fuse_vocs(
            outputs,
            voc_catalog=self.kb.vocs,
            top_n=query.top_n,
            query=query,
            zero_shot=zero_shot,
        )
        aspects = fuse_aspects(outputs)
        overrides = calibrated_weight_overrides()
        weights = {
            m.model_id: effective_weight(m, zero_shot=zero_shot, overrides=overrides)
            for m in outputs
        }
        summary = fusion_summary(query, outputs, fused_vocs, zero_shot=zero_shot)
        notes = [
            "Great Disease Prediction Stack — cutting-edge multi-model fusion "
            "(adaptive anti-dilution + epistemic UQ + zero-shot VOC projection).",
            "Research / hypothesis-generation only — not a medical device.",
        ]
        if zero_shot:
            notes.append("Zero-shot regime: mechanism/theme channels upweighted.")
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
