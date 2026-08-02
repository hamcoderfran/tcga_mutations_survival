"""Shared types for the multi-model disease prediction stack."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass
class StackQuery:
    disease: str
    location: Optional[str] = None
    comorbidities: list[str] = field(default_factory=list)
    age: Optional[float] = None
    sex: Optional[str] = None
    genes: list[str] = field(default_factory=list)
    description: Optional[str] = None
    smoking: Optional[str] = None
    top_n: int = 20
    mode: str = "hybrid"


@dataclass
class VOCSignal:
    voc_id: str
    log2fc: float
    confidence: float = 0.5
    rank: Optional[int] = None
    delta_ppb: Optional[float] = None
    evidence: list[str] = field(default_factory=list)


@dataclass
class AspectHit:
    id: str
    name: str
    score: float
    kind: str
    evidence: list[str] = field(default_factory=list)


@dataclass
class ModelOutput:
    model_id: str
    family: str
    aspect: str
    weight: float
    status: str  # ok | degraded | skipped
    voc_signals: list[VOCSignal] = field(default_factory=list)
    aspects: list[AspectHit] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


@dataclass
class FusedVOC:
    voc_id: str
    name: str
    fused_log2fc: float
    fused_delta_ppb: float
    fused_confidence: float
    rrf_score: float
    n_models_agreeing: int
    model_votes: dict[str, float]
    model_confidences: dict[str, float]
    evidence: list[str] = field(default_factory=list)
    epistemic_std: float = 0.0
    ci_low_log2fc: float = 0.0
    ci_high_log2fc: float = 0.0


@dataclass
class StackResult:
    query: StackQuery
    disease_id: str
    disease_name: str
    location: dict[str, Any]
    model_outputs: list[ModelOutput]
    fused_vocs: list[FusedVOC]
    fused_aspects: dict[str, list[AspectHit]]
    fusion_weights: dict[str, float]
    summary: dict[str, Any]
    notes: list[str] = field(default_factory=list)
    biomarker_report: Any = None
