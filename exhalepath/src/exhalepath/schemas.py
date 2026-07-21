from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


class TumorContext(BaseModel):
    """Tumor / lesion covariates that modulate pathway→VOC emission."""

    stage: Optional[str] = Field(
        default=None,
        description="AJCC / clinical stage, e.g. 'I', 'IIA', 'III', 'IV'",
    )
    primary_site: Optional[str] = Field(
        default=None,
        description="Anatomic primary site, e.g. lung, pancreas, liver",
    )
    histology: Optional[str] = Field(
        default=None,
        description="Histologic type, e.g. adenocarcinoma, squamous cell carcinoma",
    )
    tumor_type: Optional[str] = Field(
        default=None,
        description="Broad tumor class: solid, hematologic, metastatic",
    )
    laterality: Optional[str] = None
    grade: Optional[str] = None
    metastatic: bool = False
    tumor_burden_proxy: Optional[float] = Field(
        default=None,
        ge=0.0,
        le=1.0,
        description="Optional 0–1 tumor burden proxy (imaging / clinical)",
    )


class DiseaseQuery(BaseModel):
    """User-facing query for exhaled VOC prediction."""

    disease: str = Field(..., description="Disease name, alias, MONDO id, or TCGA project")
    tumor: Optional[TumorContext] = None
    mutated_genes: list[str] = Field(default_factory=list)
    pathway_overrides: dict[str, float] = Field(
        default_factory=dict,
        description="Optional pathway_id → activity score overrides",
    )
    age_years: Optional[float] = None
    sex: Optional[Literal["female", "male", "other"]] = None
    smoking_status: Optional[Literal["never", "former", "current"]] = None
    include_uncertainty: bool = True


class VOCPrediction(BaseModel):
    voc_id: str
    name: str
    cas: Optional[str] = None
    healthy_ppb: float
    predicted_ppb: float
    delta_ppb: float
    log2_fold_change: float
    fold_change: float
    ci_low_ppb: Optional[float] = None
    ci_high_ppb: Optional[float] = None
    top_pathway_drivers: list[str] = Field(default_factory=list)
    confidence: float = Field(ge=0.0, le=1.0, default=0.5)


class PathwayScore(BaseModel):
    pathway_id: str
    name: str
    score: float
    hit_genes: list[str] = Field(default_factory=list)
    disease_bias: float = 1.0


class PredictionBundle(BaseModel):
    disease_id: str
    disease_name: str
    query: DiseaseQuery
    pathway_scores: list[PathwayScore]
    predictions: list[VOCPrediction]
    model_version: str
    notes: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
