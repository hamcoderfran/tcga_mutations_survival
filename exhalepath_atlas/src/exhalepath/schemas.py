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
    # Optional single-cell informed overrides (fractions 0–1, activity ≥0)
    cell_state_fractions: dict[str, float] = Field(
        default_factory=dict,
        description="Override cell-state densities from scRNA-seq / cytometry fractions",
    )
    cell_state_activity: dict[str, float] = Field(
        default_factory=dict,
        description="Override cell-state metabolic activity scores",
    )
    alveolar_ventilation_l_per_min: Optional[float] = Field(
        default=None, gt=0, description="Override VA (L/min) for alveolar release"
    )
    cardiac_output_l_per_min: Optional[float] = Field(
        default=None, gt=0, description="Override cardiac output Q (L/min)"
    )
    mode: Literal["physiology", "hybrid", "legacy"] = Field(
        default="hybrid",
        description="physiology=cell/blood/alveolar model; legacy=prior path; hybrid=blend",
    )
    age_years: Optional[float] = None
    sex: Optional[Literal["female", "male", "other"]] = None
    smoking_status: Optional[Literal["never", "former", "current"]] = None
    comorbidities: list[str] = Field(
        default_factory=list,
        description="Comorbid conditions fused into pathway/VOC/cell-state priors",
    )
    comorbidity_weight: float = Field(
        default=0.65,
        ge=0.0,
        le=1.5,
        description="Relative weight of each comorbidity prior vs primary disease",
    )
    include_uncertainty: bool = True


class CellStateActivity(BaseModel):
    state_id: str
    name: str
    tissue: str
    density: float
    activity: float
    effective_source: float
    marker_hit_score: float = 0.0
    marker_genes_hit: list[str] = Field(default_factory=list)
    produces_chains: list[str] = Field(default_factory=list)
    disease_modulated: bool = False


class PhysiologyTrace(BaseModel):
    voc_id: str
    tissue_production: float
    blood_delivery: float
    hepatic_first_pass: float
    perfusion_fraction: float
    lambda_blood_air: float
    alveolar_fraction: float
    alveolar_ventilation_l_per_min: float
    cardiac_output_l_per_min: float
    predicted_ppb: float
    contributing_cell_states: list[str] = Field(default_factory=list)
    contributing_chains: list[str] = Field(default_factory=list)
    tissue_breakdown: dict[str, float] = Field(default_factory=dict)


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
    physiology: Optional[PhysiologyTrace] = None


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
    cell_states: list[CellStateActivity] = Field(default_factory=list)
    model_version: str
    notes: list[str] = Field(default_factory=list)
    metadata: dict[str, Any] = Field(default_factory=dict)
