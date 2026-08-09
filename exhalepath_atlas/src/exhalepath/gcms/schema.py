"""BreathVOC research schemas — patient/sample/feature/split manifests.

Designed for GC-MS / PTR diagnostic research enablement (ISA-Tab–inspired),
not as a clinical device schema.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


class SampleRecord(BaseModel):
    sample_id: str
    subject_id: str
    study_id: str
    label: Optional[str] = None  # disease | control | NA
    disease_id: Optional[str] = None
    age: Optional[float] = None
    sex: Optional[str] = None
    smoking_status: Optional[str] = None
    site_id: Optional[str] = None
    batch_id: Optional[str] = None
    collection_device: Optional[str] = None
    breath_fraction: Optional[str] = None
    paired_blank_id: Optional[str] = None
    modality: str = "gcms"
    factors_raw: Optional[str] = None
    metadata: dict[str, Any] = Field(default_factory=dict)


class FeatureAnnotation(BaseModel):
    feature_id: str
    voc_id: Optional[str] = None
    name: Optional[str] = None
    rt: Optional[float] = None
    ri: Optional[float] = None
    msi_level: Optional[Literal["1", "2", "3", "4", "unknown"]] = "unknown"
    inchikey: Optional[str] = None
    library: Optional[str] = None
    match_score: Optional[float] = None
    blank_ratio: Optional[float] = None
    below_lod: Optional[bool] = None


class SplitFold(BaseModel):
    fold_id: str
    train_subject_ids: list[str]
    test_subject_ids: list[str]
    val_subject_ids: list[str] = Field(default_factory=list)


class SplitManifest(BaseModel):
    dataset_id: str
    schema_version: str = "BreathVOC-1.0"
    seed: int
    strategy: str  # nested_kfold | loocv | predefined | stratified_kfold
    n_splits: int
    subject_ids: list[str]
    folds: list[SplitFold]
    content_sha256: str
    created_utc: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )
    notes: list[str] = Field(default_factory=list)
    forbidden_ops_on_test: list[str] = Field(
        default_factory=lambda: [
            "fit_scaler",
            "fit_batch_corrector",
            "feature_selection",
            "imputer_fit",
            "threshold_tuning_on_test",
        ]
    )


class DiagnosticMetrics(BaseModel):
    n_subjects: int
    n_positive: int
    n_negative: int
    auroc: Optional[float] = None
    auroc_ci95: Optional[tuple[float, float]] = None
    auprc: Optional[float] = None
    sensitivity: Optional[float] = None
    specificity: Optional[float] = None
    ppv: Optional[float] = None
    npv: Optional[float] = None
    youden_threshold: Optional[float] = None
    confusion: dict[str, int] = Field(default_factory=dict)
    brier: Optional[float] = None
    nested_auroc: Optional[float] = None
    non_nested_auroc: Optional[float] = None
    optimism_gap: Optional[float] = None
    notes: list[str] = Field(default_factory=list)
