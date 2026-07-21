"""Structured query slots for VOC biomarker runs (NL parse or questionnaire)."""

from __future__ import annotations

from typing import Any, Literal, Optional

from pydantic import BaseModel, Field


class QuerySlots(BaseModel):
    """Fields collected from natural language or an interactive questionnaire."""

    disease: Optional[str] = Field(None, description="Primary disease / condition")
    comorbidities: list[str] = Field(default_factory=list)
    location: Optional[str] = Field(
        None, description="Anatomic site of affected cells (e.g. brain, left lower lobe)"
    )
    age_years: Optional[float] = None
    sex: Optional[Literal["female", "male", "other"]] = None
    stage: Optional[str] = None
    genes: list[str] = Field(default_factory=list)
    smoking: Optional[Literal["never", "former", "current"]] = None
    mode: Literal["physiology", "hybrid", "legacy"] = "hybrid"
    top: int = 20
    comorbidity_weight: float = 0.65
    metastatic: bool = False
    source_text: Optional[str] = None
    parse_method: Optional[str] = None  # rules | llm:ollama | llm:openai | ask

    def missing_required(self) -> list[str]:
        miss = []
        if not (self.disease and str(self.disease).strip()):
            miss.append("disease")
        return miss

    def summary_lines(self) -> list[str]:
        lines = [
            f"disease:        {self.disease or '(required)'}",
            f"comorbidities:  {', '.join(self.comorbidities) or '(none)'}",
            f"location:       {self.location or '(default site)'}",
            f"age:            {self.age_years if self.age_years is not None else '(n/a)'}",
            f"sex:            {self.sex or '(n/a)'}",
            f"stage:          {self.stage or '(n/a)'}",
            f"genes:          {', '.join(self.genes) or '(n/a)'}",
            f"smoking:        {self.smoking or '(n/a)'}",
            f"mode:           {self.mode}",
            f"top:            {self.top}",
        ]
        if self.parse_method:
            lines.append(f"parsed via:     {self.parse_method}")
        return lines

    def to_biomarker_kwargs(self) -> dict[str, Any]:
        return {
            "disease": self.disease,
            "location": self.location,
            "top_n": int(self.top),
            "stage": self.stage,
            "genes": list(self.genes) or None,
            "mode": self.mode,
            "sex": self.sex,
            "age_years": self.age_years,
            "smoking_status": self.smoking,
            "metastatic": bool(self.metastatic),
            "comorbidities": list(self.comorbidities) or None,
            "comorbidity_weight": float(self.comorbidity_weight),
        }
