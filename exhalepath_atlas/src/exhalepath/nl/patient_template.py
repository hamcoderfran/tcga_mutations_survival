"""
Naturalistic patient input → canonical PatientTemplate.

Accepts diverse free text (vignettes, clinic notes, bullet lists, JSON/YAML-ish
key=value blocks) and normalizes into one schema the stack / biomarker engines
understand.
"""

from __future__ import annotations

import json
import re
from typing import Any, Literal, Optional

from pydantic import BaseModel, Field

from .rules import parse_rules
from .slots import QuerySlots


class PatientTemplate(BaseModel):
    """Canonical patient profile consumed by ExhalePath / Great Disease Stack."""

    # --- demographics ---
    age_years: Optional[float] = None
    sex: Optional[Literal["female", "male", "other"]] = None
    height_cm: Optional[float] = None
    weight_kg: Optional[float] = None
    bmi: Optional[float] = None

    # --- clinical ---
    primary_disease: Optional[str] = None
    comorbidities: list[str] = Field(default_factory=list)
    symptoms: list[str] = Field(default_factory=list)
    medications: list[str] = Field(default_factory=list)
    allergies: list[str] = Field(default_factory=list)
    chief_complaint: Optional[str] = None
    history: Optional[str] = None

    # --- biology / anatomy ---
    location: Optional[str] = None
    genes: list[str] = Field(default_factory=list)
    stage: Optional[str] = None
    metastatic: bool = False
    description: Optional[str] = None  # phenotype / mechanism narrative for zero-shot

    # --- lifestyle ---
    smoking: Optional[Literal["never", "former", "current"]] = None
    alcohol: Optional[str] = None  # none | social | heavy | unknown
    pregnant: Optional[bool] = None

    # --- run controls ---
    mode: Literal["physiology", "hybrid", "legacy"] = "hybrid"
    top: int = 20
    comorbidity_weight: float = 0.65

    # --- provenance ---
    source_text: Optional[str] = None
    parse_method: Optional[str] = None
    parse_confidence: float = 0.5
    unresolved: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    field_sources: dict[str, str] = Field(default_factory=dict)

    def missing_required(self) -> list[str]:
        return ["primary_disease"] if not (self.primary_disease or "").strip() else []

    def summary_lines(self) -> list[str]:
        lines = [
            f"primary_disease:  {self.primary_disease or '(required)'}",
            f"comorbidities:    {', '.join(self.comorbidities) or '(none)'}",
            f"location:         {self.location or '(default site)'}",
            f"age / sex:        {self.age_years if self.age_years is not None else 'n/a'} / {self.sex or 'n/a'}",
            f"genes:            {', '.join(self.genes) or '(n/a)'}",
            f"stage:            {self.stage or '(n/a)'}",
            f"smoking:          {self.smoking or '(n/a)'}",
            f"medications:      {', '.join(self.medications) or '(n/a)'}",
            f"symptoms:         {', '.join(self.symptoms[:8]) or '(n/a)'}",
            f"bmi:              {self.bmi if self.bmi is not None else '(n/a)'}",
            f"description:      {(self.description or '')[:120] or '(n/a)'}",
        ]
        if self.parse_method:
            lines.append(
                f"parsed via:       {self.parse_method} (conf={self.parse_confidence:.2f})"
            )
        if self.warnings:
            lines.append(f"warnings:         {'; '.join(self.warnings[:4])}")
        return lines

    def to_query_slots(self) -> QuerySlots:
        return QuerySlots(
            disease=self.primary_disease,
            comorbidities=list(self.comorbidities),
            location=self.location,
            age_years=self.age_years,
            sex=self.sex,
            stage=self.stage,
            genes=list(self.genes),
            smoking=self.smoking,
            mode=self.mode,
            top=self.top,
            comorbidity_weight=self.comorbidity_weight,
            metastatic=self.metastatic,
            source_text=self.source_text,
            parse_method=self.parse_method,
        )

    def to_stack_kwargs(self) -> dict[str, Any]:
        # Fold meds/symptoms/history into description for zero-shot mechanism heads
        bits = []
        if self.description:
            bits.append(self.description)
        if self.chief_complaint:
            bits.append(f"CC: {self.chief_complaint}")
        if self.history:
            bits.append(f"HPI: {self.history}")
        if self.symptoms:
            bits.append("symptoms: " + ", ".join(self.symptoms))
        if self.medications:
            bits.append("medications: " + ", ".join(self.medications))
        if self.allergies:
            bits.append("allergies: " + ", ".join(self.allergies))
        if self.alcohol:
            bits.append(f"alcohol: {self.alcohol}")
        if self.bmi is not None:
            bits.append(f"BMI {self.bmi:.1f}")
        if self.pregnant:
            bits.append("pregnant")
        return {
            "disease": self.primary_disease,
            "location": self.location,
            "comorbidities": list(self.comorbidities),
            "age": self.age_years,
            "sex": self.sex,
            "genes": list(self.genes),
            "description": "; ".join(bits) if bits else None,
            "smoking": self.smoking,
            "top_n": int(self.top),
        }

    def to_biomarker_kwargs(self) -> dict[str, Any]:
        kw = self.to_query_slots().to_biomarker_kwargs()
        # biomarker engine accepts description via predict()
        desc = self.to_stack_kwargs().get("description")
        if desc:
            kw["description"] = desc
        return kw


# ---------------------------------------------------------------------------
# Diverse input normalization
# ---------------------------------------------------------------------------

_KV_LINE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9_ /\-]{0,40})\s*[:=]\s*(.+?)\s*$"
)

_SYMPTOM_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bhallucin", re.I), "hallucinations"),
    (re.compile(r"\bdelusion", re.I), "delusions"),
    (re.compile(r"\bparanoia\b|\bparanoid\b", re.I), "paranoia"),
    (re.compile(r"\banhedonia\b", re.I), "anhedonia"),
    (re.compile(r"\binsomnia\b|\bcan'?t sleep\b", re.I), "insomnia"),
    (re.compile(r"\bfatigue\b|\btired\b|\blow energy\b", re.I), "fatigue"),
    (re.compile(r"\bdyspnea\b|\bshort(?:ness)? of breath\b|\bSOB\b", re.I), "dyspnea"),
    (re.compile(r"\bcough\b", re.I), "cough"),
    (re.compile(r"\bchest pain\b", re.I), "chest_pain"),
    (re.compile(r"\bweight loss\b", re.I), "weight_loss"),
    (re.compile(r"\bweight gain\b", re.I), "weight_gain"),
    (re.compile(r"\bpolyuria\b|\bpolydipsia\b", re.I), "polyuria_polydipsia"),
    (re.compile(r"\btremor\b", re.I), "tremor"),
    (re.compile(r"\bmemory loss\b|\bcognitive decline\b", re.I), "cognitive_decline"),
    (re.compile(r"\banxiety\b|\bpanic\b", re.I), "anxiety"),
    (re.compile(r"\bsuicid", re.I), "suicidal_ideation"),
]

_MED_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bolanzapine\b|\bzyprexa\b", re.I), "olanzapine"),
    (re.compile(r"\brisperidone\b|\brisperdal\b", re.I), "risperidone"),
    (re.compile(r"\bquetiapine\b|\bseroquel\b", re.I), "quetiapine"),
    (re.compile(r"\baripiprazole\b|\babilify\b", re.I), "aripiprazole"),
    (re.compile(r"\bclozapine\b", re.I), "clozapine"),
    (re.compile(r"\blithium\b", re.I), "lithium"),
    (re.compile(r"\bsertraline\b|\bzoloft\b", re.I), "sertraline"),
    (re.compile(r"\bfluoxetine\b|\bprozac\b", re.I), "fluoxetine"),
    (re.compile(r"\bescitalopram\b|\blexapro\b", re.I), "escitalopram"),
    (re.compile(r"\bmetformin\b", re.I), "metformin"),
    (re.compile(r"\binsulin\b", re.I), "insulin"),
    (re.compile(r"\blisinopril\b|\benalapril\b|\bACE.?i", re.I), "ace_inhibitor"),
    (re.compile(r"\batorvastatin\b|\bsimvastatin\b|\bstatin", re.I), "statin"),
    (re.compile(r"\balbuterol\b|\bsalbutamol\b", re.I), "albuterol"),
    (re.compile(r"\bwarfarin\b|\bapixaban\b|\beliquis\b", re.I), "anticoagulant"),
    (re.compile(r"\bprednisone\b|\bsteroids?\b", re.I), "corticosteroid"),
    (re.compile(r"\bmetoprolol\b|\bbeta.?block", re.I), "beta_blocker"),
]

_EXTRA_COMORBID: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bOSA\b|\bsleep apnea\b", re.I), "osa"),
    (re.compile(r"\bNAFLD\b|\bfatty liver\b", re.I), "nafld"),
    (re.compile(r"\bhyperlipid\b|\bhigh cholesterol\b|\bdyslipid", re.I), "dyslipidemia"),
    (re.compile(r"\bhypothyroid", re.I), "hypothyroidism"),
    (re.compile(r"\banemia\b", re.I), "anemia"),
    (re.compile(r"\bPTSD\b", re.I), "ptsd"),
    (re.compile(r"\bbipolar\b", re.I), "bipolar_disorder"),
    (re.compile(r"\banxiety disorder\b|\bGAD\b", re.I), "anxiety_disorder"),
    (re.compile(r"\brheumatoid\b|\bRA\b", re.I), "rheumatoid_arthritis"),
    (re.compile(r"\bIBD\b|\bCrohn|\bulcerative colitis\b", re.I), "ibd"),
]


def _parse_structured_block(text: str) -> dict[str, Any]:
    """Parse JSON object or key: value / key=value lines."""
    raw = (text or "").strip()
    if not raw:
        return {}
    # JSON
    if raw[0] in "{[":
        try:
            obj = json.loads(raw)
            if isinstance(obj, dict):
                return obj
            if isinstance(obj, list) and obj and isinstance(obj[0], dict):
                return obj[0]
        except json.JSONDecodeError:
            pass
    # key-value lines
    out: dict[str, Any] = {}
    for line in raw.splitlines():
        m = _KV_LINE.match(line)
        if not m:
            continue
        key = re.sub(r"[^a-z0-9]+", "_", m.group(1).strip().lower()).strip("_")
        val = m.group(2).strip().strip("\"'")
        if key in {
            "comorbidities",
            "comorbidity",
            "pmh",
            "genes",
            "medications",
            "meds",
            "symptoms",
            "allergies",
        }:
            out[key] = [x.strip() for x in re.split(r"[,;|/]", val) if x.strip()]
        else:
            out[key] = val
    return out


def _apply_structured(tpl: PatientTemplate, data: dict[str, Any]) -> None:
    alias = {
        "disease": "primary_disease",
        "dx": "primary_disease",
        "diagnosis": "primary_disease",
        "condition": "primary_disease",
        "primary": "primary_disease",
        "age": "age_years",
        "age_years": "age_years",
        "sex": "sex",
        "gender": "sex",
        "site": "location",
        "location": "location",
        "tissue": "location",
        "organ": "location",
        "stage": "stage",
        "smoking": "smoking",
        "smoker": "smoking",
        "bmi": "bmi",
        "weight": "weight_kg",
        "weight_kg": "weight_kg",
        "height": "height_cm",
        "height_cm": "height_cm",
        "cc": "chief_complaint",
        "chief_complaint": "chief_complaint",
        "hpi": "history",
        "history": "history",
        "description": "description",
        "phenotype": "description",
        "notes": "description",
        "alcohol": "alcohol",
        "mode": "mode",
        "top": "top",
    }
    list_alias = {
        "comorbidities": "comorbidities",
        "comorbidity": "comorbidities",
        "pmh": "comorbidities",
        "genes": "genes",
        "mutations": "genes",
        "medications": "medications",
        "meds": "medications",
        "drugs": "medications",
        "symptoms": "symptoms",
        "sx": "symptoms",
        "allergies": "allergies",
    }
    for k, v in data.items():
        key = re.sub(r"[^a-z0-9]+", "_", str(k).lower()).strip("_")
        if key in list_alias:
            field = list_alias[key]
            items = v if isinstance(v, list) else re.split(r"[,;|/]", str(v))
            cleaned = [str(x).strip() for x in items if str(x).strip()]
            setattr(tpl, field, list(dict.fromkeys((getattr(tpl, field) or []) + cleaned)))
            tpl.field_sources[field] = "structured"
            continue
        field = alias.get(key)
        if not field:
            continue
        if field == "age_years":
            try:
                tpl.age_years = float(re.search(r"[\d.]+", str(v)).group(0))  # type: ignore[union-attr]
                tpl.field_sources["age_years"] = "structured"
            except Exception:  # noqa: BLE001
                pass
        elif field == "sex":
            s = str(v).strip().lower()
            if s.startswith("f"):
                tpl.sex = "female"
            elif s.startswith("m") and "female" not in s:
                tpl.sex = "male"
            else:
                tpl.sex = "other"
            tpl.field_sources["sex"] = "structured"
        elif field == "smoking":
            s = str(v).strip().lower()
            if "never" in s or "non" in s:
                tpl.smoking = "never"
            elif "former" in s or "ex" in s or "quit" in s:
                tpl.smoking = "former"
            elif "current" in s or s in {"yes", "y", "smoker"}:
                tpl.smoking = "current"
            tpl.field_sources["smoking"] = "structured"
        elif field in {"bmi", "weight_kg", "height_cm"}:
            try:
                setattr(tpl, field, float(re.search(r"[\d.]+", str(v)).group(0)))  # type: ignore[union-attr]
                tpl.field_sources[field] = "structured"
            except Exception:  # noqa: BLE001
                pass
        elif field == "top":
            try:
                tpl.top = int(float(v))
            except Exception:  # noqa: BLE001
                pass
        elif field == "mode" and str(v) in {"physiology", "hybrid", "legacy"}:
            tpl.mode = str(v)  # type: ignore[assignment]
        elif field == "metastatic":
            tpl.metastatic = str(v).lower() in {"1", "true", "yes", "y"}
        else:
            setattr(tpl, field, str(v).strip())
            tpl.field_sources[field] = "structured"


def _extract_compact_demographics(text: str, tpl: PatientTemplate) -> None:
    """Handle forms like 35M, 35 F, 24yo F, Mr./Ms."""
    # strip common emoji / pictographs that break token boundaries
    text = re.sub(r"[\U0001F300-\U0001FAFF]", " ", text)
    # 35M / 35F / 35 M / 28F w/
    m = re.search(r"\b(\d{1,3})\s*([MFmf])\b(?!\w)", text)
    if m and tpl.age_years is None:
        tpl.age_years = float(m.group(1))
        tpl.sex = "male" if m.group(2).upper() == "M" else "female"
        tpl.field_sources["age_years"] = "compact"
        tpl.field_sources["sex"] = "compact"
    # BMI
    m = re.search(r"\bBMI\s*[:=]?\s*(\d{1,2}(?:\.\d+)?)\b", text, re.I)
    if m:
        tpl.bmi = float(m.group(1))
        tpl.field_sources["bmi"] = "regex"
    # weight / height
    m = re.search(r"\b(\d{2,3}(?:\.\d+)?)\s*kg\b", text, re.I)
    if m:
        tpl.weight_kg = float(m.group(1))
    m = re.search(r"\b(\d{2,3}(?:\.\d+)?)\s*cm\b", text, re.I)
    if m:
        tpl.height_cm = float(m.group(1))
    if tpl.bmi is None and tpl.weight_kg and tpl.height_cm and tpl.height_cm > 0:
        h_m = tpl.height_cm / 100.0
        tpl.bmi = tpl.weight_kg / (h_m * h_m)
        tpl.field_sources["bmi"] = "computed"
    # pregnant
    if re.search(r"\bpregnan", text, re.I):
        tpl.pregnant = True
    # alcohol
    if re.search(r"\bheavy alcohol|\balcohol(?:ism| use disorder)\b|\bEtOH abuse\b", text, re.I):
        tpl.alcohol = "heavy"
    elif re.search(r"\bsocial drinker\b|\boccasional alcohol\b", text, re.I):
        tpl.alcohol = "social"
    elif re.search(r"\bno alcohol\b|\bdoes not drink\b|\bteetotal", text, re.I):
        tpl.alcohol = "none"


def _extract_note_sections(text: str, tpl: PatientTemplate) -> None:
    """Pull CC:/HPI:/PMH:/Meds: clinic-note sections when present."""
    patterns = {
        "chief_complaint": r"(?:^|\n)\s*(?:CC|Chief Complaint)\s*[:\-]\s*(.+?)(?=\n\s*[A-Z][A-Za-z ]{1,20}:|\Z)",
        "history": r"(?:^|\n)\s*(?:HPI|History(?: of Present Illness)?)\s*[:\-]\s*(.+?)(?=\n\s*[A-Z][A-Za-z ]{1,20}:|\Z)",
        "medications_block": r"(?:^|\n)\s*(?:Meds|Medications|Current Medications)\s*[:\-]\s*(.+?)(?=\n\s*[A-Z][A-Za-z ]{1,20}:|\Z)",
        "pmh_block": r"(?:^|\n)\s*(?:PMH|Past Medical History|History)\s*[:\-]\s*(.+?)(?=\n\s*[A-Z][A-Za-z ]{1,20}:|\Z)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, re.I | re.S)
        if not m:
            continue
        val = re.sub(r"\s+", " ", m.group(1)).strip()
        if key == "chief_complaint":
            tpl.chief_complaint = val
            tpl.field_sources["chief_complaint"] = "note_section"
        elif key == "history":
            tpl.history = val
            tpl.field_sources["history"] = "note_section"
        elif key == "medications_block":
            meds = [x.strip() for x in re.split(r"[,;/\n]", val) if x.strip()]
            tpl.medications = list(dict.fromkeys(tpl.medications + meds))
            tpl.field_sources["medications"] = "note_section"
        elif key == "pmh_block":
            items = [x.strip() for x in re.split(r"[,;/\n]", val) if x.strip()]
            tpl.comorbidities = list(dict.fromkeys(tpl.comorbidities + items))
            tpl.field_sources["comorbidities"] = "note_section"


def parse_patient_template(
    text: str,
    *,
    disease_catalog: list[dict[str, Any]] | None = None,
    llm: str = "rules",
) -> PatientTemplate:
    """
    Convert naturalistic patient input into a PatientTemplate.

    Accepts:
      - free-text vignettes ("35M with schizophrenia, smokes, on olanzapine")
      - clinic-note sections (CC:/HPI:/PMH:/Meds:)
      - key: value or key=value lines
      - JSON objects with the same fields
    """
    raw = (text or "").strip()
    tpl = PatientTemplate(source_text=raw, parse_method="patient_template")
    if not raw:
        tpl.warnings.append("empty_input")
        tpl.parse_confidence = 0.0
        return tpl

    # 1) structured overlay (JSON / kv)
    structured = _parse_structured_block(raw)
    if structured:
        _apply_structured(tpl, structured)
        tpl.parse_method = "patient_template+structured"

    # 2) clinic note sections
    _extract_note_sections(raw, tpl)

    # 3) compact demographics + BMI
    _extract_compact_demographics(raw, tpl)

    # 4) reuse battle-tested rule slots for disease/location/genes/smoking/age/sex
    slots = parse_rules(raw, disease_catalog=disease_catalog)
    if slots.disease and not tpl.primary_disease:
        tpl.primary_disease = slots.disease
        tpl.field_sources["primary_disease"] = "rules"
    if slots.location and not tpl.location:
        tpl.location = slots.location
        tpl.field_sources["location"] = "rules"
    if slots.age_years is not None and tpl.age_years is None:
        tpl.age_years = slots.age_years
        tpl.field_sources["age_years"] = "rules"
    if slots.sex and not tpl.sex:
        tpl.sex = slots.sex
        tpl.field_sources["sex"] = "rules"
    if slots.smoking and not tpl.smoking:
        tpl.smoking = slots.smoking
        tpl.field_sources["smoking"] = "rules"
    if slots.stage and not tpl.stage:
        tpl.stage = slots.stage
    if slots.genes:
        tpl.genes = list(dict.fromkeys(tpl.genes + slots.genes))
    if slots.comorbidities:
        tpl.comorbidities = list(dict.fromkeys(tpl.comorbidities + slots.comorbidities))
    if slots.metastatic:
        tpl.metastatic = True

    # 5) extra comorbidities / symptoms / meds from open patterns
    for pat, name in _EXTRA_COMORBID:
        if pat.search(raw):
            tpl.comorbidities.append(name)
    tpl.comorbidities = list(dict.fromkeys(tpl.comorbidities))

    for pat, name in _SYMPTOM_PATTERNS:
        if pat.search(raw):
            tpl.symptoms.append(name)
    tpl.symptoms = list(dict.fromkeys(tpl.symptoms))

    for pat, name in _MED_PATTERNS:
        if pat.search(raw):
            tpl.medications.append(name)
    # also "on X" / "taking X"
    for m in re.finditer(
        r"\b(?:on|taking|prescribed|rx)\s+([A-Za-z][A-Za-z0-9\-]{2,30})\b", raw, re.I
    ):
        drug = m.group(1).lower()
        if drug not in {"the", "and", "with", "home", "oxygen", "therapy"}:
            tpl.medications.append(drug)
    tpl.medications = list(dict.fromkeys(tpl.medications))

    # BMI-implied obesity comorbidity
    if tpl.bmi is not None and tpl.bmi >= 30 and "obesity" not in [
        c.lower() for c in tpl.comorbidities
    ]:
        tpl.comorbidities.append("obesity")
        tpl.field_sources["comorbidities"] = tpl.field_sources.get(
            "comorbidities", "bmi_inferred"
        )

    # 6) optional LLM enrichment when requested (fills gaps only)
    if llm and llm != "rules":
        try:
            from .llm import parse_with_optional_llm

            llm_slots = parse_with_optional_llm(
                raw, llm=llm, disease_catalog=disease_catalog
            )
            if llm_slots.disease and not tpl.primary_disease:
                tpl.primary_disease = llm_slots.disease
                tpl.field_sources["primary_disease"] = llm_slots.parse_method or "llm"
            if llm_slots.location and not tpl.location:
                tpl.location = llm_slots.location
            if llm_slots.age_years is not None and tpl.age_years is None:
                tpl.age_years = llm_slots.age_years
            if llm_slots.sex and not tpl.sex:
                tpl.sex = llm_slots.sex
            if llm_slots.comorbidities:
                tpl.comorbidities = list(
                    dict.fromkeys(tpl.comorbidities + llm_slots.comorbidities)
                )
            if llm_slots.genes:
                tpl.genes = list(dict.fromkeys(tpl.genes + llm_slots.genes))
            tpl.parse_method = f"patient_template+{llm_slots.parse_method or llm}"
        except Exception as exc:  # noqa: BLE001
            tpl.warnings.append(f"llm_enrichment_failed:{exc}")

    # 7) cautious fallback — only from CC/HPI, never from bare gene/demo fragments
    if not tpl.primary_disease:
        for cand in (tpl.chief_complaint, tpl.history):
            if not cand:
                continue
            frag = re.sub(
                r"^\s*\d{1,3}\s*[yoMFmf\- ]*(?:year\s*old)?\s*", "", cand, flags=re.I
            )
            frag = re.sub(
                r"^\s*(?:a|an|the|with|presents with|c/o|complains of)\s+",
                "",
                frag,
                flags=re.I,
            ).strip(" .;")
            if re.search(
                r"\b(mutation|gene|BMI|smoker|year old|taking|metformin|insulin)\b",
                frag,
                re.I,
            ):
                continue
            if 3 <= len(frag) <= 80:
                tpl.primary_disease = frag
                tpl.field_sources["primary_disease"] = "fallback_phrase"
                tpl.warnings.append("primary_disease_from_fallback_phrase")
                break

    # strip primary from comorbidities
    if tpl.primary_disease:
        pl = tpl.primary_disease.lower().replace(" ", "_")
        tpl.comorbidities = [
            c
            for c in tpl.comorbidities
            if c.lower() not in {pl, tpl.primary_disease.lower()}
            and c.lower().replace("_", " ") != tpl.primary_disease.lower()
        ]

    # confidence heuristic
    score = 0.15
    if tpl.primary_disease:
        score += 0.35
    if tpl.location:
        score += 0.1
    if tpl.age_years is not None:
        score += 0.1
    if tpl.sex:
        score += 0.05
    if tpl.comorbidities:
        score += 0.05
    if tpl.genes:
        score += 0.05
    if tpl.medications or tpl.symptoms:
        score += 0.05
    if "fallback_phrase" in (tpl.field_sources.get("primary_disease") or ""):
        score -= 0.15
    tpl.parse_confidence = max(0.05, min(0.99, score))

    if not tpl.primary_disease:
        tpl.unresolved.append("primary_disease")
    return tpl
