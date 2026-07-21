"""Very-light natural-language → QuerySlots parser (rules; no model weights)."""

from __future__ import annotations

import re
from typing import Any

from .slots import QuerySlots

# Comorbidity / phenotype phrases → atlas disease ids / names
_COMORBID_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bobese\b|\bobesity\b|\bBMI\b|\boverweight\b", re.I), "obesity"),
    (re.compile(r"\bheart condition\b|\bheart disease\b|\bcardiac\b|\bCVD\b|\bCHF\b|\bheart failure\b|\bcoronary\b", re.I), "heart_disease"),
    (re.compile(r"\btype\s*2\s*diabetes\b|\bT2D\b|\bdiabet(?:es|ic)\b", re.I), "type_2_diabetes"),
    (re.compile(r"\bCOPD\b|\bemphysema\b|\bchronic obstructive\b", re.I), "copd"),
    (re.compile(r"\basthma\b", re.I), "asthma"),
    (re.compile(r"\bcirrhosis\b|\bliver failure\b", re.I), "cirrhosis"),
    (re.compile(r"\bCKD\b|\bchronic kidney\b|\brenal failure\b", re.I), "chronic_kidney_disease"),
    (re.compile(r"\bhypertension\b|\bhigh blood pressure\b", re.I), "hypertension"),
]

# Primary disease phrase hints (order matters — more specific first)
_DISEASE_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\blung adenocarcinoma\b|\bLUAD\b|\bNSCLC\b", re.I), "lung adenocarcinoma"),
    (re.compile(r"\bsmall cell lung\b|\bSCLC\b", re.I), "lung adenocarcinoma"),
    (re.compile(r"\blung (?:cancer|tumor|tumour|carcinoma)\b", re.I), "lung adenocarcinoma"),
    (re.compile(r"\bbreast cancer\b", re.I), "breast cancer"),
    (re.compile(r"\bpancreatic\b", re.I), "pancreatic adenocarcinoma"),
    (re.compile(r"\bcolon cancer\b|\bcolorectal\b", re.I), "colon adenocarcinoma"),
    (re.compile(r"\bmajor depressive\b|\bdepression\b|\bMDD\b|\bdepressive disorder\b", re.I), "depression"),
    (re.compile(r"\bschizophren", re.I), "schizophrenia"),
    (re.compile(r"\bAlzheimer", re.I), "alzheimer disease"),
    (re.compile(r"\bParkinson", re.I), "parkinson disease"),
    (re.compile(r"\btype\s*2\s*diabetes\b|\bT2D\b", re.I), "type 2 diabetes"),
    (re.compile(r"\bcirrhosis\b", re.I), "cirrhosis"),
    (re.compile(r"\basthma\b", re.I), "asthma"),
    (re.compile(r"\bCOPD\b", re.I), "copd"),
    (re.compile(r"\bobesity\b", re.I), "obesity"),
    (re.compile(r"\bheart (?:disease|failure|condition)\b", re.I), "heart disease"),
]

_LOCATION_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"\bleft lower lobe\b|\bLLL\b", re.I), "left lower lobe"),
    (re.compile(r"\bright lower lobe\b|\bRLL\b", re.I), "right lower lobe"),
    (re.compile(r"\bleft upper lobe\b|\bLUL\b", re.I), "left upper lobe"),
    (re.compile(r"\bright upper lobe\b|\bRUL\b", re.I), "right upper lobe"),
    (re.compile(r"\bin the lung\b|\bof the lung\b|\bpulmonary\b", re.I), "lung"),
    (re.compile(r"\bin the brain\b|\bcerebral\b|\bneural\b", re.I), "brain"),
    (re.compile(r"\bin the (?:liver|hepatic)\b", re.I), "liver"),
    (re.compile(r"\bin the (?:heart|cardiac)\b", re.I), "heart"),
    (re.compile(r"\bin the (?:pancreas|pancreatic)\b", re.I), "pancreas"),
    (re.compile(r"\bin the (?:colon|gut|intestine)\b", re.I), "colon"),
    (re.compile(r"\bsystemic\b|\bwhole body\b", re.I), "systemic"),
]

_DEFAULT_SITE = {
    "depression": "brain",
    "schizophrenia": "brain",
    "lung adenocarcinoma": "lung",
    "breast cancer": "breast",
    "type 2 diabetes": "systemic",
    "obesity": "adipose",
    "heart disease": "heart",
}


def _match_disease_catalog(text: str, catalog: list[dict[str, Any]] | None) -> str | None:
    if not catalog:
        return None
    low = text.lower()
    # Prefer longest alias hits
    hits: list[tuple[int, str]] = []
    for d in catalog:
        names = [d.get("name") or "", d.get("disease_id") or ""] + list(d.get("aliases") or [])
        for name in names:
            n = str(name).strip().lower()
            if len(n) < 4:
                continue
            if n in low:
                hits.append((len(n), d.get("name") or d.get("disease_id") or n))
    if not hits:
        return None
    hits.sort(reverse=True)
    return hits[0][1]


def parse_rules(text: str, *, disease_catalog: list[dict[str, Any]] | None = None) -> QuerySlots:
    """Parse a free-text clinical vignette into QuerySlots without any LLM."""
    raw = (text or "").strip()
    slots = QuerySlots(source_text=raw, parse_method="rules")
    if not raw:
        return slots

    # Age / sex
    m = re.search(r"\b(\d{1,3})\s*(?:year|yr|yo|y/?o)s?\s*old\b", raw, re.I)
    if not m:
        m = re.search(r"\b(\d{1,3})\s*yo\b", raw, re.I)
    if m:
        slots.age_years = float(m.group(1))

    if re.search(r"\b(female|woman|girl)\b", raw, re.I):
        slots.sex = "female"
    elif re.search(r"\b(male|man|boy)\b", raw, re.I):
        slots.sex = "male"

    # Smoking
    if re.search(r"\bcurrent smoker\b|\bsmokes\b|\bsmoking\b", raw, re.I):
        slots.smoking = "current"
    elif re.search(r"\bformer smoker\b|\bex-?smoker\b|\bquit smoking\b", raw, re.I):
        slots.smoking = "former"
    elif re.search(r"\bnever smoker\b|\bnon-?smoker\b", raw, re.I):
        slots.smoking = "never"

    # Stage
    m = re.search(r"\bstage\s*([IVX]+|\d+[A-C]?)\b", raw, re.I)
    if m:
        stage = m.group(1).upper()
        digit_map = {"1": "I", "2": "II", "3": "III", "4": "IV"}
        if stage[0].isdigit():
            stage = digit_map.get(stage[0], stage) + stage[1:]
        slots.stage = stage

    # Genes (KRAS, TP53, …)
    genes = re.findall(r"\b([A-Z][A-Z0-9]{1,7})\b", raw)
    deny = {
        "BMI", "CVD", "CHF", "CKD", "COPD", "MDD", "NSCLC", "SCLC", "LUAD", "T2D",
        "LLL", "RLL", "LUL", "RUL", "II", "III", "IV", "THE", "AND", "WITH", "FROM",
    }
    slots.genes = [g for g in genes if g not in deny and not g.isdigit()][:12]

    # Comorbidities (before primary — so "obese" doesn't become primary disease)
    comorbid = []
    for pat, name in _COMORBID_PATTERNS:
        if pat.search(raw):
            comorbid.append(name)
    slots.comorbidities = list(dict.fromkeys(comorbid))

    # Location
    for pat, loc in _LOCATION_PATTERNS:
        if pat.search(raw):
            slots.location = loc
            break

    # Primary disease
    disease = _match_disease_catalog(raw, disease_catalog)
    catalog_hit = None
    if disease and disease_catalog:
        dlow = str(disease).strip().lower()
        for d in disease_catalog:
            names = [d.get("name") or "", d.get("disease_id") or ""] + list(
                d.get("aliases") or []
            )
            if any(str(n).strip().lower() == dlow for n in names):
                catalog_hit = d
                break
    if not disease:
        for pat, name in _DISEASE_PATTERNS:
            if pat.search(raw):
                disease = name
                break
    # If obesity/heart matched only as comorbidity phrases and no other disease, keep comorbidity role
    if disease:
        # Avoid setting primary=obesity when text is "depression … obese"
        if disease.lower() in {"obesity", "heart disease"} and any(
            p.search(raw) for p, _ in _DISEASE_PATTERNS if _.lower() not in {"obesity", "heart disease"}
        ):
            # another primary pattern exists — re-resolve without those
            for pat, name in _DISEASE_PATTERNS:
                if name.lower() in {"obesity", "heart disease"}:
                    continue
                if pat.search(raw):
                    disease = name
                    break
        slots.disease = disease

    # Strip primary from comorbidities if duplicated
    if slots.disease:
        did = slots.disease.lower().replace(" ", "_")
        slots.comorbidities = [
            c
            for c in slots.comorbidities
            if c.lower() not in {did, slots.disease.lower().replace(" ", "_")}
            and c.lower().replace("_", " ") != slots.disease.lower()
        ]

    # Default location from disease if unset
    if slots.disease and not slots.location:
        if catalog_hit and catalog_hit.get("default_site"):
            slots.location = catalog_hit["default_site"]
        else:
            key = slots.disease.lower()
            slots.location = _DEFAULT_SITE.get(key) or _DEFAULT_SITE.get(
                key.replace("_", " ")
            )
            if not slots.location:
                if "depress" in key:
                    slots.location = "brain"
                elif "schizophren" in key:
                    slots.location = "brain"
                elif "lung" in key:
                    slots.location = "lung"

    # Metastatic hint
    if re.search(r"\bmetastat", raw, re.I):
        slots.metastatic = True

    return slots
