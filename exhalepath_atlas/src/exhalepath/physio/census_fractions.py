from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from ..config import DATA_DIR, KNOWLEDGE_DIR

CENSUS_DIR = DATA_DIR / "census"

# exhalepath disease_id → us_disease_id used in harvest tables
EXHALE_TO_US = {
    "lung_adenocarcinoma": "cancer_lung",
    "breast_invasive_carcinoma": "cancer_breast",
    "colon_adenocarcinoma": "cancer_colorectal",
    "pancreatic_adenocarcinoma": "cancer_pancreas",
    "hepatocellular_carcinoma": "cancer_liver",
    "ovarian_cancer": "cancer_ovary",
    "glioblastoma": "glioblastoma",
    "alzheimer_disease": "alzheimer",
    "parkinson_disease": "parkinson",
    "schizophrenia": "schizophrenia",
    "major_depressive_disorder": "depression",
    "epilepsy": "epilepsy",
    "multiple_sclerosis": "multiple_sclerosis",
    "type_2_diabetes": "type2_diabetes",
    "type_1_diabetes": "type1_diabetes",
    "chronic_kidney_disease": "ckd",
    "chronic_liver_disease": "cirrhosis",
    "inflammatory_bowel_disease": "ibd",
    "autism_spectrum_disorder": "autism",
    "helicobacter_pylori_infection": "hpylori",
    "clostridioides_difficile_infection": "cdiff",
    "copd": "copd",
    "asthma": "copd",  # closest pulmonary census proxy when asthma absent
    "covid19": "covid19",
    "covid_19": "covid19",
    "heart_disease": "heart_disease",
    "coronary_artery_disease": "heart_disease",
    "rheumatoid_arthritis": "rheumatoid_arthritis",
    "systemic_lupus_erythematosus": "lupus",
    "lupus": "lupus",
    "influenza": "influenza",
    "hiv": "hiv",
    "cystic_fibrosis": "cystic_fibrosis",
    "als": "als",
    "amyotrophic_lateral_sclerosis": "als",
    "bipolar": "bipolar",
    "bipolar_disorder": "bipolar",
    "lewy_body": "lewy_body",
    "lewy_body_dementia": "lewy_body",
    "frontotemporal": "frontotemporal",
    "frontotemporal_dementia": "frontotemporal",
    "neuroblastoma": "neuroblastoma",
    "wilms": "wilms",
    "retinoblastoma": "retinoblastoma",
    "periodontitis": "periodontitis",
    "ild": "ild",
    "interstitial_lung_disease": "ild",
    "acute_kidney_injury": "acute_kidney_injury",
    "pneumonia_bacterial": "pneumonia_bacterial",
    "sjogren": "sjogren",
    "macular_degeneration": "macular_degeneration",
    "glaucoma": "glaucoma",
    "down_syndrome": "down_syndrome",
    "cervical_cancer": "cervical_cancer",
    "esophageal_cancer": "esophageal_cancer",
    "head_neck_cancer": "head_neck_cancer",
    "cancer_lung": "cancer_lung",
    "cancer_breast": "cancer_breast",
    "cancer_colorectal": "cancer_colorectal",
    "cancer_pancreas": "cancer_pancreas",
    "cancer_liver": "cancer_liver",
    "cancer_ovary": "cancer_ovary",
    "cancer_prostate": "cancer_prostate",
    "cancer_kidney": "cancer_kidney",
    "cancer_stomach": "cancer_stomach",
    "cancer_skin_melanoma": "cancer_skin_melanoma",
    "cancer_blood": "cancer_blood",
    "type2_diabetes": "type2_diabetes",
    "type1_diabetes": "type1_diabetes",
    "ckd": "ckd",
    "ibd": "ibd",
    "alzheimer": "alzheimer",
    "parkinson": "parkinson",
}


@lru_cache(maxsize=1)
def _disease_fractions_table(census_dir: str) -> pd.DataFrame | None:
    path = Path(census_dir) / "disease_exhalepath_state_fractions.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


@lru_cache(maxsize=1)
def _us_id_set(census_dir: str) -> set[str]:
    table = _disease_fractions_table(census_dir)
    if table is None or table.empty:
        return set()
    return set(table["us_disease_id"].astype(str))


@lru_cache(maxsize=1)
def _alias_to_us() -> dict[str, str]:
    """Build fuzzy alias map from top100 list + disease atlas."""
    import json
    import re

    def norm(s: str) -> str:
        s = s.strip().lower()
        s = re.sub(r"[^a-z0-9]+", " ", s)
        return re.sub(r"\s+", " ", s).strip()

    out: dict[str, str] = {}
    top_path = KNOWLEDGE_DIR / "top100_us_diseases.json"
    if top_path.exists():
        for d in json.loads(top_path.read_text()).get("diseases") or []:
            uid = d.get("id")
            if not uid or uid == "normal_baseline":
                continue
            keys = [uid, d.get("name", "")] + list(d.get("aliases") or [])
            for k in keys:
                if k:
                    out[norm(str(k))] = uid
    # Explicit exhalepath map
    for eid, uid in EXHALE_TO_US.items():
        out[norm(eid)] = uid
        out[norm(eid.replace("_", " "))] = uid
    return out


def resolve_us_disease_id(disease_id: str, *, census_dir: Path | None = None) -> str | None:
    """Map an ExhalePath / free-text disease id to a Census harvest us_disease_id."""
    if not disease_id:
        return None
    cdir = str(census_dir or CENSUS_DIR)
    available = _us_id_set(cdir)
    if disease_id in EXHALE_TO_US:
        uid = EXHALE_TO_US[disease_id]
        return uid if uid in available else uid
    if disease_id in available:
        return disease_id
    # alias fuzzy
    import re

    key = re.sub(r"[^a-z0-9]+", " ", disease_id.strip().lower())
    key = re.sub(r"\s+", " ", key).strip()
    aliases = _alias_to_us()
    if key in aliases:
        return aliases[key]
    # containment
    for alias, uid in aliases.items():
        if key in alias or alias in key:
            if uid in available or True:
                return uid
    return None


def census_cell_state_fractions_for_disease(
    disease_id: str,
    *,
    census_dir: Path | None = None,
) -> dict[str, float]:
    """Return ExhalePath state fractions from Census harvest for a disease, if available."""
    us_id = resolve_us_disease_id(disease_id, census_dir=census_dir)
    if not us_id:
        return {}
    table = _disease_fractions_table(str(census_dir or CENSUS_DIR))
    if table is None or table.empty:
        return {}
    sub = table[table["us_disease_id"] == us_id]
    if sub.empty:
        return {}
    g = sub.groupby("exhalepath_state")["n_cells"].sum()
    total = float(g.sum()) or 1.0
    return {str(k): float(v) / total for k, v in g.items()}
