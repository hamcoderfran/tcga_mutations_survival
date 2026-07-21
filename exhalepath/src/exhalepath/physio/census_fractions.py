from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd

from ..config import DATA_DIR

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
    "chronic_kidney_disease": "ckd",
    "chronic_liver_disease": "cirrhosis",
    "inflammatory_bowel_disease": "ibd",
    "autism_spectrum_disorder": "autism",
    "helicobacter_pylori_infection": "hpylori",
    "clostridioides_difficile_infection": "cdiff",
}


@lru_cache(maxsize=1)
def _disease_fractions_table(census_dir: str) -> pd.DataFrame | None:
    path = Path(census_dir) / "disease_exhalepath_state_fractions.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def census_cell_state_fractions_for_disease(
    disease_id: str,
    *,
    census_dir: Path | None = None,
) -> dict[str, float]:
    """Return ExhalePath state fractions from Census harvest for a disease, if available."""
    us_id = EXHALE_TO_US.get(disease_id)
    if not us_id:
        return {}
    table = _disease_fractions_table(str(census_dir or CENSUS_DIR))
    if table is None or table.empty:
        return {}
    sub = table[table["us_disease_id"] == us_id]
    if sub.empty:
        return {}
    # If multiple rows (shouldn't), weight by n_cells
    g = sub.groupby("exhalepath_state")["n_cells"].sum()
    total = float(g.sum()) or 1.0
    return {str(k): float(v) / total for k, v in g.items()}
