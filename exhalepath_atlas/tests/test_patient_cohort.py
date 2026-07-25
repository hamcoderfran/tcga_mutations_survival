"""Smoke tests for 1000-patient cohort builder / eval."""

from __future__ import annotations

from pathlib import Path

from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge
from exhalepath.eval.patient_cohort import run_cohort


def test_bronchitis_and_lusc_resolve():
    clear_knowledge_cache()
    kb = default_knowledge()
    assert kb.resolve_disease("bronchitis")["disease_id"] == "chronic_bronchitis"
    assert kb.resolve_disease("COPD")["disease_id"] == "copd"
    assert kb.resolve_disease("LUSC")["disease_id"] == "lung_squamous_cell_carcinoma"
    assert kb.resolve_disease("lung cancer")["disease_id"] == "lung_adenocarcinoma"


def test_cohort_smoke(tmp_path: Path):
    # Tiny synthetic cohort
    cohort = {
        "version": "test",
        "n_patients": 6,
        "focus_triad": [
            "copd",
            "chronic_bronchitis",
            "lung_adenocarcinoma",
            "lung_squamous_cell_carcinoma",
        ],
        "patients": [
            {
                "patient_id": "T1",
                "stratum": "focus",
                "site_match": "matched",
                "query": {
                    "disease": "copd",
                    "location": "lung",
                    "genes": ["SERPINA1", "HHIP"],
                    "smoking_status": "current",
                    "age_years": 68,
                    "sex": "male",
                },
                "meta": {"disease_id": "copd", "category": "pulmonary"},
            },
            {
                "patient_id": "T2",
                "stratum": "focus",
                "site_match": "matched",
                "query": {
                    "disease": "chronic_bronchitis",
                    "location": "lung",
                    "genes": ["TNF", "IL6"],
                    "smoking_status": "former",
                    "age_years": 61,
                    "sex": "female",
                },
                "meta": {"disease_id": "chronic_bronchitis", "category": "pulmonary"},
            },
            {
                "patient_id": "T3",
                "stratum": "focus",
                "site_match": "matched",
                "query": {
                    "disease": "lung_adenocarcinoma",
                    "location": "lung",
                    "genes": ["KRAS", "TP53", "EGFR"],
                    "smoking_status": "current",
                    "stage": "IV",
                    "metastatic": True,
                    "age_years": 64,
                    "sex": "male",
                },
                "meta": {"disease_id": "lung_adenocarcinoma", "category": "cancer"},
            },
            {
                "patient_id": "T4",
                "stratum": "focus",
                "site_match": "matched",
                "query": {
                    "disease": "lung_squamous_cell_carcinoma",
                    "location": "lung",
                    "genes": ["TP53", "CDKN2A"],
                    "smoking_status": "current",
                    "stage": "III",
                    "age_years": 70,
                    "sex": "male",
                },
                "meta": {
                    "disease_id": "lung_squamous_cell_carcinoma",
                    "category": "cancer",
                },
            },
            {
                "patient_id": "T5",
                "stratum": "diverse",
                "site_match": "matched",
                "query": {
                    "disease": "type_2_diabetes",
                    "location": "systemic",
                    "smoking_status": "never",
                    "age_years": 55,
                    "sex": "female",
                },
                "meta": {"disease_id": "type_2_diabetes", "category": "metabolic"},
            },
            {
                "patient_id": "T6",
                "stratum": "diverse",
                "site_match": "mismatched",
                "query": {
                    "disease": "asthma",
                    "location": "brain",
                    "genes": ["IL13"],
                    "smoking_status": "never",
                    "age_years": 28,
                    "sex": "female",
                },
                "meta": {"disease_id": "asthma", "category": "pulmonary"},
            },
        ],
    }
    path = tmp_path / "cohort.json"
    path.write_text(__import__("json").dumps(cohort))
    out = tmp_path / "out"
    report = run_cohort(out_dir=out, cohort_path=path)
    assert report["overall"]["n_ok"] == 6
    assert (out / "patients_dense.csv").exists()
    assert (out / "patient_voc_matrix.csv").exists()
    assert (out / "patient_gene_shift_matrix.csv").exists()
    assert (out / "figures" / "voc_pca.png").exists()
    assert (out / "figures" / "voc_pca_multidisease.png").exists()
    assert (out / "disease_relatability.json").exists()
    assert (out / "COHORT_1000.md").exists()
