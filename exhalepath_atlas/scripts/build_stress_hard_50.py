#!/usr/bin/env python3
"""Build 50 adversarial / difficult stress cases for ExhalePath."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTS = [
    ROOT / "data" / "knowledge" / "stress_hard_50.json",
    ROOT / "src" / "exhalepath" / "data" / "knowledge" / "stress_hard_50.json",
]


def case(
    id: str,
    *,
    category: str,
    difficulty: str,
    disease: str,
    location: str | None = None,
    expect: dict,
    genes: list[str] | None = None,
    comorbidities: list[str] | None = None,
    comorbidity_weight: float | None = None,
    age_years: float | None = None,
    sex: str | None = None,
    smoking_status: str | None = None,
    stage: str | None = None,
    metastatic: bool = False,
    notes: str = "",
) -> dict:
    q: dict = {"disease": disease, "location": location}
    if genes is not None:
        q["genes"] = genes
    if comorbidities is not None:
        q["comorbidities"] = comorbidities
    if comorbidity_weight is not None:
        q["comorbidity_weight"] = comorbidity_weight
    if age_years is not None:
        q["age_years"] = age_years
    if sex is not None:
        q["sex"] = sex
    if smoking_status is not None:
        q["smoking_status"] = smoking_status
    if stage is not None:
        q["stage"] = stage
    if metastatic:
        q["metastatic"] = True
    return {
        "id": id,
        "category": category,
        "difficulty": difficulty,
        "query": q,
        "expect": expect,
        "notes": notes,
    }


def build() -> dict:
    cases: list[dict] = []

    # --- A. Crash / robustness (should never throw; finite outputs) ---
    cases += [
        case(
            "empty_disease",
            category="robustness",
            difficulty="hard",
            disease="",
            location="lung",
            expect={"no_exception": True, "near_healthy": True, "max_abs_log2fc": 0.4, "unresolved": True},
            notes="Empty disease must not crash; stay near-healthy custom::",
        ),
        case(
            "whitespace_disease",
            category="robustness",
            difficulty="hard",
            disease="   ",
            location="brain",
            expect={"no_exception": True, "near_healthy": True, "max_abs_log2fc": 0.4, "unresolved": True},
        ),
        case(
            "nonsense_disease",
            category="robustness",
            difficulty="hard",
            disease="zxqy_not_a_real_disease_999",
            location="liver",
            expect={"no_exception": True, "near_healthy": True, "max_abs_log2fc": 0.45, "unresolved": True},
        ),
        case(
            "symbol_location",
            category="robustness",
            difficulty="hard",
            disease="asthma",
            location="???",
            expect={
                "no_exception": True,
                "location_unmatched": True,
                "must_elevate": ["pentane", "ethane"],
            },
            notes="??? must not fuzzy-match a real tissue",
        ),
        case(
            "comorbidity_weight_overflow",
            category="robustness",
            difficulty="hard",
            disease="depression",
            location="brain",
            comorbidities=["obesity"],
            comorbidity_weight=3.5,
            expect={"no_exception": True, "must_elevate": ["indole", "acetone"]},
            notes="Weight >1.5 must clamp, not ValidationError",
        ),
        case(
            "comorbidity_weight_negative",
            category="robustness",
            difficulty="hard",
            disease="depression",
            location="brain",
            comorbidities=["obesity"],
            comorbidity_weight=-1.0,
            expect={"no_exception": True},
        ),
        case(
            "twenty_comorbidities",
            category="robustness",
            difficulty="hard",
            disease="depression",
            location="brain",
            comorbidities=[
                "obesity",
                "hypertension",
                "type_2_diabetes",
                "heart_disease",
                "copd",
                "asthma",
                "chronic_kidney_disease",
                "nafld",
                "thyroid",
                "rheumatoid_arthritis",
            ]
            * 2,
            expect={"no_exception": True, "finite": True},
        ),
    ]

    # --- B. Gene-as-disease traps ---
    cases += [
        case(
            "gene_BRCA1_as_disease",
            category="gene_trap",
            difficulty="hard",
            disease="BRCA1",
            location="breast",
            expect={
                "no_exception": True,
                "unresolved_or_not_breast": True,
                "near_healthy": True,
                "max_abs_log2fc": 0.5,
            },
            notes="BRCA1 is a gene, not breast cancer",
        ),
        case(
            "gene_KRAS_as_disease",
            category="gene_trap",
            difficulty="hard",
            disease="KRAS",
            location="lung",
            expect={"no_exception": True, "unresolved": True, "near_healthy": True, "max_abs_log2fc": 0.5},
        ),
        case(
            "gene_TP53_as_disease",
            category="gene_trap",
            difficulty="hard",
            disease="TP53",
            location="lung",
            expect={"no_exception": True, "unresolved": True, "near_healthy": True, "max_abs_log2fc": 0.5},
        ),
        case(
            "luad_with_driver_genes",
            category="gene_trap",
            difficulty="medium",
            disease="lung adenocarcinoma",
            location="lung",
            genes=["KRAS", "TP53", "EGFR"],
            stage="IV",
            smoking_status="current",
            metastatic=True,
            expect={"must_elevate": ["hexanal", "pentane"], "min_fold": {"hexanal": 1.15}},
        ),
    ]

    # --- C. Alias / resolve collisions ---
    cases += [
        case(
            "bare_pneumonia",
            category="alias",
            difficulty="hard",
            disease="pneumonia",
            location="lung",
            expect={"resolved_disease_id": "pneumonia_bacterial", "must_elevate": ["hexanal", "pentane"]},
        ),
        case(
            "bare_diabetes",
            category="alias",
            difficulty="medium",
            disease="diabetes",
            location="systemic",
            expect={"resolved_in": ["type_2_diabetes", "type1_diabetes"], "must_elevate": ["acetone"]},
        ),
        case(
            "heart_failure_not_cad",
            category="alias",
            difficulty="hard",
            disease="heart failure",
            location="heart",
            expect={"resolved_disease_id": "heart_failure", "must_elevate": ["acetone", "pentane"]},
        ),
        case(
            "nafld_not_cirrhosis",
            category="alias",
            difficulty="hard",
            disease="nafld",
            location="liver",
            expect={"resolved_disease_id": "nafld", "must_elevate": ["acetone"]},
        ),
        case(
            "autism_not_heart_failure",
            category="alias",
            difficulty="hard",
            disease="autism",
            location="brain",
            expect={"resolved_disease_id": "autism_spectrum_disorder"},
            notes="MONDO collision must not map autism→HF",
        ),
        case(
            "bare_cancer_unresolved",
            category="alias",
            difficulty="hard",
            disease="cancer",
            location="lung",
            expect={"unresolved": True, "near_healthy": True, "max_abs_log2fc": 0.45},
        ),
    ]

    # --- D. Location adversaries ---
    cases += [
        case(
            "head_neck_exact",
            category="location",
            difficulty="hard",
            disease="head_neck_cancer",
            location="head_neck",
            expect={
                "location_tissue_id": "head_neck",
                "must_elevate": ["ethanol", "hexanal"],
            },
        ),
        case(
            "pharynx_osa",
            category="location",
            difficulty="medium",
            disease="sleep_apnea",
            location="pharynx",
            expect={"location_tissue_id": "pharynx", "must_elevate": ["acetone", "isoprene"]},
        ),
        case(
            "joint_ra",
            category="location",
            difficulty="medium",
            disease="rheumatoid_arthritis",
            location="joint",
            expect={"location_tissue_id": "joint", "must_elevate": ["pentane", "hexanal"]},
        ),
        case(
            "left_breast_free_text",
            category="location",
            difficulty="medium",
            disease="breast cancer",
            location="left breast upper outer",
            expect={"location_contains": "breast", "must_elevate": ["hexanal"]},
        ),
        case(
            "retroperitoneal_mass",
            category="location",
            difficulty="hard",
            disease="ovarian cancer",
            location="retroperitoneal soft tissue mass",
            expect={"location_matched": True, "must_elevate": ["hexanal"]},
        ),
    ]

    # --- E. Clinical hard / known weak spots ---
    cases += [
        case(
            "sz_suppress_acetone",
            category="clinical_hard",
            difficulty="extreme",
            disease="schizophrenia",
            location="brain",
            age_years=18,
            sex="male",
            expect={
                "must_elevate": ["pentane", "ethane"],
                "must_suppress": ["acetone", "isoprene"],
            },
            notes="Magdeburg PTR-MS suppress panel",
        ),
        case(
            "sz_heart_suppress_tma",
            category="clinical_hard",
            difficulty="extreme",
            disease="schizophrenia",
            location="brain",
            comorbidities=["heart_disease"],
            age_years=18,
            sex="male",
            expect={
                "must_elevate": ["pentane"],
                "must_suppress": ["acetone", "trimethylamine"],
            },
        ),
        case(
            "mdd_obesity_ketones",
            category="clinical_hard",
            difficulty="medium",
            disease="depression",
            location="brain",
            comorbidities=["obesity"],
            age_years=24,
            sex="male",
            expect={"must_elevate": ["acetone", "ammonia", "indole"]},
        ),
        case(
            "t2d_acetone_gate",
            category="clinical_hard",
            difficulty="medium",
            disease="type 2 diabetes",
            location="systemic",
            age_years=55,
            sex="male",
            expect={"must_elevate": ["acetone", "2_butanone"], "min_fold": {"acetone": 1.8}},
        ),
        case(
            "copd_ethane_hexanal",
            category="clinical_hard",
            difficulty="medium",
            disease="copd",
            location="lung",
            age_years=68,
            smoking_status="former",
            expect={"must_elevate": ["ethane", "hexanal"], "min_fold": {"ethane": 1.5}},
        ),
        case(
            "malaria_benzene_acetone",
            category="clinical_hard",
            difficulty="medium",
            disease="malaria",
            location="systemic",
            expect={"must_elevate": ["benzene", "acetone", "pentane"], "must_suppress": ["isoprene"]},
        ),
        case(
            "cirrhosis_sulfur",
            category="clinical_hard",
            difficulty="hard",
            disease="cirrhosis",
            location="liver",
            expect={"must_elevate": ["dms", "ammonia"], "min_fold": {"dms": 1.4}},
        ),
        case(
            "ibd_h2s",
            category="clinical_hard",
            difficulty="hard",
            disease="inflammatory bowel disease",
            location="colon",
            expect={"must_elevate": ["hydrogen_sulfide", "pentane"]},
        ),
        case(
            "sibo_fermentation",
            category="clinical_hard",
            difficulty="hard",
            disease="SIBO",
            location="small_intestine",
            expect={"must_elevate": ["hydrogen_sulfide", "ethanol"]},
        ),
    ]

    # --- F. Zero-shot / near-healthy adversaries ---
    cases += [
        case(
            "zero_shot_freckling",
            category="zero_shot",
            difficulty="hard",
            disease="benign essential freckling syndrome ZX-0",
            location="skin",
            expect={"near_healthy": True, "max_abs_log2fc": 0.4, "unresolved": True},
        ),
        case(
            "zero_shot_madeup_syndrome",
            category="zero_shot",
            difficulty="hard",
            disease="Quigley-Hart metabolic vapor syndrome",
            location="systemic",
            expect={"near_healthy": True, "max_abs_log2fc": 0.45, "unresolved": True},
        ),
        case(
            "cjd_oxidative",
            category="zero_shot",
            difficulty="hard",
            disease="Creutzfeldt-Jakob disease",
            location="brain",
            age_years=62,
            sex="female",
            expect={"must_elevate": ["hexanal", "pentane"]},
        ),
        case(
            "epilepsy_brain_energy",
            category="zero_shot",
            difficulty="medium",
            disease="epilepsy",
            location="brain",
            expect={"must_elevate": ["acetone", "hexanal"]},
        ),
    ]

    # --- G. Conflicting / multi-signal pressure ---
    cases += [
        case(
            "obesity_alone_ketones",
            category="conflict",
            difficulty="medium",
            disease="obesity",
            location="systemic",
            expect={"must_elevate": ["acetone", "2_butanone"]},
        ),
        case(
            "paad_kras_tp53",
            category="conflict",
            difficulty="hard",
            disease="pancreatic adenocarcinoma",
            location="pancreas",
            genes=["KRAS", "TP53", "CDKN2A"],
            stage="III",
            expect={"must_elevate": ["2_butanone", "acetone", "hexanal"], "min_fold": {"2_butanone": 1.2}},
        ),
        case(
            "hcc_sulfur_not_only_ketone",
            category="conflict",
            difficulty="hard",
            disease="hepatocellular carcinoma",
            location="liver",
            genes=["TP53", "CTNNB1"],
            expect={"must_elevate": ["dms", "ammonia"]},
        ),
        case(
            "covid_aldehydes",
            category="conflict",
            difficulty="hard",
            disease="covid19",
            location="lung",
            expect={"must_elevate": ["nonanal", "octanal", "heptanal"]},
        ),
        case(
            "tb_oxidative",
            category="conflict",
            difficulty="medium",
            disease="tuberculosis",
            location="lung",
            expect={"must_elevate": ["pentane", "ethane", "benzene"]},
        ),
        case(
            "ards_ventilator",
            category="conflict",
            difficulty="hard",
            disease="ards",
            location="lung",
            comorbidities=["pneumonia_bacterial"],
            expect={"must_elevate": ["acetaldehyde", "pentane"], "must_suppress": ["isoprene"]},
        ),
        case(
            "thyroid_isoprene",
            category="conflict",
            difficulty="hard",
            disease="thyroid",
            location="thyroid",
            expect={"must_elevate": ["acetone", "isoprene"]},
        ),
        case(
            "sickle_oxidative",
            category="conflict",
            difficulty="hard",
            disease="sickle cell",
            location="blood",
            expect={"must_elevate": ["pentane", "hexanal"]},
        ),
        case(
            "parkinson_oxidative",
            category="conflict",
            difficulty="medium",
            disease="Parkinson disease",
            location="brain",
            genes=["SNCA", "PRKN"],
            expect={"must_elevate": ["pentane", "hexanal"]},
        ),
    ]

    # --- H. Metastatic / stage / sex edge ---
    cases += [
        case(
            "luad_brain_met",
            category="metastatic",
            difficulty="hard",
            disease="lung adenocarcinoma",
            location="brain",
            genes=["TP53", "KRAS"],
            stage="IV",
            metastatic=True,
            expect={"must_elevate": ["hexanal", "pentane"]},
            notes="Primary lung biology queried at brain met site",
        ),
        case(
            "breast_male",
            category="metastatic",
            difficulty="hard",
            disease="breast cancer",
            location="breast",
            sex="male",
            genes=["BRCA2"],
            expect={"must_elevate": ["hexanal"]},
        ),
        case(
            "prostate_stage_high",
            category="metastatic",
            difficulty="medium",
            disease="cancer_prostate",
            location="prostate",
            stage="IV",
            sex="male",
            age_years=72,
            expect={"must_elevate": ["hexanal", "pentane"]},
        ),
        case(
            "glioblastoma_warburg",
            category="metastatic",
            difficulty="medium",
            disease="glioblastoma",
            location="brain",
            genes=["EGFR", "PTEN", "TP53"],
            expect={"must_elevate": ["acetone", "hexanal"]},
        ),
        case(
            "sepsis_icu",
            category="metastatic",
            difficulty="hard",
            disease="sepsis",
            location="systemic",
            comorbidities=["pneumonia_bacterial"],
            age_years=72,
            expect={"must_elevate": ["hexanal", "ammonia"]},
        ),
        case(
            "influenza_not_bacterial",
            category="metastatic",
            difficulty="hard",
            disease="influenza",
            location="lung",
            expect={
                "resolved_disease_id": "influenza",
                "must_elevate": ["hexanal", "pentane"],
            },
        ),
    ]

    assert len(cases) == 50, len(cases)
    cats: dict[str, int] = {}
    for c in cases:
        cats[c["category"]] = cats.get(c["category"], 0) + 1
    return {
        "version": "1.0.0",
        "description": (
            "50 adversarial / difficult stress cases designed to break ExhalePath: "
            "empty inputs, gene-as-disease traps, alias collisions, bogus locations, "
            "comorbidity extremes, Magdeburg SZ suppressions, and zero-shot near-healthy."
        ),
        "n_cases": 50,
        "category_counts": cats,
        "cases": cases,
    }


def main() -> None:
    doc = build()
    text = json.dumps(doc, indent=2) + "\n"
    for p in OUTS:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        print(f"Wrote {p} ({doc['n_cases']} cases; {doc['category_counts']})")


if __name__ == "__main__":
    main()
