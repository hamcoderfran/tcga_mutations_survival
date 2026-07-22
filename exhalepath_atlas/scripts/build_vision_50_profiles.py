#!/usr/bin/env python3
"""Build the 50-disease vision evaluation patient profiles + ground-truth panels.

Evidence grades
---------------
A — quantified breath/metabolome fold changes (MW, Sci Data, or explicit lit ratios)
B — directional literature / priority panels with strong consensus
C — held-out / zero-shot literature expectations (may still overlap atlas priors)
D — atlas-prior–aligned approximate GT (mechanistic direction only; circularity risk)

Run: python scripts/build_vision_50_profiles.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT_PATHS = [
    ROOT / "data" / "knowledge" / "vision_50_profiles.json",
    ROOT / "src" / "exhalepath" / "data" / "knowledge" / "vision_50_profiles.json",
]


def _l2_to_fold(x: float) -> float:
    return float(2.0 ** float(x))


def _must_from_log2(measured: dict[str, float], *, elev_thr: float = 0.45, supp_thr: float = -0.35):
    must_e = [k for k, v in measured.items() if v >= elev_thr]
    should_e = [k for k, v in measured.items() if 0.25 <= v < elev_thr]
    must_s = [k for k, v in measured.items() if v <= supp_thr]
    should_s = [k for k, v in measured.items() if supp_thr < v <= -0.2]
    min_fold = {
        k: round(max(1.15, min(2.5, _l2_to_fold(v) * 0.85)), 3)
        for k, v in measured.items()
        if v >= 0.7
    }
    return must_e, should_e, must_s, should_s, min_fold


# Atlas pathway IDs only (see KnowledgeBase.pathways).
PATHWAY_MAP = {
    "lipid_peroxidation": "lipid_peroxidation",
    "oxidative_stress": "lipid_peroxidation",
    "type2_airway_inflammation": "neuroinflammation",
    "aldehyde_detoxification": "cytochrome_p450_detox",
    "ketone_body_metabolism": "ketone_body_metabolism",
    "fatty_acid_oxidation": "fatty_acid_oxidation",
    "apoptosis_necrosis": "apoptosis_necrosis",
    "one_carbon_metabolism": "one_carbon_folate",
    "mevalonate_cholesterol": "mevalonate_cholesterol",
    "microbial_fermentation": "gut_microbiome_fermentation",
    "gastric_microbiome": "gut_microbiome_fermentation",
    "microbial_proteolysis_putrefaction": "microbial_proteolysis_putrefaction",
    "neuroinflammation": "neuroinflammation",
    "brain_energy_metabolism": "brain_energy_metabolism",
    "neurotransmitter_metabolism": "neurotransmitter_metabolism",
    "urea_cycle": "urea_cycle",
    "methionine_transsulfuration": "methionine_transsulfuration",
    "glycolysis_warburg": "glycolysis_warburg",
    "mitochondrial_dysfunction": "brain_energy_metabolism",
    "hepatic_congestion": "urea_cycle",
    "tumor_lipid_peroxidation": "lipid_peroxidation",
    "alcohol_dehydrogenase_metabolism": "cytochrome_p450_detox",
    "sepsis_inflammation": "apoptosis_necrosis",
    "macrophage_oxidative_burst": "lipid_peroxidation",
    "neutrophil_oxidative_burst": "lipid_peroxidation",
    "chronic_bacterial_infection": "gut_microbiome_fermentation",
    "ventilator_induced_lung_injury": "lipid_peroxidation",
    "Plasmodium_isoprenoid_metabolism": "mevalonate_cholesterol",
    "mosquito_attractant_cue_balance": "mevalonate_cholesterol",
    "sympathetic_activation": "brain_energy_metabolism",
    "xenobiotic_retention": "cytochrome_p450_detox",
    "microbial_or_host_carbonyl_metabolism": "cytochrome_p450_detox",
    "microbial_mycobacterium_voc": "lipid_peroxidation",
}


def _map_pathways(raw: list[str] | set[str]) -> list[str]:
    out: list[str] = []
    for p in raw:
        mid = PATHWAY_MAP.get(p, p if p in PATHWAY_MAP.values() else None)
        if mid and mid not in out:
            out.append(mid)
    return out[:4]


def profile(
    *,
    id: str,
    name: str,
    location: str,
    grade: str,
    evidence_basis: str,
    patient: dict,
    must_elevate: list[str],
    should_elevate: list[str] | None = None,
    must_suppress: list[str] | None = None,
    should_suppress: list[str] | None = None,
    min_fold: dict[str, float] | None = None,
    top_k_should_include: list[str] | None = None,
    pathways: list[str] | None = None,
    cells: list[str] | None = None,
    sites: list[str] | None = None,
    query_name: str | None = None,
    genes: list[str] | None = None,
) -> dict:
    pe = patient.copy()
    if genes is not None:
        pe["genes"] = genes
    gt = {
        "must_elevate": must_elevate,
        "should_elevate": should_elevate or [],
        "must_suppress": must_suppress or [],
        "should_suppress": should_suppress or [],
        "min_fold": min_fold or {},
        "top_k_should_include": top_k_should_include or must_elevate[:3],
        "pathways_should_include": pathways or [],
        "cell_states_should_include": cells or [],
        "sites_should_include": sites or ([location] if location not in ("systemic",) else []),
    }
    return {
        "id": id,
        "name": name,
        "query_name": query_name or name,
        "location": location,
        "grade": grade,
        "evidence_basis": evidence_basis,
        "patient": pe,
        "ground_truth": gt,
    }


def _from_priority10() -> list[dict]:
    path = ROOT / "data" / "real_breath" / "literature_panels" / "priority10_voc_panels.json"
    panels = json.loads(path.read_text())["panels"]
    # Patient vignettes for priority diseases
    patients = {
        "asthma": {"age": 34, "sex": "female", "tempo": "uncontrolled", "comorbidities": []},
        "copd": {"age": 68, "sex": "male", "tempo": "severe", "smoking_status": "former", "comorbidities": []},
        "covid19": {"age": 52, "sex": "male", "tempo": "acute_hospitalized", "comorbidities": []},
        "pneumonia_bacterial": {"age": 71, "sex": "female", "tempo": "community_acquired", "comorbidities": []},
        "tuberculosis": {"age": 41, "sex": "male", "tempo": "pulmonary_active", "comorbidities": []},
        "cystic_fibrosis": {"age": 22, "sex": "female", "tempo": "chronic_infection", "comorbidities": []},
        "sleep_apnea": {"age": 55, "sex": "male", "tempo": "moderate_osa", "comorbidities": ["obesity"]},
        "cancer_stomach": {"age": 62, "sex": "male", "tempo": "stage_III", "comorbidities": []},
        "head_neck_cancer": {"age": 59, "sex": "male", "tempo": "stage_III", "smoking_status": "former", "comorbidities": []},
        "cancer_prostate": {"age": 67, "sex": "male", "tempo": "localized_high_grade", "comorbidities": []},
        "heart_failure": {"age": 74, "sex": "female", "tempo": "HFrEF_NYHA_III", "comorbidities": []},
        "malaria": {"age": 28, "sex": "male", "tempo": "falciparum_acute", "comorbidities": []},
        "ards": {"age": 58, "sex": "female", "tempo": "ventilated", "comorbidities": ["pneumonia_bacterial"]},
    }
    grades = {
        "asthma": "A",
        "copd": "A",
        "heart_failure": "A",
        "malaria": "A",
        "covid19": "B",
        "pneumonia_bacterial": "B",
        "tuberculosis": "B",
        "cystic_fibrosis": "B",
        "sleep_apnea": "B",
        "cancer_stomach": "B",
        "head_neck_cancer": "B",
        "cancer_prostate": "B",
        "ards": "B",
    }
    evidence = {
        "asthma": "Quantified ethane elevation (Olopade 2000) + Sci Data 2024 directional panel",
        "copd": "Quantified ethane + EBC aldehydes (hexanal/heptanal) + Sci Data 2024",
        "heart_failure": "Metabolomics Workbench ST000587 exhaled VOC HF vs control",
        "malaria": "Metabolomics Workbench ST000883 malaria breath panel",
        "covid19": "Priority-10 literature directional breath carbonyl/ketone panel",
        "pneumonia_bacterial": "Priority-10 literature directional pneumonia breath panel",
        "tuberculosis": "Priority-10 literature TB breath VOC panel",
        "cystic_fibrosis": "Priority-10 literature CF breath VOC panel",
        "sleep_apnea": "Priority-10 literature OSA breath VOC panel",
        "cancer_stomach": "Priority-10 gastric cancer breath VOC panel",
        "head_neck_cancer": "Priority-10 H&N cancer breath VOC panel",
        "cancer_prostate": "Priority-10 prostate cancer breath VOC panel",
        "ards": "Priority-10 ARDS / ventilator exhalate VOC panel",
    }
    locations = {
        "asthma": "lung",
        "copd": "lung",
        "covid19": "lung",
        "pneumonia_bacterial": "lung",
        "tuberculosis": "lung",
        "cystic_fibrosis": "lung",
        "sleep_apnea": "pharynx",
        "cancer_stomach": "stomach",
        "head_neck_cancer": "head_neck",
        "cancer_prostate": "prostate",
        "heart_failure": "heart",
        "malaria": "systemic",
        "ards": "lung",
    }
    sites = {
        "asthma": ["lung"],
        "copd": ["lung"],
        "covid19": ["lung"],
        "pneumonia_bacterial": ["lung"],
        "tuberculosis": ["lung"],
        "cystic_fibrosis": ["lung"],
        "sleep_apnea": ["pharynx", "lung"],
        "cancer_stomach": ["stomach", "gut"],
        "head_neck_cancer": ["head_neck", "oral"],
        "cancer_prostate": ["prostate"],
        "heart_failure": ["heart", "liver"],
        "malaria": ["blood", "liver"],
        "ards": ["lung"],
    }
    out = []
    for pan in panels:
        did = pan["disease_id"]
        m = pan["measured_log2fc"]
        must_e, should_e, must_s, should_s, min_fold = _must_from_log2(m)
        # Tighten A-grade quantified anchors
        if did == "asthma":
            must_e = ["ethane", "pentane", "2_pentanone", "hexanal"]
            should_e = ["benzaldehyde", "ethyl_acetate", "octane"]
            must_s = ["butane", "hexane", "isopropanol"]
            min_fold = {"ethane": 1.8, "pentane": 1.3}
        elif did == "copd":
            must_e = ["ethane", "hexanal", "heptanal", "pentane"]
            should_e = ["benzene", "nonanal", "isopropanol"]
            must_s = ["2_pentanone", "decanal", "isoprene"]
            min_fold = {"ethane": 2.0, "hexanal": 2.0}
        elif did == "heart_failure":
            must_e = ["acetone", "pentane", "isopropanol", "hexanal"]
            must_s = ["isoprene"]
            min_fold = {"acetone": 1.4, "pentane": 1.25}
        elif did == "malaria":
            must_e = ["benzene", "acetone", "pentane"]
            must_s = ["isoprene"]
            min_fold = {"acetone": 1.3, "pentane": 1.2}
        raw_pw = []
        for vs in (pan.get("pathway_links") or {}).values():
            raw_pw.extend(vs)
        pathways = _map_pathways(raw_pw)
        # Prefer atlas-native pathways when mapping empty
        atlas_pw = {
            "lipid_peroxidation",
            "glycolysis_warburg",
            "ketone_body_metabolism",
            "mevalonate_cholesterol",
            "fatty_acid_oxidation",
            "methionine_transsulfuration",
            "urea_cycle",
            "one_carbon_folate",
            "cytochrome_p450_detox",
            "apoptosis_necrosis",
            "Kras_mapk_proliferation",
            "pi3k_akt_mtor",
            "gut_microbiome_fermentation",
            "microbial_proteolysis_putrefaction",
            "neuroinflammation",
            "neurotransmitter_metabolism",
            "brain_energy_metabolism",
        }
        pathways = [p for p in pathways if p in atlas_pw]
        if not pathways:
            pathways = ["lipid_peroxidation"]
        cells = list(pan.get("cell_states") or [])[:3]
        # Only pass driver genes for solid tumors — immune/pathway gene lists on
        # infectious disease panels resolve to custom:: bundles and wipe atlas priors.
        cancer_ids = {"cancer_stomach", "head_neck_cancer", "cancer_prostate"}
        gene_list = list(pan.get("genes") or [])[:6] if did in cancer_ids else None
        out.append(
            profile(
                id=did,
                name=pan.get("disease_name") or did,
                query_name=did,
                location=locations.get(did, "lung"),
                grade=grades[did],
                evidence_basis=evidence[did],
                patient=patients[did],
                must_elevate=must_e[:6],
                should_elevate=should_e[:4],
                must_suppress=must_s[:4],
                should_suppress=should_s[:3],
                min_fold=min_fold,
                top_k_should_include=must_e[:3],
                pathways=pathways,
                cells=cells,
                sites=sites.get(did, [locations.get(did, "lung")]),
                genes=gene_list,
            )
        )
    return out


def _curated_extra() -> list[dict]:
    """Non-priority diseases with literature / clinical / prior-approximate GT."""
    rows = [
        profile(
            id="type_2_diabetes",
            name="type 2 diabetes",
            location="systemic",
            grade="A",
            evidence_basis="Breath acetone literature gate (PMID:21903721); ketone panel",
            patient={"age": 55, "sex": "male", "tempo": "poorly_controlled", "comorbidities": []},
            must_elevate=["acetone", "isopropanol", "2_butanone"],
            should_elevate=["pentane", "2_pentanone"],
            min_fold={"acetone": 1.8, "2_butanone": 1.2},
            top_k_should_include=["acetone", "2_butanone"],
            pathways=["ketone_body_metabolism", "fatty_acid_oxidation"],
            cells=["Ketogenic hepatocyte", "Lipolytic adipocyte"],
            sites=["liver", "adipose"],
        ),
        profile(
            id="type1_diabetes",
            name="type 1 diabetes",
            location="systemic",
            grade="B",
            evidence_basis="Ketotic breath acetone axis (stronger than T2D when uncontrolled)",
            patient={"age": 19, "sex": "female", "tempo": "ketotic", "comorbidities": []},
            must_elevate=["acetone", "2_butanone", "isopropanol"],
            should_elevate=["2_pentanone"],
            min_fold={"acetone": 2.0},
            pathways=["ketone_body_metabolism", "fatty_acid_oxidation"],
            cells=["Ketogenic hepatocyte"],
            sites=["liver"],
        ),
        profile(
            id="lung_adenocarcinoma",
            name="lung adenocarcinoma",
            location="lung",
            grade="A",
            evidence_basis="Lung cancer breath aldehyde/alkane panel (DOI:10.1038/bjc.2011.277)",
            patient={"age": 64, "sex": "male", "tempo": "stage_IV", "smoking_status": "current", "comorbidities": []},
            must_elevate=["hexanal", "heptanal", "nonanal", "pentane", "2_butanone"],
            must_suppress=["isoprene"],
            min_fold={"hexanal": 1.2, "pentane": 1.15},
            pathways=["lipid_peroxidation", "glycolysis_warburg", "apoptosis_necrosis"],
            cells=["Oxidative-stress", "Tumor"],
            sites=["lung"],
            genes=["TP53", "KRAS", "EGFR"],
        ),
        profile(
            id="breast_invasive_carcinoma",
            name="breast cancer",
            query_name="breast cancer",
            location="breast",
            grade="B",
            evidence_basis="Breast cancer breath DMTS/hexanal literature (DOI:10.3233/CBM-160587)",
            patient={"age": 58, "sex": "female", "tempo": "stage_IIA", "comorbidities": []},
            must_elevate=["dmts", "hexanal", "2_butanone"],
            min_fold={"dmts": 1.5, "hexanal": 1.15},
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Tumor", "Oxidative"],
            sites=["breast"],
            genes=["PIK3CA", "TP53"],
        ),
        profile(
            id="pancreatic_adenocarcinoma",
            name="pancreatic adenocarcinoma",
            location="pancreas",
            grade="B",
            evidence_basis="PAAD breath ketone/aldehyde literature",
            patient={"age": 66, "sex": "male", "tempo": "stage_III", "comorbidities": []},
            must_elevate=["2_butanone", "acetone", "hexanal", "acetaldehyde"],
            min_fold={"2_butanone": 1.3, "hexanal": 1.2},
            pathways=["ketone_body_metabolism", "lipid_peroxidation", "glycolysis_warburg"],
            cells=["Tumor", "Ketogenic"],
            sites=["pancreas", "liver"],
            genes=["KRAS", "TP53", "CDKN2A"],
        ),
        profile(
            id="colon_adenocarcinoma",
            name="colon adenocarcinoma",
            location="colon",
            grade="B",
            evidence_basis="CRC breath aldehyde + gut sulfur VOC literature",
            patient={"age": 61, "sex": "female", "tempo": "stage_II", "comorbidities": []},
            must_elevate=["hexanal", "hydrogen_sulfide", "acetone"],
            min_fold={"hexanal": 1.2},
            pathways=["lipid_peroxidation", "gut_microbiome_fermentation", "apoptosis_necrosis"],
            cells=["Tumor", "Colonocyte"],
            sites=["colon", "gut"],
            genes=["APC", "KRAS", "TP53"],
        ),
        profile(
            id="hepatocellular_carcinoma",
            name="hepatocellular carcinoma",
            location="liver",
            grade="B",
            evidence_basis="HCC sulfur/ammonia/limonene breath literature",
            patient={"age": 63, "sex": "male", "tempo": "stage_III", "comorbidities": ["chronic_liver_disease"]},
            must_elevate=["dms", "ammonia", "limonene", "acetone"],
            min_fold={"dms": 1.5, "ammonia": 1.3},
            pathways=["methionine_transsulfuration", "urea_cycle", "ketone_body_metabolism"],
            cells=["Sulfur-metabolizing hepatocyte", "Ketogenic hepatocyte"],
            sites=["liver"],
            genes=["TP53", "CTNNB1"],
        ),
        profile(
            id="ovarian_cancer",
            name="ovarian cancer",
            location="ovary",
            grade="B",
            evidence_basis="Ovarian cancer breath aldehyde panel literature",
            patient={"age": 57, "sex": "female", "tempo": "stage_III", "comorbidities": []},
            must_elevate=["hexanal", "nonanal", "decanal", "2_butanone"],
            min_fold={"hexanal": 1.2},
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Tumor", "Oxidative"],
            sites=["ovary"],
            genes=["TP53", "BRCA1"],
        ),
        profile(
            id="glioblastoma",
            name="glioblastoma",
            location="brain",
            grade="B",
            evidence_basis="GBM breath ketone/aldehyde literature",
            patient={"age": 54, "sex": "male", "tempo": "newly_diagnosed", "comorbidities": []},
            must_elevate=["acetone", "hexanal", "2_butanone", "ethanol"],
            should_elevate=["pentane"],
            pathways=["lipid_peroxidation", "glycolysis_warburg", "brain_energy_metabolism"],
            cells=["Metabolically stressed neuron", "Oxidative"],
            sites=["brain"],
            genes=["EGFR", "PTEN", "TP53"],
        ),
        profile(
            id="chronic_liver_disease",
            name="cirrhosis",
            query_name="cirrhosis",
            location="liver",
            grade="C",
            evidence_basis="Zero-shot cirrhosis sulfur/ammonia/limonene literature",
            patient={"age": 58, "sex": "male", "tempo": "decompensated", "comorbidities": []},
            must_elevate=["dms", "ammonia", "limonene"],
            min_fold={"dms": 1.5, "ammonia": 1.4},
            pathways=["methionine_transsulfuration", "urea_cycle"],
            cells=["Sulfur-metabolizing hepatocyte"],
            sites=["liver"],
        ),
        profile(
            id="chronic_kidney_disease",
            name="chronic kidney disease",
            location="kidney",
            grade="C",
            evidence_basis="CKD breath ammonia literature (zero-shot style)",
            patient={"age": 70, "sex": "female", "tempo": "stage_4", "comorbidities": []},
            must_elevate=["ammonia", "acetone", "dms"],
            min_fold={"ammonia": 1.4},
            pathways=["urea_cycle", "ketone_body_metabolism"],
            cells=["Hepatocyte", "Kidney"],
            sites=["kidney", "liver"],
        ),
        profile(
            id="inflammatory_bowel_disease",
            name="inflammatory bowel disease",
            location="gut",
            grade="C",
            evidence_basis="IBD breath H2S/pentane literature",
            patient={"age": 36, "sex": "female", "tempo": "active_flare", "comorbidities": []},
            must_elevate=["hydrogen_sulfide", "pentane"],
            should_elevate=["hexanal"],
            min_fold={"hydrogen_sulfide": 1.5},
            pathways=["gut_microbiome_fermentation", "lipid_peroxidation", "neuroinflammation"],
            cells=["Colonocyte", "Neutrophil", "Microbiome"],
            sites=["gut", "colon"],
        ),
        profile(
            id="gut_dysbiosis",
            name="gut dysbiosis",
            location="gut",
            grade="C",
            evidence_basis="Dysbiosis putrefaction VOC literature",
            patient={"age": 42, "sex": "male", "tempo": "chronic", "comorbidities": []},
            must_elevate=["indole", "phenol", "hydrogen_sulfide", "ammonia"],
            min_fold={"indole": 1.5, "hydrogen_sulfide": 1.4},
            pathways=["microbial_proteolysis_putrefaction", "gut_microbiome_fermentation", "urea_cycle"],
            cells=["Gut", "Microbiome"],
            sites=["gut"],
        ),
        profile(
            id="sibo",
            name="SIBO",
            location="small_intestine",
            grade="C",
            evidence_basis="SIBO fermentation VOC literature",
            patient={"age": 45, "sex": "female", "tempo": "symptomatic", "comorbidities": []},
            must_elevate=["hydrogen_sulfide", "ethanol", "indole"],
            min_fold={"hydrogen_sulfide": 1.8, "ethanol": 1.4},
            pathways=["gut_microbiome_fermentation", "microbial_proteolysis_putrefaction"],
            cells=["Gut", "Microbiome"],
            sites=["small_intestine", "gut"],
        ),
        profile(
            id="alzheimer_disease",
            name="Alzheimer disease",
            location="brain",
            grade="B",
            evidence_basis="AD breath oxidative VOC literature",
            patient={"age": 78, "sex": "female", "tempo": "moderate", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            min_fold={"hexanal": 1.3, "pentane": 1.25},
            pathways=["lipid_peroxidation", "neuroinflammation", "brain_energy_metabolism"],
            cells=["Activated microglia", "Metabolically stressed neuron"],
            sites=["brain"],
            genes=["APP", "PSEN1"],
        ),
        profile(
            id="parkinson_disease",
            name="Parkinson disease",
            location="brain",
            grade="B",
            evidence_basis="PD breath oxidative VOC literature",
            patient={"age": 69, "sex": "male", "tempo": "mid_stage", "comorbidities": []},
            must_elevate=["pentane", "hexanal"],
            should_elevate=["acetone"],
            min_fold={"pentane": 1.3, "hexanal": 1.3},
            pathways=["lipid_peroxidation", "neuroinflammation", "brain_energy_metabolism"],
            cells=["Metabolically stressed neuron", "Activated microglia"],
            sites=["brain"],
            genes=["SNCA", "PRKN"],
        ),
        profile(
            id="creutzfeldt_jakob",
            name="Creutzfeldt-Jakob disease",
            location="brain",
            grade="C",
            evidence_basis="CJD clinical oxidative tempo profile (no large breath gold standard)",
            patient={"age": 62, "sex": "female", "tempo": "early_clinical", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "nonanal"],
            should_elevate=["heptanal", "ethane"],
            min_fold={"hexanal": 1.25, "pentane": 1.2},
            pathways=["neuroinflammation", "lipid_peroxidation", "apoptosis_necrosis", "brain_energy_metabolism"],
            cells=["Activated microglia", "Metabolically stressed neuron"],
            sites=["brain"],
        ),
        profile(
            id="major_depressive_disorder",
            name="depression",
            query_name="depression",
            location="brain",
            grade="B",
            evidence_basis="MDD breath + ST003181 clinical comorbidity benchmarks",
            patient={"age": 24, "sex": "male", "tempo": "moderate", "comorbidities": []},
            must_elevate=["indole", "ammonia", "pentane", "hexanal"],
            should_elevate=["trimethylamine"],
            must_suppress=["isoprene"],
            pathways=["neurotransmitter_metabolism", "neuroinflammation", "lipid_peroxidation"],
            cells=["Neuron", "Microglia"],
            sites=["brain"],
        ),
        profile(
            id="schizophrenia",
            name="schizophrenia",
            location="brain",
            grade="B",
            evidence_basis="Magdeburg PTR-MS schizophrenia breath cohort",
            patient={"age": 18, "sex": "male", "tempo": "first_episode", "comorbidities": []},
            must_elevate=["pentane", "ethane", "ammonia"],
            must_suppress=["acetone", "isoprene"],
            pathways=["lipid_peroxidation", "neuroinflammation", "brain_energy_metabolism"],
            cells=["Neuron", "Microglia"],
            sites=["brain"],
        ),
        profile(
            id="obesity",
            name="obesity",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned metabolic ketone axis (approximate GT)",
            patient={"age": 48, "sex": "female", "tempo": "class_II", "comorbidities": []},
            must_elevate=["acetone", "2_butanone", "isopropanol"],
            should_elevate=["ammonia"],
            pathways=["ketone_body_metabolism", "fatty_acid_oxidation"],
            cells=["Lipolytic adipocyte", "Ketogenic hepatocyte"],
            sites=["adipose", "liver"],
        ),
        profile(
            id="nafld",
            name="NAFLD",
            location="liver",
            grade="D",
            evidence_basis="Atlas prior–aligned hepatic ketone/ammonia axis (approximate GT)",
            patient={"age": 51, "sex": "male", "tempo": "NASH_suspect", "comorbidities": ["obesity"]},
            must_elevate=["acetone", "2_butanone", "ammonia"],
            pathways=["ketone_body_metabolism", "fatty_acid_oxidation", "urea_cycle"],
            cells=["Ketogenic hepatocyte"],
            sites=["liver"],
        ),
        profile(
            id="sepsis",
            name="sepsis",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned oxidative/inflammatory VOC direction (approximate GT)",
            patient={"age": 72, "sex": "male", "tempo": "ICU", "comorbidities": ["pneumonia_bacterial"]},
            must_elevate=["hexanal", "ammonia", "acetone"],
            should_elevate=["acetonitrile"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Oxidative", "Neutrophil", "Macrophage"],
            sites=["lung", "liver"],
        ),
        profile(
            id="influenza",
            name="influenza",
            location="lung",
            grade="D",
            evidence_basis="Atlas prior–aligned viral respiratory oxidative panel (approximate GT)",
            patient={"age": 39, "sex": "female", "tempo": "acute", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Airway", "Alveolar"],
            sites=["lung"],
        ),
        profile(
            id="hiv",
            name="HIV",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned chronic infection VOC direction (approximate GT)",
            patient={"age": 44, "sex": "male", "tempo": "untreated", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Immune", "Macrophage"],
            sites=["systemic", "lung"],
        ),
        profile(
            id="helicobacter_pylori_infection",
            name="Helicobacter pylori infection",
            location="stomach",
            grade="D",
            evidence_basis="Atlas prior–aligned gastric microbiome VOC direction (approximate GT)",
            patient={"age": 47, "sex": "female", "tempo": "chronic_gastritis", "comorbidities": []},
            must_elevate=["hydrogen_sulfide", "ammonia", "acetone"],
            pathways=["gut_microbiome_fermentation", "urea_cycle"],
            cells=["Gastric", "Epithelium"],
            sites=["stomach", "gut"],
        ),
        profile(
            id="autism_spectrum_disorder",
            name="autism spectrum disorder",
            location="brain",
            grade="D",
            evidence_basis="Atlas prior–aligned neuro/gut VOC direction (approximate GT)",
            patient={"age": 12, "sex": "male", "tempo": "school_age", "comorbidities": []},
            must_elevate=["indole", "pentane", "hexanal"],
            pathways=["neurotransmitter_metabolism", "gut_microbiome_fermentation", "neuroinflammation"],
            cells=["Neuron", "Gut"],
            sites=["brain", "gut"],
        ),
        profile(
            id="multiple_sclerosis",
            name="multiple sclerosis",
            location="brain",
            grade="D",
            evidence_basis="Atlas prior–aligned neuroinflammatory oxidative panel (approximate GT)",
            patient={"age": 38, "sex": "female", "tempo": "relapsing", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            pathways=["neuroinflammation", "lipid_peroxidation", "brain_energy_metabolism"],
            cells=["Activated microglia", "Oligodendrocyte"],
            sites=["brain"],
        ),
        profile(
            id="epilepsy",
            name="epilepsy",
            location="brain",
            grade="D",
            evidence_basis="Atlas prior–aligned brain energy / oxidative panel (approximate GT)",
            patient={"age": 31, "sex": "male", "tempo": "drug_resistant", "comorbidities": []},
            must_elevate=["acetone", "hexanal", "pentane"],
            pathways=["brain_energy_metabolism", "lipid_peroxidation", "neuroinflammation"],
            cells=["Metabolically stressed neuron"],
            sites=["brain"],
        ),
        profile(
            id="rheumatoid_arthritis",
            name="rheumatoid arthritis",
            location="joint",
            grade="D",
            evidence_basis="Atlas prior–aligned inflammatory oxidative panel (approximate GT)",
            patient={"age": 56, "sex": "female", "tempo": "active", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "ethane"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis", "neuroinflammation"],
            cells=["Macrophage", "Synovial", "Oxidative"],
            sites=["joint", "systemic"],
        ),
        profile(
            id="heart_disease",
            name="heart disease",
            location="heart",
            grade="D",
            evidence_basis="Atlas prior–aligned cardiac ischemia / ketone axis (approximate GT)",
            patient={"age": 65, "sex": "male", "tempo": "CAD", "comorbidities": ["hypertension"]},
            must_elevate=["acetone", "pentane", "hexanal"],
            pathways=["ketone_body_metabolism", "lipid_peroxidation", "fatty_acid_oxidation"],
            cells=["Cardiomyocyte", "Endothelial"],
            sites=["heart"],
        ),
        profile(
            id="hypertension",
            name="hypertension",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned vascular oxidative panel (approximate GT)",
            patient={"age": 60, "sex": "female", "tempo": "uncontrolled", "comorbidities": []},
            must_elevate=["pentane", "hexanal", "acetone"],
            pathways=["lipid_peroxidation", "cytochrome_p450_detox"],
            cells=["Endothelial", "Oxidative"],
            sites=["systemic", "heart"],
        ),
        profile(
            id="endometriosis",
            name="endometriosis",
            location="pelvis",
            grade="D",
            evidence_basis="Atlas prior–aligned inflammatory oxidative panel (approximate GT)",
            patient={"age": 33, "sex": "female", "tempo": "moderate", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Oxidative", "Immune"],
            sites=["pelvis", "ovary"],
        ),
        profile(
            id="pcos",
            name="PCOS",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned metabolic ketone axis (approximate GT)",
            patient={"age": 29, "sex": "female", "tempo": "insulin_resistant", "comorbidities": ["obesity"]},
            must_elevate=["acetone", "2_butanone", "isopropanol"],
            pathways=["ketone_body_metabolism", "fatty_acid_oxidation"],
            cells=["Lipolytic adipocyte", "Ketogenic hepatocyte"],
            sites=["ovary", "adipose", "liver"],
        ),
        profile(
            id="thyroid",
            name="thyroid disease",
            query_name="thyroid",
            location="thyroid",
            grade="D",
            evidence_basis="Atlas prior–aligned thyroid metabolic VOC direction (approximate GT)",
            patient={"age": 46, "sex": "female", "tempo": "hyperthyroid", "comorbidities": []},
            must_elevate=["acetone", "isoprene", "pentane"],
            pathways=["mevalonate_cholesterol", "ketone_body_metabolism", "fatty_acid_oxidation"],
            cells=["Thyroid", "Hepatocyte"],
            sites=["thyroid", "systemic"],
        ),
        profile(
            id="sickle_cell",
            name="sickle cell disease",
            location="systemic",
            grade="D",
            evidence_basis="Atlas prior–aligned hemolytic / oxidative panel (approximate GT)",
            patient={"age": 27, "sex": "male", "tempo": "steady_state", "comorbidities": []},
            must_elevate=["pentane", "hexanal", "acetone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis"],
            cells=["Erythrocyte", "Oxidative"],
            sites=["blood", "systemic"],
        ),
        profile(
            id="traumatic_brain_injury",
            name="traumatic brain injury",
            location="brain",
            grade="D",
            evidence_basis="Atlas prior–aligned acute brain oxidative panel (approximate GT)",
            patient={"age": 35, "sex": "male", "tempo": "acute", "comorbidities": []},
            must_elevate=["hexanal", "pentane", "acetone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis", "brain_energy_metabolism"],
            cells=["Metabolically stressed neuron", "Activated microglia"],
            sites=["brain"],
        ),
        profile(
            id="esophageal_cancer",
            name="esophageal cancer",
            location="esophagus",
            grade="D",
            evidence_basis="Atlas prior–aligned GI tumor aldehyde panel (approximate GT)",
            patient={"age": 64, "sex": "male", "tempo": "stage_III", "smoking_status": "former", "comorbidities": []},
            must_elevate=["hexanal", "acetone", "2_butanone"],
            pathways=["lipid_peroxidation", "apoptosis_necrosis", "glycolysis_warburg"],
            cells=["Tumor", "Epithelium"],
            sites=["esophagus"],
            genes=["TP53"],
        ),
    ]
    return rows


def main() -> None:
    profiles = _from_priority10() + _curated_extra()
    # Ensure uniqueness and exactly 50
    by_id = {}
    for p in profiles:
        by_id[p["id"]] = p
    profiles = list(by_id.values())
    if len(profiles) < 50:
        raise SystemExit(f"Need 50 profiles, got {len(profiles)}")
    # Prefer curated order: keep first 50 after stable sort by grade then id
    grade_rank = {"A": 0, "B": 1, "C": 2, "D": 3}
    profiles = sorted(profiles, key=lambda p: (grade_rank.get(p["grade"], 9), p["id"]))
    if len(profiles) > 50:
        # Drop lowest-priority extras if over: prefer keeping all A/B/C
        keep = [p for p in profiles if p["grade"] in ("A", "B", "C")]
        d_only = [p for p in profiles if p["grade"] == "D"]
        need = 50 - len(keep)
        profiles = keep + d_only[:need]
    assert len(profiles) == 50, len(profiles)

    grades = {}
    for p in profiles:
        grades[p["grade"]] = grades.get(p["grade"], 0) + 1

    doc = {
        "version": "1.0.0",
        "description": (
            "50 patient profiles for ExhalePath vision evaluation: disease + location + "
            "demographics → ranked VOC Δppb with pathway/cell/site explainability, scored "
            "against measured or high-accuracy approximate ground truth."
        ),
        "scoring_notes": {
            "vision_fidelity_pct": "Mean composite score across all 50 profiles (includes Grade D prior-aligned cases).",
            "evidence_backed_pct": "Mean composite over Grade A+B only (measured/directional literature).",
            "held_out_style_pct": "Mean composite over Grade A+B+C (excludes prior-only Grade D).",
            "grade_D_circularity": (
                "Grade D ground truth is derived from atlas mechanistic priors; high scores "
                "there measure prior consistency, not independent clinical accuracy."
            ),
        },
        "n_profiles": 50,
        "grade_counts": grades,
        "profiles": profiles,
    }
    text = json.dumps(doc, indent=2) + "\n"
    for out in OUT_PATHS:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        print(f"Wrote {out} ({len(profiles)} profiles; grades={grades})")


if __name__ == "__main__":
    main()
