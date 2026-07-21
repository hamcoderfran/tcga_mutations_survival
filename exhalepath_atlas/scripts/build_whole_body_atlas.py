#!/usr/bin/env python3
"""Build whole-body tissue map, 50-VOC panel, and 100-disease prior atlas."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KNOW = ROOT / "data" / "knowledge"
CENSUS = ROOT / "data" / "census"

# Additional VOCs to reach a 50-biomarker exhaled panel (literature-typical healthy ppb)
EXTRA_VOCS = [
    ("propanol", "1-Propanol", "71-23-8", "C3H8O", 12.0, 1.0, 80.0, ["microbiome", "lipid_oxidation"]),
    ("isopropanol", "2-Propanol", "67-63-0", "C3H8O", 25.0, 2.0, 150.0, ["ketone_body_metabolism"]),
    ("butane", "n-Butane", "106-97-8", "C4H10", 3.0, 0.3, 20.0, ["lipid_peroxidation"]),
    ("hexane", "n-Hexane", "110-54-3", "C6H14", 2.5, 0.2, 15.0, ["lipid_peroxidation"]),
    ("octane", "n-Octane", "111-65-9", "C8H18", 1.2, 0.1, 8.0, ["lipid_peroxidation"]),
    ("styrene", "Styrene", "100-42-5", "C8H8", 0.6, 0.05, 5.0, ["aromatic_hydrocarbon_metabolism"]),
    ("ethylbenzene", "Ethylbenzene", "100-41-4", "C8H10", 0.7, 0.05, 6.0, ["aromatic_hydrocarbon_metabolism"]),
    ("xylene", "Xylene", "1330-20-7", "C8H10", 0.9, 0.1, 8.0, ["aromatic_hydrocarbon_metabolism"]),
    ("cyclohexane", "Cyclohexane", "110-82-7", "C6H12", 1.5, 0.1, 10.0, ["lipid_peroxidation"]),
    ("methyl_acetate", "Methyl acetate", "79-20-9", "C3H6O2", 2.0, 0.2, 15.0, ["ester_metabolism"]),
    ("ethyl_acetate", "Ethyl acetate", "141-78-6", "C4H8O2", 3.5, 0.3, 25.0, ["ester_metabolism", "microbiome"]),
    ("acetonitrile", "Acetonitrile", "75-05-8", "C2H3N", 1.0, 0.1, 10.0, ["smoking", "nitrogen_metabolism"]),
    ("furan", "Furan", "110-00-9", "C4H4O", 0.8, 0.05, 6.0, ["lipid_peroxidation"]),
    ("2_pentanone", "2-Pentanone", "107-87-9", "C5H10O", 1.4, 0.1, 10.0, ["fatty_acid_oxidation"]),
    ("3_methylbutanal", "3-Methylbutanal", "590-86-3", "C5H10O", 1.8, 0.2, 12.0, ["protein_catabolism"]),
    ("benzaldehyde", "Benzaldehyde", "100-52-7", "C7H6O", 1.2, 0.1, 10.0, ["aromatic_hydrocarbon_metabolism"]),
    ("octanal", "Octanal", "124-13-0", "C8H16O", 2.2, 0.2, 12.0, ["lipid_peroxidation"]),
    ("undecane", "Undecane", "1120-21-4", "C11H24", 1.0, 0.1, 8.0, ["lipid_peroxidation"]),
    ("dodecane", "Dodecane", "112-40-3", "C12H26", 0.9, 0.1, 7.0, ["lipid_peroxidation"]),
    ("carbon_disulfide", "Carbon disulfide", "75-15-0", "CS2", 0.5, 0.05, 4.0, ["methionine_metabolism"]),
    ("methyl_mercaptan", "Methanethiol", "74-93-1", "CH4S", 0.4, 0.02, 5.0, ["methionine_metabolism", "microbiome"]),
    ("trimethylamine", "Trimethylamine", "75-50-3", "C3H9N", 2.0, 0.2, 20.0, ["microbiome", "choline_metabolism"]),
    ("dimethyl_amine", "Dimethylamine", "124-40-3", "C2H7N", 3.0, 0.3, 25.0, ["protein_catabolism"]),
    ("pyrrole", "Pyrrole", "109-97-7", "C4H5N", 0.5, 0.05, 4.0, ["microbiome"]),
    ("allyl_methyl_sulfide", "Allyl methyl sulfide", "10152-76-8", "C4H8S", 0.3, 0.02, 3.0, ["microbiome", "diet"]),
    ("propionaldehyde", "Propanal", "123-38-6", "C3H6O", 2.5, 0.2, 15.0, ["lipid_peroxidation"]),
    ("crotonaldehyde", "Crotonaldehyde", "4170-30-3", "C4H6O", 0.6, 0.05, 5.0, ["lipid_peroxidation", "smoking"]),
]

CATEGORY_TEMPLATES = {
    "cancer": {
        "pathway_bias": {
            "glycolysis_warburg": 1.4,
            "lipid_peroxidation": 1.3,
            "Kras_mapk_proliferation": 1.25,
            "apoptosis_necrosis": 1.15,
            "cytochrome_p450_detox": 1.1,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.7,
            "heptanal": 0.5,
            "nonanal": 0.5,
            "octanal": 0.45,
            "pentane": 0.45,
            "butane": 0.3,
            "hexane": 0.35,
            "2_butanone": 0.55,
            "acetaldehyde": 0.4,
            "ethanol": 0.3,
            "acetone": 0.25,
            "propionaldehyde": 0.35,
            "furan": 0.25,
            "styrene": 0.2,
            "ethylbenzene": 0.2,
        },
    },
    "neurological": {
        "pathway_bias": {
            "neuroinflammation": 1.5,
            "brain_energy_metabolism": 1.35,
            "lipid_peroxidation": 1.3,
            "neurotransmitter_metabolism": 1.25,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.65,
            "pentane": 0.6,
            "ethane": 0.4,
            "nonanal": 0.4,
            "octanal": 0.35,
            "acetone": 0.3,
            "propanol": 0.2,
            "isoprene": -0.15,
        },
    },
    "metabolic": {
        "pathway_bias": {
            "ketone_body_metabolism": 1.5,
            "fatty_acid_oxidation": 1.35,
            "urea_cycle": 1.2,
        },
        "voc_log2fc_prior": {
            "acetone": 1.0,
            "2_butanone": 0.6,
            "isopropanol": 0.5,
            "2_pentanone": 0.4,
            "ammonia": 0.35,
            "trimethylamine": 0.25,
        },
    },
    "microbiome": {
        "pathway_bias": {
            "gut_microbiome_fermentation": 1.6,
            "microbial_proteolysis_putrefaction": 1.5,
        },
        "voc_log2fc_prior": {
            "indole": 1.0,
            "phenol": 0.8,
            "hydrogen_sulfide": 0.9,
            "trimethylamine": 0.7,
            "dimethyl_amine": 0.45,
            "methyl_mercaptan": 0.55,
            "propanol": 0.4,
            "ethyl_acetate": 0.35,
            "ethanol": 0.5,
            "pyrrole": 0.3,
        },
    },
    "inflammatory": {
        "pathway_bias": {
            "lipid_peroxidation": 1.35,
            "apoptosis_necrosis": 1.2,
            "gut_microbiome_fermentation": 1.2,
        },
        "voc_log2fc_prior": {
            "pentane": 0.55,
            "hexanal": 0.5,
            "octanal": 0.35,
            "hydrogen_sulfide": 0.45,
            "3_methylbutanal": 0.3,
        },
    },
    "infectious": {
        "pathway_bias": {
            "lipid_peroxidation": 1.25,
            "apoptosis_necrosis": 1.2,
            "cytochrome_p450_detox": 1.15,
        },
        "voc_log2fc_prior": {
            "acetonitrile": 0.4,
            "hexanal": 0.45,
            "acetone": 0.3,
            "ammonia": 0.35,
            "crotonaldehyde": 0.25,
        },
    },
    "cardiovascular": {
        "pathway_bias": {
            "lipid_peroxidation": 1.3,
            "fatty_acid_oxidation": 1.2,
            "mevalonate_cholesterol": 1.15,
        },
        "voc_log2fc_prior": {
            "pentane": 0.5,
            "hexanal": 0.45,
            "isoprene": -0.2,
            "acetone": 0.25,
            "2_pentanone": 0.25,
        },
    },
    "pulmonary": {
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "cytochrome_p450_detox": 1.25,
            "glycolysis_warburg": 1.15,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.7,
            "pentane": 0.55,
            "benzene": 0.35,
            "toluene": 0.3,
            "xylene": 0.25,
            "2_butanone": 0.4,
            "furan": 0.25,
            "styrene": 0.2,
        },
    },
    "default": {
        "pathway_bias": {"lipid_peroxidation": 1.15, "glycolysis_warburg": 1.1},
        "voc_log2fc_prior": {
            "hexanal": 0.3,
            "acetone": 0.2,
            "pentane": 0.25,
            "octanal": 0.2,
        },
    },
}

# Heuristic category from top100 id / name
def infer_category(d: dict) -> str:
    i = d["id"]
    n = (d.get("name") or "").lower()
    if i.startswith("cancer_") or "cancer" in n or i in {"glioblastoma", "neuroblastoma", "mesothelioma", "sarcoma", "wilms", "retinoblastoma"}:
        return "cancer"
    if i in {"alzheimer", "parkinson", "schizophrenia", "depression", "bipolar", "epilepsy", "multiple_sclerosis", "als", "stroke", "autism", "adhd", "migraine", "huntington", "lewy_body", "frontotemporal", "ptsd", "anxiety"}:
        return "neurological"
    if i in {"type2_diabetes", "type1_diabetes", "obesity", "hyperlipidemia", "thyroid", "gout", "gout_ckd_overlap", "nafld"}:
        return "metabolic"
    if i in {"ckd", "cirrhosis", "acute_kidney_injury"}:
        return "metabolic"
    if i in {"ibd", "celiac", "gerd", "hpylori", "cdiff", "periodontitis"}:
        return "microbiome" if i in {"hpylori", "cdiff"} else "inflammatory"
    if i in {"covid19", "influenza", "tuberculosis", "hiv", "cmv", "sepsis", "uti", "pneumonia_bacterial"}:
        return "infectious"
    if i in {"heart_disease", "hypertension", "arrhythmia", "peripheral_artery", "venous_thromboembolism", "pulmonary_hypertension"}:
        return "cardiovascular"
    if i in {"copd", "asthma", "ild", "sleep_apnea"}:
        return "pulmonary"
    if i in {"lupus", "rheumatoid_arthritis", "psoriasis", "atopic_dermatitis", "sjogren", "scleroderma", "vasculitis", "myasthenia"}:
        return "inflammatory"
    return "default"


def default_site(d: dict, category: str) -> str | None:
    i = d["id"]
    mapping = {
        "cancer_lung": "lung",
        "cancer_breast": "breast",
        "cancer_colorectal": "colon",
        "cancer_prostate": "prostate gland",
        "cancer_pancreas": "pancreas",
        "cancer_liver": "liver",
        "cancer_ovary": "ovary",
        "cancer_stomach": "stomach",
        "cancer_kidney": "kidney",
        "cancer_bladder": "bladder organ",
        "cancer_skin_melanoma": "skin of body",
        "cancer_blood": "blood",
        "glioblastoma": "brain",
        "neuroblastoma": "adrenal gland",
        "alzheimer": "brain",
        "parkinson": "brain",
        "schizophrenia": "brain",
        "depression": "brain",
        "type2_diabetes": "systemic",
        "ckd": "kidney",
        "cirrhosis": "liver",
        "ibd": "intestine",
        "covid19": "lung",
        "copd": "lung",
        "asthma": "lung",
        "heart_disease": "heart",
        "hpylori": "stomach",
        "cdiff": "colon",
    }
    if i in mapping:
        return mapping[i]
    if category == "neurological":
        return "brain"
    if category == "pulmonary":
        return "lung"
    if category == "cardiovascular":
        return "heart"
    if category == "microbiome":
        return "gut"
    return "systemic"


def main():
    # --- VOC catalog ---
    voc_doc = json.loads((KNOW / "voc_catalog.json").read_text())
    existing = {v["voc_id"] for v in voc_doc["vocs"]}
    for voc_id, name, cas, formula, med, lo, hi, sources in EXTRA_VOCS:
        if voc_id in existing:
            continue
        voc_doc["vocs"].append(
            {
                "voc_id": voc_id,
                "name": name,
                "cas": cas,
                "formula": formula,
                "healthy_ppb_median": med,
                "healthy_ppb_low": lo,
                "healthy_ppb_high": hi,
                "sources": sources,
                "refs": ["PMID:21903721", "DOI:10.1016/j.cpt.2024.12.004"],
            }
        )
    voc_doc["version"] = "2.0.0"
    voc_doc["description"] = (
        "50-biomarker exhaled VOC panel with literature-typical healthy ppb baselines "
        "for whole-body disease biomarker prediction."
    )
    (KNOW / "voc_catalog.json").write_text(json.dumps(voc_doc, indent=2))
    print(f"VOCs: {len(voc_doc['vocs'])}")

    # --- Whole-body tissue map ---
    tissues = []
    if (CENSUS / "healthy_tissue_summary.csv").exists():
        ts = pd.read_csv(CENSUS / "healthy_tissue_summary.csv")
        for r in ts.itertuples():
            tissues.append(
                {
                    "tissue_id": str(r.tissue_general).lower().replace(" ", "_"),
                    "name": str(r.tissue_general),
                    "aliases": [str(r.tissue_general), str(r.tissue_general).lower()],
                    "n_census_healthy_cells": int(r.n_cells),
                    "compartment": (
                        "cns"
                        if "brain" in str(r.tissue_general).lower() or "nervous" in str(r.tissue_general).lower()
                        else "vascular"
                        if str(r.tissue_general).lower() in {"blood", "vasculature", "heart"}
                        else "epithelial"
                    ),
                }
            )
    # Ensure common clinical aliases exist even without census
    extras = [
        ("systemic", "Systemic / whole body", ["systemic", "whole body", "body"]),
        ("gut", "Gut", ["gut", "gi tract", "gastrointestinal"]),
        ("small_intestine", "Small intestine", ["small intestine", "small bowel"]),
        ("bone", "Bone", ["bone", "skeleton"]),
        ("lymph_node", "Lymph node", ["lymph node", "lymph"]),
        ("peritoneum", "Peritoneum", ["peritoneum", "peritoneal"]),
        ("spinal_cord", "Spinal cord", ["spinal cord"]),
        ("bladder", "Bladder", ["bladder"]),
    ]
    have = {t["tissue_id"] for t in tissues}
    for tid, name, aliases in extras:
        if tid not in have:
            tissues.append(
                {
                    "tissue_id": tid,
                    "name": name,
                    "aliases": aliases,
                    "n_census_healthy_cells": 0,
                    "compartment": "other",
                }
            )
    tissue_doc = {
        "version": "1.0.0",
        "description": "Whole-body anatomic locations for affected-cell placement (Census tissues + clinical aliases).",
        "tissues": tissues,
    }
    (KNOW / "whole_body_tissues.json").write_text(json.dumps(tissue_doc, indent=2))
    print(f"Tissues: {len(tissues)}")

    # --- 100-disease atlas (merge with existing curated priors) ---
    curated = json.loads((KNOW / "disease_voc_priors.json").read_text())
    curated_by_id = {d["disease_id"]: d for d in curated["diseases"]}
    # alias index from curated
    alias_to_curated = {}
    for d in curated["diseases"]:
        for a in [d["disease_id"], d["name"], *(d.get("aliases") or [])]:
            if a:
                alias_to_curated[a.strip().lower()] = d["disease_id"]

    top100 = json.loads((KNOW / "top100_us_diseases.json").read_text())["diseases"]
    merged = list(curated["diseases"])
    seen = set(curated_by_id)

    # Explicit map top100 id → curated disease_id
    top_to_curated = {
        "cancer_lung": "lung_adenocarcinoma",
        "cancer_breast": "breast_invasive_carcinoma",
        "cancer_colorectal": "colon_adenocarcinoma",
        "cancer_pancreas": "pancreatic_adenocarcinoma",
        "cancer_liver": "hepatocellular_carcinoma",
        "cancer_ovary": "ovarian_cancer",
        "glioblastoma": "glioblastoma",
        "alzheimer": "alzheimer_disease",
        "parkinson": "parkinson_disease",
        "schizophrenia": "schizophrenia",
        "depression": "major_depressive_disorder",
        "epilepsy": "epilepsy",
        "multiple_sclerosis": "multiple_sclerosis",
        "type2_diabetes": "type_2_diabetes",
        "ckd": "chronic_kidney_disease",
        "cirrhosis": "chronic_liver_disease",
        "ibd": "inflammatory_bowel_disease",
        "autism": "autism_spectrum_disorder",
        "hpylori": "helicobacter_pylori_infection",
        "cdiff": "clostridioides_difficile_infection",
    }

    for d in top100:
        if d["id"] == "normal_baseline":
            continue
        if d["id"] in top_to_curated:
            # extend aliases on curated entry
            cid = top_to_curated[d["id"]]
            entry = next(x for x in merged if x["disease_id"] == cid)
            aliases = set(entry.get("aliases") or [])
            aliases.update(d.get("aliases") or [])
            aliases.add(d["name"])
            entry["aliases"] = sorted(aliases)
            continue
        # skip if already covered by alias
        if any(a.lower() in alias_to_curated for a in [d["name"], *(d.get("aliases") or [])]):
            continue
        cat = infer_category(d)
        tmpl = CATEGORY_TEMPLATES.get(cat, CATEGORY_TEMPLATES["default"])
        did = d["id"]
        if did in seen:
            continue
        entry = {
            "disease_id": did,
            "name": d["name"],
            "aliases": list(dict.fromkeys([d["name"], *(d.get("aliases") or []), did.replace("_", " ")])),
            "mondo_id": None,
            "category": cat if cat != "default" else "unspecified",
            "default_site": default_site(d, cat),
            "pathway_bias": dict(tmpl["pathway_bias"]),
            "voc_log2fc_prior": dict(tmpl["voc_log2fc_prior"]),
            "atlas_source": "top100_template",
        }
        merged.append(entry)
        seen.add(did)

    # Pad every disease with category-template VOC priors for the 50-panel
    # (never overwrite stronger curated values).
    for entry in merged:
        cat = (entry.get("category") or "default").lower()
        if cat in {"neurological", "neurodegenerative"}:
            tmpl_cat = "neurological"
        elif cat == "cancer":
            tmpl_cat = "cancer"
        elif cat in CATEGORY_TEMPLATES:
            tmpl_cat = cat
        else:
            tmpl_cat = "default"
        tmpl = CATEGORY_TEMPLATES[tmpl_cat]["voc_log2fc_prior"]
        priors = dict(entry.get("voc_log2fc_prior") or {})
        for vid, val in tmpl.items():
            if vid not in priors:
                priors[vid] = float(val) * 0.85
        entry["voc_log2fc_prior"] = priors

    out_diseases = {
        "version": "2.0.0",
        "description": (
            "Whole-body disease atlas (~100 US high-burden conditions). "
            "Curated entries retained; others filled from category templates + Census-informed sites."
        ),
        "diseases": merged,
    }
    (KNOW / "disease_voc_priors.json").write_text(json.dumps(out_diseases, indent=2))
    print(f"Diseases: {len(merged)}")

    # --- Extend pathway VOC effects for new panel members ---
    pw = json.loads((KNOW / "pathway_voc_map.json").read_text())
    extras_effects = {
        "lipid_peroxidation": {
            "butane": 0.35,
            "hexane": 0.4,
            "octane": 0.35,
            "cyclohexane": 0.3,
            "octanal": 0.55,
            "undecane": 0.3,
            "dodecane": 0.28,
            "propionaldehyde": 0.45,
            "crotonaldehyde": 0.4,
            "furan": 0.25,
        },
        "fatty_acid_oxidation": {"2_pentanone": 0.45, "isopropanol": 0.3},
        "ketone_body_metabolism": {"isopropanol": 0.55, "2_pentanone": 0.35},
        "microbial_proteolysis_putrefaction": {
            "trimethylamine": 0.7,
            "dimethyl_amine": 0.45,
            "pyrrole": 0.35,
            "methyl_mercaptan": 0.55,
            "allyl_methyl_sulfide": 0.3,
        },
        "gut_microbiome_fermentation": {
            "propanol": 0.4,
            "ethyl_acetate": 0.35,
            "trimethylamine": 0.4,
        },
        "methionine_transsulfuration": {
            "carbon_disulfide": 0.45,
            "methyl_mercaptan": 0.6,
        },
        "cytochrome_p450_detox": {
            "styrene": 0.35,
            "ethylbenzene": 0.35,
            "xylene": 0.35,
            "benzaldehyde": 0.3,
        },
        "apoptosis_necrosis": {"3_methylbutanal": 0.35, "dimethyl_amine": 0.25},
    }
    for p in pw["pathways"]:
        pid = p["pathway_id"]
        if pid in extras_effects:
            p.setdefault("voc_effects", {}).update(extras_effects[pid])
    pw["version"] = "2.0.0"
    (KNOW / "pathway_voc_map.json").write_text(json.dumps(pw, indent=2))

    # --- Physio constants for new VOCs ---
    phys = json.loads((KNOW / "physio_constants.json").read_text())
    lam = phys.setdefault("blood_air_partition_lambda", {})
    defaults_lam = {
        "propanol": 900,
        "isopropanol": 800,
        "butane": 0.5,
        "hexane": 0.6,
        "octane": 0.8,
        "styrene": 20,
        "ethylbenzene": 18,
        "xylene": 20,
        "cyclohexane": 1.2,
        "methyl_acetate": 100,
        "ethyl_acetate": 120,
        "acetonitrile": 300,
        "furan": 5,
        "2_pentanone": 90,
        "3_methylbutanal": 40,
        "benzaldehyde": 80,
        "octanal": 70,
        "undecane": 1.5,
        "dodecane": 2.0,
        "carbon_disulfide": 8,
        "methyl_mercaptan": 6,
        "trimethylamine": 50,
        "dimethyl_amine": 60,
        "pyrrole": 30,
        "allyl_methyl_sulfide": 10,
        "propionaldehyde": 35,
        "crotonaldehyde": 40,
    }
    for k, v in defaults_lam.items():
        lam.setdefault(k, v)
    fp = phys.setdefault("hepatic_first_pass", {})
    for k in defaults_lam:
        fp.setdefault(k, 0.2)
    (KNOW / "physio_constants.json").write_text(json.dumps(phys, indent=2))
    print("Updated pathway map + physio constants")

    # --- Extra pathway chains for 50-VOC physiology coverage ---
    chains_doc = json.loads((KNOW / "voc_pathway_chains.json").read_text())
    have_c = {c["chain_id"] for c in chains_doc["chains"]}
    extra_chains = [
        {
            "chain_id": "octanal_lipid_peroxidation",
            "voc_id": "octanal",
            "name": "PUFA peroxidation → octanal",
            "primary_tissues": ["tumor", "lung", "brain"],
            "steps": [
                {
                    "step": 1,
                    "from": "membrane_pufa",
                    "to": "octanal",
                    "enzymes": ["ALOX15", "CYBB"],
                    "pathway_ids": ["lipid_peroxidation"],
                    "spontaneous": True,
                }
            ],
            "stoichiometry": 0.9,
        },
        {
            "chain_id": "trimethylamine_choline",
            "voc_id": "trimethylamine",
            "name": "Microbial choline → TMA",
            "primary_tissues": ["gut"],
            "steps": [
                {
                    "step": 1,
                    "from": "choline",
                    "to": "trimethylamine",
                    "enzymes": ["FMO3"],
                    "pathway_ids": [
                        "microbial_proteolysis_putrefaction",
                        "gut_microbiome_fermentation",
                    ],
                }
            ],
            "stoichiometry": 1.0,
        },
        {
            "chain_id": "isopropanol_ketone",
            "voc_id": "isopropanol",
            "name": "Ketone reduction → isopropanol",
            "primary_tissues": ["liver"],
            "steps": [
                {
                    "step": 1,
                    "from": "acetone",
                    "to": "isopropanol",
                    "enzymes": ["ADH1B"],
                    "pathway_ids": ["ketone_body_metabolism", "fatty_acid_oxidation"],
                }
            ],
            "stoichiometry": 0.8,
        },
        {
            "chain_id": "propanol_fermentation",
            "voc_id": "propanol",
            "name": "Microbial fermentation → 1-propanol",
            "primary_tissues": ["gut"],
            "steps": [
                {
                    "step": 1,
                    "from": "propionate",
                    "to": "propanol",
                    "enzymes": [],
                    "spontaneous": True,
                    "pathway_ids": ["gut_microbiome_fermentation"],
                }
            ],
            "stoichiometry": 0.85,
        },
        {
            "chain_id": "xylene_cyp",
            "voc_id": "xylene",
            "name": "Aromatic hydrocarbon / CYP axis → xylene signal",
            "primary_tissues": ["lung", "liver"],
            "steps": [
                {
                    "step": 1,
                    "from": "aromatic_precursors",
                    "to": "xylene",
                    "enzymes": ["CYP1A1", "CYP2E1"],
                    "pathway_ids": ["cytochrome_p450_detox"],
                }
            ],
            "stoichiometry": 0.7,
        },
        {
            "chain_id": "furan_peroxidation",
            "voc_id": "furan",
            "name": "Oxidative stress → furan",
            "primary_tissues": ["lung", "tumor"],
            "steps": [
                {
                    "step": 1,
                    "from": "lipid_radicals",
                    "to": "furan",
                    "enzymes": [],
                    "spontaneous": True,
                    "pathway_ids": ["lipid_peroxidation"],
                }
            ],
            "stoichiometry": 0.6,
        },
        {
            "chain_id": "2_pentanone_fao",
            "voc_id": "2_pentanone",
            "name": "Fatty-acid oxidation → 2-pentanone",
            "primary_tissues": ["liver", "adipose"],
            "steps": [
                {
                    "step": 1,
                    "from": "fatty_acids",
                    "to": "2_pentanone",
                    "enzymes": ["ACADM", "HADHA"],
                    "pathway_ids": ["fatty_acid_oxidation", "ketone_body_metabolism"],
                    "spontaneous": True,
                }
            ],
            "stoichiometry": 0.85,
        },
        {
            "chain_id": "methyl_mercaptan_sulfur",
            "voc_id": "methyl_mercaptan",
            "name": "Methionine / microbial sulfur → methanethiol",
            "primary_tissues": ["liver", "gut"],
            "steps": [
                {
                    "step": 1,
                    "from": "methionine",
                    "to": "methyl_mercaptan",
                    "enzymes": ["CTH", "CBS"],
                    "pathway_ids": [
                        "methionine_transsulfuration",
                        "microbial_proteolysis_putrefaction",
                    ],
                }
            ],
            "stoichiometry": 0.9,
        },
    ]
    for c in extra_chains:
        if c["chain_id"] not in have_c:
            chains_doc["chains"].append(c)
    chains_doc["version"] = "2.0.0"
    (KNOW / "voc_pathway_chains.json").write_text(json.dumps(chains_doc, indent=2))
    print(f"Pathway chains: {len(chains_doc['chains'])}")


if __name__ == "__main__":
    main()
