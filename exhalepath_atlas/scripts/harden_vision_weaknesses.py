#!/usr/bin/env python3
"""Harden weak vision-eval diseases: tissues, cell atlas, priors, alias collisions.

Updates both data/knowledge and src/exhalepath/data/knowledge mirrors.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KDIRS = [
    ROOT / "data" / "knowledge",
    ROOT / "src" / "exhalepath" / "data" / "knowledge",
]


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


def _dump(path: Path, doc: dict) -> None:
    path.write_text(json.dumps(doc, indent=2) + "\n")


NEW_TISSUES = [
    {
        "tissue_id": "head_neck",
        "name": "head and neck",
        "aliases": ["head and neck", "head_neck", "hnscc", "oropharynx", "larynx"],
        "n_census_healthy_cells": 0,
        "compartment": "epithelial",
    },
    {
        "tissue_id": "pharynx",
        "name": "pharynx",
        "aliases": ["pharynx", "throat", "upper airway"],
        "n_census_healthy_cells": 0,
        "compartment": "epithelial",
    },
    {
        "tissue_id": "joint",
        "name": "joint / synovium",
        "aliases": ["joint", "synovium", "synovial", "articular"],
        "n_census_healthy_cells": 0,
        "compartment": "connective",
    },
    {
        "tissue_id": "pelvis",
        "name": "pelvis",
        "aliases": ["pelvis", "pelvic", "pelvic cavity"],
        "n_census_healthy_cells": 0,
        "compartment": "reproductive",
    },
    {
        "tissue_id": "thyroid",
        "name": "thyroid",
        "aliases": ["thyroid", "thyroid gland"],
        "n_census_healthy_cells": 0,
        "compartment": "endocrine",
    },
    {
        "tissue_id": "oral",
        "name": "oral cavity",
        "aliases": ["oral", "oral cavity", "mouth", "buccal"],
        "n_census_healthy_cells": 0,
        "compartment": "epithelial",
    },
    {
        "tissue_id": "prostate",
        "name": "prostate",
        "aliases": ["prostate"],
        "n_census_healthy_cells": 0,
        "compartment": "epithelial",
    },
]


NEW_CELL_STATES = [
    {
        "state_id": "airway_epithelium_stressed",
        "name": "Stressed airway epithelium",
        "tissue": "lung",
        "markers": ["SCGB1A1", "FOXJ1", "MUC5AC", "IL13", "TSLP"],
        "baseline_density": 0.08,
        "baseline_activity": 0.3,
        "produces_chains": ["hexanal_lipid_peroxidation", "pentane_peroxidation"],
        "disease_density_mult": {
            "asthma": 2.2,
            "copd": 2.4,
            "covid19": 2.0,
            "pneumonia_bacterial": 2.1,
            "tuberculosis": 1.8,
            "cystic_fibrosis": 2.3,
            "sleep_apnea": 2.0,
            "ards": 2.5,
            "influenza": 1.9,
            "hiv": 1.3,
            "head_neck_cancer": 1.4,
        },
        "disease_activity_mult": {
            "asthma": 2.0,
            "copd": 2.2,
            "covid19": 2.1,
            "pneumonia_bacterial": 2.0,
            "tuberculosis": 1.8,
            "cystic_fibrosis": 2.1,
            "sleep_apnea": 1.9,
            "ards": 2.4,
            "influenza": 1.9,
            "hiv": 1.4,
            "head_neck_cancer": 1.5,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "macrophage_activated",
        "name": "Activated macrophage",
        "tissue": "multi",
        "markers": ["CD68", "CD163", "TNF", "IL1B", "NOS2", "CYBB"],
        "baseline_density": 0.06,
        "baseline_activity": 0.25,
        "produces_chains": ["hexanal_lipid_peroxidation", "pentane_peroxidation"],
        "disease_density_mult": {
            "rheumatoid_arthritis": 2.4,
            "tuberculosis": 2.5,
            "pneumonia_bacterial": 2.2,
            "hiv": 2.0,
            "influenza": 1.8,
            "malaria": 1.6,
            "inflammatory_bowel_disease": 1.7,
            "endometriosis": 1.6,
            "sepsis": 2.0,
            "cancer_prostate": 1.4,
            "head_neck_cancer": 1.5,
            "cancer_stomach": 1.4,
        },
        "disease_activity_mult": {
            "rheumatoid_arthritis": 2.2,
            "tuberculosis": 2.3,
            "pneumonia_bacterial": 2.0,
            "hiv": 1.9,
            "influenza": 1.7,
            "malaria": 1.6,
            "inflammatory_bowel_disease": 1.7,
            "endometriosis": 1.6,
            "sepsis": 2.0,
            "cancer_prostate": 1.4,
            "head_neck_cancer": 1.5,
            "cancer_stomach": 1.4,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "erythrocyte_oxidative",
        "name": "Oxidatively stressed erythrocyte",
        "tissue": "blood",
        "markers": ["HBA1", "HBB", "SOD1", "PRDX2", "NFE2L2"],
        "baseline_density": 0.05,
        "baseline_activity": 0.2,
        "produces_chains": ["pentane_peroxidation", "hexanal_lipid_peroxidation"],
        "disease_density_mult": {
            "malaria": 3.0,
            "sickle_cell": 2.8,
            "sepsis": 1.5,
            "hiv": 1.3,
        },
        "disease_activity_mult": {
            "malaria": 2.5,
            "sickle_cell": 2.4,
            "sepsis": 1.5,
            "hiv": 1.3,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "cardiomyocyte_stressed",
        "name": "Stressed cardiomyocyte",
        "tissue": "heart",
        "markers": ["MYH7", "TNNT2", "NPPA", "NPPB", "PPARA"],
        "baseline_density": 0.05,
        "baseline_activity": 0.3,
        "produces_chains": ["pentane_peroxidation", "2_butanone_fao", "acetone_ketogenesis"],
        "disease_density_mult": {
            "heart_failure": 2.4,
            "heart_disease": 2.2,
            "hypertension": 1.6,
            "sleep_apnea": 1.4,
        },
        "disease_activity_mult": {
            "heart_failure": 2.2,
            "heart_disease": 2.0,
            "hypertension": 1.5,
            "sleep_apnea": 1.4,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "endothelial_oxidative",
        "name": "Oxidative endothelial cell",
        "tissue": "vasculature",
        "markers": ["PECAM1", "VWF", "NOS3", "SELE", "ICAM1"],
        "baseline_density": 0.04,
        "baseline_activity": 0.25,
        "produces_chains": ["pentane_peroxidation", "hexanal_lipid_peroxidation"],
        "disease_density_mult": {
            "hypertension": 2.2,
            "heart_disease": 2.0,
            "sickle_cell": 1.8,
            "diabetes": 1.5,
            "type_2_diabetes": 1.4,
            "sleep_apnea": 1.5,
        },
        "disease_activity_mult": {
            "hypertension": 2.0,
            "heart_disease": 1.9,
            "sickle_cell": 1.7,
            "type_2_diabetes": 1.4,
            "sleep_apnea": 1.5,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "synovial_inflammatory",
        "name": "Inflammatory synovial cell",
        "tissue": "joint",
        "markers": ["IL6", "TNF", "MMP3", "CXCL8", "PTGS2"],
        "baseline_density": 0.03,
        "baseline_activity": 0.25,
        "produces_chains": ["hexanal_lipid_peroxidation", "pentane_peroxidation"],
        "disease_density_mult": {
            "rheumatoid_arthritis": 3.0,
            "osteoarthritis": 1.8,
            "lupus": 1.5,
        },
        "disease_activity_mult": {
            "rheumatoid_arthritis": 2.5,
            "osteoarthritis": 1.7,
            "lupus": 1.5,
        },
        "census_calibrated": False,
    },
    {
        "state_id": "thyroid_follicular_hypermetabolic",
        "name": "Hypermetabolic thyroid follicular cell",
        "tissue": "thyroid",
        "markers": ["TG", "TPO", "TSHR", "DIO2", "SLC5A5"],
        "baseline_density": 0.04,
        "baseline_activity": 0.3,
        "produces_chains": ["isoprene_mevalonate", "acetone_ketogenesis"],
        "disease_density_mult": {"thyroid": 2.5, "thyroid_cancer": 2.0},
        "disease_activity_mult": {"thyroid": 2.2, "thyroid_cancer": 1.8},
        "census_calibrated": False,
    },
    {
        "state_id": "colonocyte_inflamed",
        "name": "Inflamed colonocyte",
        "tissue": "colon",
        "markers": ["MUC2", "LGR5", "TNF", "IL17A", "NOD2"],
        "baseline_density": 0.05,
        "baseline_activity": 0.3,
        "produces_chains": ["hexanal_lipid_peroxidation", "pentane_peroxidation"],
        "disease_density_mult": {
            "inflammatory_bowel_disease": 2.6,
            "colon_adenocarcinoma": 1.8,
            "gut_dysbiosis": 1.5,
            "sibo": 1.3,
        },
        "disease_activity_mult": {
            "inflammatory_bowel_disease": 2.3,
            "colon_adenocarcinoma": 1.7,
            "gut_dysbiosis": 1.4,
            "sibo": 1.3,
        },
        "census_calibrated": False,
    },
]


# Specialize cloned/weak disease priors (VOC + pathway + site).
PRIOR_PATCHES = {
    "rheumatoid_arthritis": {
        "default_site": "joint",
        "pathway_bias": {
            "lipid_peroxidation": 1.55,
            "apoptosis_necrosis": 1.35,
            "neuroinflammation": 1.3,
        },
        "voc_log2fc_prior": {
            "pentane": 0.7,
            "hexanal": 0.65,
            "ethane": 0.55,
            "octanal": 0.35,
            "hydrogen_sulfide": 0.25,
        },
        "atlas_source": "vision_hardened",
    },
    "thyroid": {
        "default_site": "thyroid",
        "pathway_bias": {
            "mevalonate_cholesterol": 1.45,
            "ketone_body_metabolism": 1.35,
            "fatty_acid_oxidation": 1.25,
        },
        "voc_log2fc_prior": {
            "acetone": 0.7,
            "isoprene": 0.55,
            "2_butanone": 0.4,
            "isopropanol": 0.35,
            "pentane": 0.35,
            "ammonia": 0.25,
        },
        "atlas_source": "vision_hardened",
    },
    "hiv": {
        "default_site": "systemic",
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "apoptosis_necrosis": 1.35,
            "cytochrome_p450_detox": 1.15,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.55,
            "pentane": 0.5,
            "acetone": 0.35,
            "ammonia": 0.35,
            "acetonitrile": 0.3,
        },
        "atlas_source": "vision_hardened",
    },
    "influenza": {
        "default_site": "lung",
        "aliases": [
            "Influenza",
            "influenza",
            "flu",
            "viral influenza",
            "Influenza / viral pneumonia",
            "viral pneumonia",
        ],
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "apoptosis_necrosis": 1.3,
            "cytochrome_p450_detox": 1.15,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.55,
            "pentane": 0.5,
            "acetone": 0.35,
            "ammonia": 0.3,
            "acetonitrile": 0.3,
        },
        "atlas_source": "vision_hardened",
    },
    "endometriosis": {
        "default_site": "pelvis",
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "apoptosis_necrosis": 1.3,
            "glycolysis_warburg": 1.15,
        },
        "voc_log2fc_prior": {
            "hexanal": 0.5,
            "pentane": 0.45,
            "acetone": 0.35,
            "octanal": 0.25,
        },
        "atlas_source": "vision_hardened",
    },
    "sickle_cell": {
        "default_site": "blood",
        "aliases": [
            "Sickle cell disease",
            "sickle cell disease",
            "sickle cell",
            "SCD",
        ],
        "pathway_bias": {
            "lipid_peroxidation": 1.45,
            "apoptosis_necrosis": 1.35,
            "glycolysis_warburg": 1.15,
        },
        "voc_log2fc_prior": {
            "pentane": 0.6,
            "hexanal": 0.55,
            "acetone": 0.3,
            "ethane": 0.35,
        },
        "atlas_source": "vision_hardened",
    },
    "hypertension": {
        "default_site": "heart",
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "cytochrome_p450_detox": 1.25,
            "fatty_acid_oxidation": 1.15,
        },
        "voc_log2fc_prior": {
            "pentane": 0.55,
            "hexanal": 0.5,
            "acetone": 0.3,
            "isoprene": -0.2,
        },
        "atlas_source": "vision_hardened",
    },
    "heart_disease": {
        "default_site": "heart",
        "aliases": [
            "Heart disease",
            "heart disease",
            "coronary artery disease",
            "CAD",
            "atherosclerosis",
            "myocardial infarction",
            "cardiomyopathy",
            "dilated cardiomyopathy",
        ],
        "pathway_bias": {
            "lipid_peroxidation": 1.4,
            "fatty_acid_oxidation": 1.3,
            "ketone_body_metabolism": 1.25,
            "mevalonate_cholesterol": 1.15,
        },
        "voc_log2fc_prior": {
            "pentane": 0.55,
            "hexanal": 0.5,
            "acetone": 0.4,
            "trimethylamine": 0.45,
            "dimethyl_amine": 0.4,
            "isoprene": -0.2,
        },
        "atlas_source": "vision_hardened",
    },
    "sleep_apnea": {
        "default_site": "pharynx",
        "pathway_bias": {
            "lipid_peroxidation": 1.45,
            "ketone_body_metabolism": 1.3,
            "mevalonate_cholesterol": 1.25,
            "brain_energy_metabolism": 1.2,
            "cytochrome_p450_detox": 1.15,
        },
    },
    "head_neck_cancer": {
        "default_site": "head_neck",
        "pathway_bias": {
            "glycolysis_warburg": 1.55,
            "lipid_peroxidation": 1.4,
            "apoptosis_necrosis": 1.35,
            "cytochrome_p450_detox": 1.2,
            "ketone_body_metabolism": 1.15,
        },
    },
    "cancer_prostate": {
        "pathway_bias": {
            "glycolysis_warburg": 1.4,
            "lipid_peroxidation": 1.35,
            "one_carbon_folate": 1.25,
            "cytochrome_p450_detox": 1.2,
            "apoptosis_necrosis": 1.2,
        },
    },
    "cancer_stomach": {
        "pathway_bias": {
            "glycolysis_warburg": 1.55,
            "lipid_peroxidation": 1.35,
            "ketone_body_metabolism": 1.25,
            "gut_microbiome_fermentation": 1.25,
            "mevalonate_cholesterol": 1.2,
            "apoptosis_necrosis": 1.3,
        },
    },
    "inflammatory_bowel_disease": {
        "default_site": "colon",
        "pathway_bias": {
            "gut_microbiome_fermentation": 1.5,
            "lipid_peroxidation": 1.4,
            "neuroinflammation": 1.25,
            "apoptosis_necrosis": 1.2,
        },
    },
    "anemia": {
        "aliases": [
            "Anemia",
            "anemia",
            "iron deficiency anemia",
        ],
    },
}


def patch_tissues(doc: dict) -> dict:
    existing = {t["tissue_id"] for t in doc["tissues"]}
    for t in NEW_TISSUES:
        if t["tissue_id"] in existing:
            # merge aliases
            for old in doc["tissues"]:
                if old["tissue_id"] == t["tissue_id"]:
                    aliases = list(dict.fromkeys(list(old.get("aliases") or []) + list(t["aliases"])))
                    old["aliases"] = aliases
                    break
        else:
            doc["tissues"].append(copy.deepcopy(t))
    # enrich head / lung aliases carefully (do NOT add head_neck as head alias)
    for t in doc["tissues"]:
        if t["tissue_id"] == "prostate_gland":
            aliases = list(t.get("aliases") or [])
            for a in ["prostate", "prostate gland"]:
                if a not in aliases:
                    aliases.append(a)
            t["aliases"] = aliases
        if t["tissue_id"] == "lung":
            aliases = list(t.get("aliases") or [])
            for a in ["lung", "pulmonary", "airway"]:
                if a not in aliases:
                    aliases.append(a)
            t["aliases"] = aliases
        if t["tissue_id"] == "blood":
            aliases = list(t.get("aliases") or [])
            for a in ["blood", "hematologic", "erythrocyte"]:
                if a not in aliases:
                    aliases.append(a)
            t["aliases"] = aliases
    return doc


def patch_cells(doc: dict) -> dict:
    existing = {c["state_id"] for c in doc["cell_states"]}
    for st in NEW_CELL_STATES:
        if st["state_id"] in existing:
            # replace with hardened definition
            doc["cell_states"] = [c for c in doc["cell_states"] if c["state_id"] != st["state_id"]]
        doc["cell_states"].append(copy.deepcopy(st))
    # boost existing oxidative / tumor disease mults for weak diseases
    boosts = {
        "oxidative_stress_cell": {
            "rheumatoid_arthritis": (1.8, 2.0),
            "thyroid": (1.3, 1.4),
            "hiv": (1.6, 1.8),
            "influenza": (1.7, 1.9),
            "endometriosis": (1.6, 1.7),
            "sickle_cell": (1.8, 2.0),
            "hypertension": (1.5, 1.6),
            "heart_disease": (1.6, 1.8),
            "inflammatory_bowel_disease": (1.8, 2.0),
        },
        "tumor_epithelial_warburg": {
            "head_neck_cancer": (2.6, 2.3),
            "cancer_prostate": (2.4, 2.0),
            "cancer_stomach": (2.5, 2.2),
            "esophageal_cancer": (2.3, 2.0),
        },
        "hepatocyte_ketogenic": {
            "heart_disease": (1.3, 1.5),
            "thyroid": (1.2, 1.4),
            "sleep_apnea": (1.25, 1.4),
        },
        "gastric_urease_niche": {
            "cancer_stomach": (1.8, 1.6),
            "helicobacter_pylori_infection": (2.5, 2.2),
        },
        "gut_fermentative_microbe": {
            "cancer_stomach": (1.5, 1.4),
            "inflammatory_bowel_disease": (2.0, 1.8),
            "head_neck_cancer": (1.3, 1.2),
        },
    }
    for st in doc["cell_states"]:
        bid = st["state_id"]
        if bid not in boosts:
            continue
        dens = dict(st.get("disease_density_mult") or {})
        act = dict(st.get("disease_activity_mult") or {})
        for did, (d, a) in boosts[bid].items():
            dens[did] = max(float(dens.get(did, 1.0)), d)
            act[did] = max(float(act.get(did, 1.0)), a)
        st["disease_density_mult"] = dens
        st["disease_activity_mult"] = act
    return doc


def patch_priors(doc: dict) -> dict:
    by_id = {d["disease_id"]: d for d in doc["diseases"]}
    for did, patch in PRIOR_PATCHES.items():
        if did not in by_id:
            continue
        d = by_id[did]
        for k, v in patch.items():
            d[k] = copy.deepcopy(v)
    # ensure pneumonia_bacterial keeps pneumonia alias; influenza no longer steals bare 'pneumonia'
    return doc


def main() -> None:
    for kdir in KDIRS:
        tissues_path = kdir / "whole_body_tissues.json"
        cells_path = kdir / "cell_state_atlas.json"
        priors_path = kdir / "disease_voc_priors.json"
        tissues = patch_tissues(_load(tissues_path))
        cells = patch_cells(_load(cells_path))
        priors = patch_priors(_load(priors_path))
        _dump(tissues_path, tissues)
        _dump(cells_path, cells)
        _dump(priors_path, priors)
        print(f"Patched {kdir}")
        print(f"  tissues={len(tissues['tissues'])} cells={len(cells['cell_states'])} diseases={len(priors['diseases'])}")


if __name__ == "__main__":
    main()
