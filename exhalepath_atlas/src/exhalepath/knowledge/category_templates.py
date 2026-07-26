"""Category pathway templates + token cues for zero-shot custom diseases.

Pathway bias only — never VOC log2fc priors (those would leak breath labels).
"""

from __future__ import annotations

from typing import Any

# Pathway-only templates (mirrors scripts/build_whole_body_atlas.CATEGORY_TEMPLATES)
CATEGORY_PATHWAY_TEMPLATES: dict[str, dict[str, float]] = {
    "cancer": {
        "glycolysis_warburg": 1.4,
        "lipid_peroxidation": 1.3,
        "Kras_mapk_proliferation": 1.25,
        "apoptosis_necrosis": 1.15,
        "cytochrome_p450_detox": 1.1,
    },
    "neurological": {
        "neuroinflammation": 1.5,
        "brain_energy_metabolism": 1.35,
        "lipid_peroxidation": 1.3,
        "neurotransmitter_metabolism": 1.25,
    },
    "metabolic": {
        "ketone_body_metabolism": 1.5,
        "fatty_acid_oxidation": 1.35,
        "urea_cycle": 1.2,
    },
    "microbiome": {
        "gut_microbiome_fermentation": 1.6,
        "microbial_proteolysis_putrefaction": 1.5,
    },
    "inflammatory": {
        "lipid_peroxidation": 1.35,
        "apoptosis_necrosis": 1.2,
        "gut_microbiome_fermentation": 1.2,
    },
    "infectious": {
        "lipid_peroxidation": 1.25,
        "apoptosis_necrosis": 1.2,
        "cytochrome_p450_detox": 1.15,
    },
    "cardiovascular": {
        "lipid_peroxidation": 1.3,
        "fatty_acid_oxidation": 1.2,
        "mevalonate_cholesterol": 1.15,
    },
    "pulmonary": {
        "lipid_peroxidation": 1.4,
        "cytochrome_p450_detox": 1.25,
        "glycolysis_warburg": 1.15,
    },
    "default": {
        "lipid_peroxidation": 1.15,
        "glycolysis_warburg": 1.1,
    },
}

CATEGORY_DEFAULT_SITE: dict[str, str | None] = {
    "cancer": None,
    "neurological": "brain",
    "metabolic": "systemic",
    "microbiome": "intestine",
    "inflammatory": "systemic",
    "infectious": None,
    "cardiovascular": "heart",
    "pulmonary": "lung",
    "default": "systemic",
}

CATEGORY_DRIVER_GENES: dict[str, list[str]] = {
    "cancer": ["KRAS", "TP53", "MYC", "HIF1A", "LDHA", "HK2", "EGFR", "PIK3CA"],
    "neurological": ["APOE", "SNCA", "MAPT", "PSEN1", "SOD1", "GRIN1", "IDO1"],
    "neurodegenerative": ["APOE", "SNCA", "MAPT", "PSEN1", "SOD1"],
    "metabolic": ["HMGCS2", "CPT1A", "INSR", "PPARG", "GCK"],
    "microbiome": ["IDO1", "TDO2", "AHR", "NOS2", "FMO3"],
    "inflammatory": ["TNF", "IL6", "PTGS2", "NOS2", "CYBB", "NFE2L2"],
    "infectious": ["TLR4", "MYD88", "NOS2", "CYBB"],
    "cardiovascular": ["APOE", "PCSK9", "NOS3", "HMGCR"],
    "pulmonary": ["CYP1A1", "CYP2E1", "NOS2", "CYBB", "HIF1A", "GPX4"],
    "default": ["NFE2L2", "GPX4"],
}

# Category keyword → category (whole-word / stem match against normalized text)
CATEGORY_KEYWORDS: dict[str, tuple[str, ...]] = {
    "cancer": (
        "carcinoma",
        "adenocarcinoma",
        "sarcoma",
        "lymphoma",
        "leukemia",
        "leukaemia",
        "blastoma",
        "neoplasm",
        "malignant",
        "tumor",
        "tumour",
        "cancer",
        "melanoma",
        "myeloma",
    ),
    "neurological": (
        "neurodegenerat",
        "encephal",
        "neuropath",
        "parkinson",
        "alzheimer",
        "dementia",
        "epilep",
        "seizure",
        "demyelin",
        "ataxia",
        "myasthenia",
        "psychos",
        "schizophren",
        "depressi",
        "autism",
        "migraine",
    ),
    "metabolic": (
        "metabol",
        "ketoacid",
        "diabet",
        "hyperlipid",
        "dyslipid",
        "glycogen",
        "lysosomal",
        "urea cycle",
        "aminoacidopath",
        "mitochondri",
        "oxphos",
        "storage disease",
        "obesity",
        "nafld",
        "nash",
    ),
    "microbiome": (
        "dysbiosis",
        "clostrid",
        "helicobacter",
        "microbiome",
        "gut flora",
    ),
    "inflammatory": (
        "inflam",
        "arthritis",
        "colitis",
        "crohn",
        "lupus",
        "vasculitis",
        "psoriasis",
        "scleroderma",
    ),
    "infectious": (
        "infect",
        "viral",
        "bacteri",
        "fungal",
        "sepsis",
        "tubercul",
        "pneumonia",
        "hepatitis",
        "malaria",
        "hiv",
    ),
    "cardiovascular": (
        "cardio",
        "myocard",
        "atheroscler",
        "hypertens",
        "heart failure",
        "ischemi",
        "coronary",
        "arrhythmi",
        "stroke",
    ),
    "pulmonary": (
        "pulmon",
        "bronch",
        "asthma",
        "copd",
        "emphysema",
        "interstitial lung",
        "fibrosis lung",
        "pulmonary fibrosis",
        "respiratory",
        "sarcoid",
    ),
}

# Specific mechanism token cues (higher confidence than bare category keywords)
# (match_stems, category, extra_pathway_bias, site_hint, driver_genes)
TOKEN_CUES: list[tuple[tuple[str, ...], str, dict[str, float], str | None, list[str]]] = [
    (
        ("mitochondri", "oxphos", "leber hereditary", "lhond", "mt dna", "mtdna"),
        "metabolic",
        {
            "ketone_body_metabolism": 1.35,
            "fatty_acid_oxidation": 1.3,
            "brain_energy_metabolism": 1.25,
            "lipid_peroxidation": 1.2,
        },
        "brain",
        ["MT-ND1", "MT-ND4", "POLG", "SURF1"],
    ),
    (
        (
            "cholestat",
            "cholangio",
            "cholangitis",
            "biliary cirrhosis",
            "biliary cholang",
            "primary biliary",
            "sclerosing cholang",
        ),
        "metabolic",
        {
            "cytochrome_p450_detox": 1.45,
            "methionine_transsulfuration": 1.35,
            "urea_cycle": 1.25,
            "lipid_peroxidation": 1.2,
        },
        "liver",
        ["ABCB4", "NR1H4", "ATP8B1"],
    ),
    (
        (
            "neuroinflam",
            "microglia",
            "demyelin",
            "leukodystroph",
            "krabbe",
            "prion",
            "creutzfeldt",
        ),
        "neurological",
        {
            "neuroinflammation": 1.55,
            "lipid_peroxidation": 1.35,
            "brain_energy_metabolism": 1.25,
        },
        "brain",
        ["TNF", "IL6", "IDO1", "NOS2"],
    ),
    (
        ("ketotic", "ketogenesis", "ketoacidosis", "maple syrup", "branched chain"),
        "metabolic",
        {"ketone_body_metabolism": 1.7, "fatty_acid_oxidation": 1.4},
        "systemic",
        ["HMGCS2", "BCKDHA", "BCKDHB", "DBT"],
    ),
    (
        ("urea cycle", "hyperammon", "ornithine", "citrullin"),
        "metabolic",
        {"urea_cycle": 1.7, "ketone_body_metabolism": 1.2},
        "liver",
        ["OTC", "ASS1", "ASL", "CPS1"],
    ),
    (
        ("oxidative stress", "lipid peroxid", "ferroptos"),
        "inflammatory",
        {"lipid_peroxidation": 1.55, "apoptosis_necrosis": 1.25},
        None,
        ["GPX4", "ACSL4", "NFE2L2", "CYBB"],
    ),
    (
        ("gut dysbiosis", "putrefaction", "small intestinal bacterial"),
        "microbiome",
        {
            "gut_microbiome_fermentation": 1.65,
            "microbial_proteolysis_putrefaction": 1.55,
        },
        "intestine",
        ["FMO3", "IDO1", "TDO2"],
    ),
    (
        ("warburg", "glycolytic tumor", "hypoxia inducible"),
        "cancer",
        {"glycolysis_warburg": 1.55, "Kras_mapk_proliferation": 1.3},
        None,
        ["HK2", "LDHA", "HIF1A", "KRAS"],
    ),
    (
        ("cystic fibrosis", "bronchiectasis", "airway inflam"),
        "pulmonary",
        {"lipid_peroxidation": 1.4, "cytochrome_p450_detox": 1.2},
        "lung",
        ["CFTR", "CYP1A1", "NOS2"],
    ),
    (
        ("atheroscler", "coronary artery", "myocardial infarct"),
        "cardiovascular",
        {
            "lipid_peroxidation": 1.35,
            "fatty_acid_oxidation": 1.25,
            "mevalonate_cholesterol": 1.25,
        },
        "heart",
        ["APOE", "PCSK9", "LDLR"],
    ),
    (
        ("phenylketonuria", "pku", "phenylalanine"),
        "metabolic",
        {"one_carbon_folate": 1.3, "neurotransmitter_metabolism": 1.25, "urea_cycle": 1.15},
        "systemic",
        ["PAH", "QDPR"],
    ),
    (
        ("wilson", "copper overload", "hepatolenticular"),
        "metabolic",
        {
            "cytochrome_p450_detox": 1.35,
            "lipid_peroxidation": 1.3,
            "methionine_transsulfuration": 1.2,
        },
        "liver",
        ["ATP7B"],
    ),
    (
        ("fabry", "gaucher", "niemann pick", "pompe", "lysosomal storage"),
        "metabolic",
        {"fatty_acid_oxidation": 1.3, "lipid_peroxidation": 1.25, "apoptosis_necrosis": 1.15},
        "systemic",
        ["GLA", "GBA", "NPC1", "GAA"],
    ),
    (
        ("cystinosis", "nephropathic", "glomerulopath", "nephritis"),
        "metabolic",
        {"urea_cycle": 1.35, "lipid_peroxidation": 1.2},
        "kidney",
        ["CTNS", "UMOD"],
    ),
    (
        ("pancreatit", "exocrine pancreas"),
        "inflammatory",
        {"lipid_peroxidation": 1.3, "apoptosis_necrosis": 1.25, "glycolysis_warburg": 1.15},
        "pancreas",
        ["PRSS1", "SPINK1", "CFTR"],
    ),
    (
        ("hemoly", "g6pd", "sickle", "thalassem"),
        "inflammatory",
        {"lipid_peroxidation": 1.5, "apoptosis_necrosis": 1.35, "glycolysis_warburg": 1.15},
        "blood",
        ["HBB", "G6PD", "CYBB", "NFE2L2"],
    ),
    (
        ("thyroidit", "hypothyroid", "hyperthyroid", "hashimoto", "graves"),
        "metabolic",
        {
            "fatty_acid_oxidation": 1.3,
            "ketone_body_metabolism": 1.25,
            "mevalonate_cholesterol": 1.2,
            "lipid_peroxidation": 1.15,
        },
        "endocrine gland",
        ["TPO", "TG", "TSHR"],
    ),
    (
        ("addison", "adrenal insuffici", "cushing", "cortisol"),
        "metabolic",
        {"ketone_body_metabolism": 1.4, "fatty_acid_oxidation": 1.3, "glycolysis_warburg": 1.15},
        "adrenal gland",
        ["CYP21A2", "MC2R", "NR3C1"],
    ),
    (
        ("uremi", "uremic", "end stage renal", "stage 5 ckd"),
        "metabolic",
        {
            "urea_cycle": 1.6,
            "microbial_proteolysis_putrefaction": 1.35,
            "gut_microbiome_fermentation": 1.25,
        },
        "kidney",
        ["UMOD", "PKD1"],
    ),
    (
        ("septic", "sepsis", "endotox"),
        "infectious",
        {
            "lipid_peroxidation": 1.4,
            "apoptosis_necrosis": 1.35,
            "glycolysis_warburg": 1.25,
            "urea_cycle": 1.15,
        },
        "systemic",
        ["TLR4", "MYD88", "TNF", "IL6"],
    ),
    (
        ("steatohepat", "hepatic steatos", "alcohol use", "alcoholic liver"),
        "metabolic",
        {
            "fatty_acid_oxidation": 1.4,
            "lipid_peroxidation": 1.35,
            "cytochrome_p450_detox": 1.35,
            "ketone_body_metabolism": 1.2,
        },
        "liver",
        ["PNPLA3", "CYP2E1", "ADH1B"],
    ),
    (
        ("carcinoid", "neuroendocrine tumor", "pheochromocytoma", "paraganglioma"),
        "neurological",
        {
            "neurotransmitter_metabolism": 1.5,
            "gut_microbiome_fermentation": 1.3,
            "lipid_peroxidation": 1.2,
        },
        "intestine",
        ["TPH1", "DDC", "CHGA"],
    ),
    (
        ("porphyria", "heme biosynth"),
        "metabolic",
        {
            "cytochrome_p450_detox": 1.4,
            "one_carbon_folate": 1.25,
            "urea_cycle": 1.2,
            "lipid_peroxidation": 1.2,
        },
        "liver",
        ["HMBS", "PPOX", "ALAD"],
    ),
    (
        ("toxicity", "toxic metabol", "lactic acidosis", "metformin"),
        "metabolic",
        {"ketone_body_metabolism": 1.45, "urea_cycle": 1.3, "glycolysis_warburg": 1.25},
        "systemic",
        ["HMGCS2", "LDHA", "OTC"],
    ),
    (
        ("scurvy", "ascorbate defici", "vitamin c defici"),
        "inflammatory",
        {"lipid_peroxidation": 1.45, "apoptosis_necrosis": 1.3},
        "systemic",
        ["SOD2", "GPX4", "CYBB"],
    ),
    (
        ("beriberi", "thiamine defici", "wernicke"),
        "metabolic",
        {
            "brain_energy_metabolism": 1.45,
            "glycolysis_warburg": 1.3,
            "lipid_peroxidation": 1.2,
        },
        "heart",
        ["PDHA1", "DLD", "OGDH"],
    ),
    (
        ("preeclamp", "eclampsia", "endometriosis"),
        "inflammatory",
        {"lipid_peroxidation": 1.4, "apoptosis_necrosis": 1.25},
        "uterus",
        ["VEGFA", "TNF", "PTGS2"],
    ),
    (
        ("pulmonary fibrosis", "interstitial lung", "idiopathic pulmonary"),
        "pulmonary",
        {"lipid_peroxidation": 1.45, "cytochrome_p450_detox": 1.25, "glycolysis_warburg": 1.2},
        "lung",
        ["MUC5B", "TERT", "SFTPC"],
    ),
    (
        ("long covid", "post covid", "dysautonomia"),
        "infectious",
        {
            "lipid_peroxidation": 1.3,
            "neuroinflammation": 1.25,
            "brain_energy_metabolism": 1.2,
            "glycolysis_warburg": 1.15,
        },
        "systemic",
        ["ACE2", "TNF", "IL6"],
    ),
    (
        ("vasculitis", "granulomatosis", "wegener", "behcet", "behçet"),
        "inflammatory",
        {"lipid_peroxidation": 1.4, "apoptosis_necrosis": 1.3},
        "systemic",
        ["PRTN3", "MPO", "TNF"],
    ),
]

# Anatomic tokens in disease text → default_site
SITE_CUES: list[tuple[tuple[str, ...], str]] = [
    (("lung", "pulmon", "bronch", "airway", "alveol"), "lung"),
    (("brain", "cerebr", "cortical", "optic neuropath", "encephal"), "brain"),
    (("liver", "hepat", "biliary", "cholang", "cirrhos"), "liver"),
    (("kidney", "renal", "nephro", "glomerul"), "kidney"),
    (("heart", "cardio", "myocard", "coronary"), "heart"),
    (("colon", "intestin", "bowel", "gut", "enteric", "ileum"), "intestine"),
    (("pancrea",), "pancreas"),
    (("breast", "mammary"), "breast"),
    (("prostate",), "prostate gland"),
    (("skin", "dermat", "cutane"), "skin of body"),
    (("stomach", "gastric"), "stomach"),
    (("ovary", "ovarian"), "ovary"),
    (("bladder",), "bladder organ"),
    (("adipose", "obesity", "adipocyte"), "adipose tissue"),
    (("eye", "optic", "retina", "ocular"), "eye"),
    (("thyroid",), "endocrine gland"),
]

# Bare tokens too generic for mechanism activation alone (keep near-healthy)
GENERIC_DISEASE_TOKENS = frozenset(
    {
        "cancer",
        "tumor",
        "tumour",
        "carcinoma",
        "neoplasm",
        "disease",
        "syndrome",
        "disorder",
        "infection",
        "illness",
        "condition",
    }
)


def merge_pathway_bias(*biases: dict[str, float] | None) -> dict[str, float]:
    out: dict[str, float] = {}
    for bias in biases:
        if not bias:
            continue
        for pid, val in bias.items():
            out[pid] = max(float(out.get(pid, 0.0)), float(val))
    return out


def template_for_category(category: str | None) -> dict[str, float]:
    cat = (category or "default").lower()
    if cat == "neurodegenerative":
        cat = "neurological"
    return dict(CATEGORY_PATHWAY_TEMPLATES.get(cat) or CATEGORY_PATHWAY_TEMPLATES["default"])


def infer_site_from_text(text: str) -> str | None:
    t = text.lower()
    for stems, site in SITE_CUES:
        if any(s in t for s in stems):
            return site
    return None


def match_category_keywords(text: str) -> list[str]:
    t = text.lower()
    hits: list[str] = []
    for cat, stems in CATEGORY_KEYWORDS.items():
        if any(s in t for s in stems):
            hits.append(cat)
    return hits


def match_token_cues(text: str) -> list[dict[str, Any]]:
    t = text.lower()
    hits: list[dict[str, Any]] = []
    for stems, cat, bias, site, genes in TOKEN_CUES:
        matched = [s for s in stems if s in t]
        if not matched:
            continue
        # Fix gene typos with spaces (ABC B4 → ABCB4)
        clean_genes = [g.replace(" ", "") for g in genes]
        hits.append(
            {
                "matched": matched,
                "category": cat,
                "pathway_bias": dict(bias),
                "default_site": site,
                "driver_genes": clean_genes,
                "strength": 0.55 + 0.05 * min(len(matched), 3),
            }
        )
    return hits
