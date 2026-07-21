from __future__ import annotations

from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PACKAGE_ROOT / "data"
KNOWLEDGE_DIR = DATA_DIR / "knowledge"
CACHE_DIR = DATA_DIR / "cache"
PROCESSED_DIR = DATA_DIR / "processed"
MODELS_DIR = DATA_DIR / "models"

GDC_API = "https://api.gdc.cancer.gov"
REACTOME_CONTENT = "https://reactome.org/ContentService"
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"

# Broad TCGA / TARGET cancer project panel for multi-million mutation rows
DEFAULT_GDC_PROJECTS = [
    "TCGA-LUAD",
    "TCGA-LUSC",
    "TCGA-BRCA",
    "TCGA-PAAD",
    "TCGA-COAD",
    "TCGA-READ",
    "TCGA-LIHC",
    "TCGA-STAD",
    "TCGA-ESCA",
    "TCGA-OV",
    "TCGA-UCEC",
    "TCGA-CESC",
    "TCGA-PRAD",
    "TCGA-KIRC",
    "TCGA-KIRP",
    "TCGA-KICH",
    "TCGA-BLCA",
    "TCGA-THCA",
    "TCGA-GBM",
    "TCGA-LGG",
    "TCGA-SKCM",
    "TCGA-HNSC",
    "TCGA-SARC",
    "TCGA-LAML",
    "TCGA-DLBC",
    "TCGA-THYM",
    "TCGA-MESO",
    "TCGA-UVM",
    "TCGA-ACC",
    "TCGA-PCPG",
    "TCGA-TGCT",
    "TCGA-CHOL",
]

USER_AGENT = "ExhalePath/0.1 (research; pathway-VOC breath modeling)"
