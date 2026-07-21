from __future__ import annotations

import os
from pathlib import Path


def _resolve_data_dir() -> Path:
    """
    Resolve the atlas data directory for both editable repo installs and wheels.

    Priority:
      1. EXHALEPATH_DATA / VOC_DATA env override
      2. Packaged data next to this module (``exhalepath/data``) — pip wheel
      3. Repo layout ``exhalepath_atlas/data`` (editable / source checkout)
    """
    for key in ("EXHALEPATH_DATA", "VOC_DATA"):
        env = os.environ.get(key)
        if env:
            p = Path(env).expanduser().resolve()
            if p.exists():
                return p

    pkg_data = Path(__file__).resolve().parent / "data"
    if (pkg_data / "knowledge" / "voc_catalog.json").exists():
        return pkg_data

    # src/exhalepath/config.py → parents[2] == exhalepath_atlas/
    repo_data = Path(__file__).resolve().parents[2] / "data"
    if (repo_data / "knowledge" / "voc_catalog.json").exists():
        return repo_data

    # Last resort: create writable user cache and point there
    user = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share")) / "voc-breath"
    user.mkdir(parents=True, exist_ok=True)
    return user


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = _resolve_data_dir()
KNOWLEDGE_DIR = DATA_DIR / "knowledge"
CACHE_DIR = DATA_DIR / "cache"
PROCESSED_DIR = DATA_DIR / "processed"
MODELS_DIR = DATA_DIR / "models"

GDC_API = "https://api.gdc.cancer.gov"
REACTOME_CONTENT = "https://reactome.org/ContentService"
OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_DIR = DATA_DIR / "chembl"

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

USER_AGENT = "voc-breath/1.1 (ExhalePath Atlas; research exhaled VOC prediction)"
