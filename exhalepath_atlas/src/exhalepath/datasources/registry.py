"""Registry and integrator for priority 1–14 datasources."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .ds01_metabolomics import MetabolomicsReposSource
from .ds02_hmdb import HMDBSource
from .ds03_partition import PartitionLambdaSource
from .ds04_mvoc import MVOCSource
from .ds05_pubchem import PubChemSource
from .ds06_reactome import ReactomeSource
from .ds07_gtex import GTExSource
from .ds08_opentargets import OpenTargetsSource
from .ds09_gdc import GDCSource
from .ds10_bindingdb import BindingDBSource
from .ds11_blood_proxy import BloodProxySource
from .ds12_nist import NISTSource
from .ds13_hbdb import HBDBSource
from .ds14_kegg import KEGGSource

SOURCE_CLASSES = [
    MetabolomicsReposSource,
    HMDBSource,
    PartitionLambdaSource,
    MVOCSource,
    PubChemSource,
    ReactomeSource,
    GTExSource,
    OpenTargetsSource,
    GDCSource,
    BindingDBSource,
    BloodProxySource,
    NISTSource,
    HBDBSource,
    KEGGSource,
]


def DATASOURCE_REGISTRY(root: Path) -> list:
    return [cls(root) for cls in SOURCE_CLASSES]


def integrate_all_datasources(
    *,
    root: Path | None = None,
    offline: bool = False,
    priorities: list[int] | None = None,
) -> dict[str, Any]:
    """
    Harvest + fuse priority 1–14 datasources into exhalepath_atlas knowledge/.

    offline=True skips live HTTP where possible (still writes curated tables).
    """
    root = Path(root or Path(__file__).resolve().parents[3])
    knowledge_dir = root / "data" / "knowledge"
    knowledge_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for src in DATASOURCE_REGISTRY(root):
        if priorities and src.priority not in priorities:
            continue
        harvest_paths = src.harvest(offline=offline)
        fuse_info = src.fuse(knowledge_dir)
        results.append(
            {
                "priority": src.priority,
                "key": src.key,
                "title": src.title,
                "harvest": {k: str(v) for k, v in harvest_paths.items()},
                "fuse": fuse_info,
            }
        )

    # Master integration manifest
    manifest = {
        "version": "1.0.0",
        "product": "exhalepath_atlas",
        "n_datasources": len(results),
        "offline": offline,
        "datasources": results,
        "note": (
            "ChEMBL chemogenomic harvest remains available via "
            "`python -m exhalepath harvest-chembl` (already bundled in this atlas tree)."
        ),
    }
    man_path = root / "data" / "datasources" / "INTEGRATION_MANIFEST.json"
    man_path.parent.mkdir(parents=True, exist_ok=True)
    man_path.write_text(json.dumps(manifest, indent=2))

    # Capability card for 100+ diseases
    priors = knowledge_dir / "disease_voc_priors.json"
    n_diseases = 0
    if priors.exists():
        n_diseases = len(json.loads(priors.read_text()).get("diseases") or [])
    capability = {
        "n_diseases": n_diseases,
        "n_datasources_integrated": len(results),
        "priorities": [r["priority"] for r in results],
        "knowledge_fragments": sorted(
            p.name for p in knowledge_dir.glob("datasource_*.json")
        ),
    }
    cap_path = knowledge_dir / "atlas_capability.json"
    cap_path.write_text(json.dumps(capability, indent=2))
    manifest["capability"] = capability
    man_path.write_text(json.dumps(manifest, indent=2))
    return manifest
