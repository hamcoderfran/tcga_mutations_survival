"""Priority 7 — GTEx-like tissue expression priors for VOC-pathway enzymes."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .base import DataSource

# Relative expression priors (0–1) by tissue for key VOC enzymes (GTEx-informed literature)
ENZYME_TISSUE = {
    "HMGCS2": {"liver": 1.0, "kidney": 0.35, "colon": 0.15, "lung": 0.05, "brain": 0.05},
    "HMGCL": {"liver": 0.95, "kidney": 0.4, "brain": 0.1},
    "CPT1A": {"liver": 0.9, "heart": 0.7, "muscle": 0.65, "adipose": 0.5},
    "LDHA": {"tumor": 0.95, "lung": 0.55, "muscle": 0.7, "brain": 0.4},
    "HK2": {"tumor": 0.9, "lung": 0.5, "brain": 0.35},
    "ALOX15": {"lung": 0.7, "blood": 0.6, "brain": 0.3},
    "CYP2E1": {"liver": 1.0, "lung": 0.35, "kidney": 0.25},
    "CYP1A1": {"lung": 0.8, "liver": 0.5},
    "IDO1": {"gut": 0.7, "lung": 0.4, "brain": 0.35, "blood": 0.5},
    "GPX4": {"brain": 0.7, "liver": 0.6, "lung": 0.55, "tumor": 0.5},
    "FMO3": {"liver": 1.0, "lung": 0.2},
    "CBS": {"liver": 0.85, "brain": 0.5, "kidney": 0.4},
    "HMGCR": {"liver": 0.9, "brain": 0.35, "muscle": 0.3},
}


class GTExSource(DataSource):
    priority = 7
    key = "gtex"
    title = "GTEx tissue expression priors"
    description = "Enzyme×tissue expression priors for VOC biosynthetic localization"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        # offline-capable curated matrix; optional live GTEx median query skipped (auth/size)
        genes = sorted(ENZYME_TISSUE)
        tissues = sorted({t for m in ENZYME_TISSUE.values() for t in m})
        matrix = []
        for g in genes:
            for t in tissues:
                matrix.append(
                    {
                        "gene": g,
                        "tissue": t,
                        "rel_expression": float(ENZYME_TISSUE[g].get(t, 0.05)),
                        "source": "gtex_informed_curated",
                    }
                )
        doc = {
            "version": "1.0.0",
            "n_genes": len(genes),
            "n_tissues": len(tissues),
            "expression": matrix,
            "note": "Relative units 0–1; live GTEx bulk download optional later",
        }
        path = self.write_json("gtex_enzyme_tissue.json", doc)
        return {"expression": path, "manifest": self.write_manifest(n_genes=len(genes))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "gtex_enzyme_tissue.json").read_text())
        out = knowledge_dir / "datasource_gtex.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_rows": len(doc["expression"])}
