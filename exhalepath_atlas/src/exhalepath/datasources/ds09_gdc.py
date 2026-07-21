"""Priority 9 — TCGA/GDC molecular feature summaries for cancer diseases."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# Representative frequent drivers by TCGA project (public knowledge summaries)
TCGA_DRIVERS = {
    "lung_adenocarcinoma": {"project": "TCGA-LUAD", "drivers": ["KRAS", "TP53", "STK11", "KEAP1", "EGFR"]},
    "breast_invasive_carcinoma": {"project": "TCGA-BRCA", "drivers": ["PIK3CA", "TP53", "GATA3", "CDH1", "MAP3K1"]},
    "colon_adenocarcinoma": {"project": "TCGA-COAD", "drivers": ["APC", "TP53", "KRAS", "PIK3CA", "SMAD4"]},
    "pancreatic_adenocarcinoma": {"project": "TCGA-PAAD", "drivers": ["KRAS", "TP53", "CDKN2A", "SMAD4"]},
    "hepatocellular_carcinoma": {"project": "TCGA-LIHC", "drivers": ["TP53", "CTNNB1", "ALB", "AXIN1"]},
    "ovarian_cancer": {"project": "TCGA-OV", "drivers": ["TP53", "BRCA1", "BRCA2", "CDK12"]},
    "glioblastoma": {"project": "TCGA-GBM", "drivers": ["EGFR", "PTEN", "TP53", "CDKN2A"]},
    "cancer_prostate": {"project": "TCGA-PRAD", "drivers": ["TP53", "SPOP", "FOXA1", "PTEN"]},
    "cancer_kidney": {"project": "TCGA-KIRC", "drivers": ["VHL", "PBRM1", "SETD2", "BAP1"]},
    "cancer_stomach": {"project": "TCGA-STAD", "drivers": ["TP53", "ARID1A", "CDH1", "PIK3CA"]},
    "cancer_skin_melanoma": {"project": "TCGA-SKCM", "drivers": ["BRAF", "NRAS", "TP53", "CDKN2A"]},
    "cancer_bladder": {"project": "TCGA-BLCA", "drivers": ["TP53", "FGFR3", "RB1", "KDM6A"]},
    "esophageal_cancer": {"project": "TCGA-ESCA", "drivers": ["TP53", "CDKN2A", "NFE2L2"]},
    "head_neck_cancer": {"project": "TCGA-HNSC", "drivers": ["TP53", "CDKN2A", "PIK3CA", "NOTCH1"]},
    "cervical_cancer": {"project": "TCGA-CESC", "drivers": ["PIK3CA", "EP300", "FBXW7"]},
}


class GDCSource(DataSource):
    priority = 9
    key = "gdc"
    title = "TCGA/GDC molecular drivers"
    description = "Cancer project driver-gene panels for mutation-conditioned VOC prediction"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        rows = []
        for did, meta in TCGA_DRIVERS.items():
            rows.append(
                {
                    "disease_id": did,
                    "gdc_project": meta["project"],
                    "driver_genes": meta["drivers"],
                    "gdc_portal": f"https://portal.gdc.cancer.gov/projects/{meta['project']}",
                }
            )
        doc = {"version": "1.0.0", "n_projects": len(rows), "projects": rows}
        path = self.write_json("gdc_disease_drivers.json", doc)
        return {"drivers": path, "manifest": self.write_manifest(n_projects=len(rows))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "gdc_disease_drivers.json").read_text())
        out = knowledge_dir / "datasource_gdc.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_projects": doc["n_projects"]}
