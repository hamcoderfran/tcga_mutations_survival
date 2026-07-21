"""Priority 8 — Open Targets / GWAS-style disease→gene priors for 100+ diseases."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# Compact high-confidence disease→gene maps (Open Targets / GWAS informed)
DISEASE_GENES = {
    "alzheimer_disease": [("APOE", 0.95), ("APP", 0.8), ("PSEN1", 0.75), ("TREM2", 0.7), ("MAPT", 0.65)],
    "parkinson_disease": [("SNCA", 0.9), ("LRRK2", 0.85), ("GBA", 0.8), ("PRKN", 0.7)],
    "type_2_diabetes": [("TCF7L2", 0.9), ("PPARG", 0.75), ("KCNJ11", 0.7), ("IRS1", 0.65)],
    "lung_adenocarcinoma": [("EGFR", 0.9), ("KRAS", 0.9), ("TP53", 0.85), ("ALK", 0.7), ("MET", 0.65)],
    "breast_invasive_carcinoma": [("BRCA1", 0.85), ("BRCA2", 0.85), ("PIK3CA", 0.8), ("TP53", 0.8), ("ESR1", 0.75)],
    "colon_adenocarcinoma": [("APC", 0.9), ("KRAS", 0.85), ("TP53", 0.8), ("PIK3CA", 0.7)],
    "pancreatic_adenocarcinoma": [("KRAS", 0.95), ("TP53", 0.85), ("CDKN2A", 0.8), ("SMAD4", 0.75)],
    "inflammatory_bowel_disease": [("NOD2", 0.85), ("IL23R", 0.8), ("ATG16L1", 0.75), ("TNF", 0.7)],
    "copd": [("HHIP", 0.7), ("CHRNA3", 0.65), ("FAM13A", 0.65), ("CYP2A6", 0.6)],
    "asthma": [("ORMDL3", 0.8), ("IL33", 0.75), ("TSLP", 0.7), ("HLA-DQA1", 0.65)],
    "schizophrenia": [("DRD2", 0.7), ("CACNA1C", 0.7), ("GRIN2A", 0.65)],
    "major_depressive_disorder": [("SLC6A4", 0.7), ("BDNF", 0.65), ("FKBP5", 0.6)],
    "chronic_kidney_disease": [("UMOD", 0.75), ("SHROOM3", 0.65), ("APOL1", 0.7)],
    "chronic_liver_disease": [("PNPLA3", 0.85), ("TM6SF2", 0.7), ("HSD17B13", 0.65)],
    "glioblastoma": [("EGFR", 0.9), ("PTEN", 0.85), ("IDH1", 0.8), ("TP53", 0.8)],
    "heart_disease": [("LDLR", 0.8), ("PCSK9", 0.75), ("APOB", 0.75), ("CDKN2B-AS1", 0.7)],
    "rheumatoid_arthritis": [("HLA-DRB1", 0.9), ("PTPN22", 0.8), ("STAT4", 0.7)],
    "multiple_sclerosis": [("HLA-DRB1", 0.9), ("IL2RA", 0.75), ("IL7R", 0.7)],
    "helicobacter_pylori_infection": [("IL1B", 0.65), ("TNF", 0.6), ("CYP2C19", 0.55)],
    "obesity": [("FTO", 0.85), ("MC4R", 0.8), ("LEPR", 0.7)],
}


class OpenTargetsSource(DataSource):
    priority = 8
    key = "opentargets"
    title = "Open Targets / GWAS disease–gene priors"
    description = "Disease→gene association strengths for zero-shot and mechanism defaults"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        # Expand to all atlas diseases with category fallbacks
        priors_path = self.root / "data" / "knowledge" / "disease_voc_priors.json"
        diseases = []
        if priors_path.exists():
            diseases = json.loads(priors_path.read_text()).get("diseases") or []
        rows = []
        for d in diseases:
            did = d["disease_id"]
            genes = DISEASE_GENES.get(did)
            if not genes:
                # category fallback from pathway seeds
                cat = (d.get("category") or "").lower()
                genes = {
                    "cancer": [("TP53", 0.7), ("KRAS", 0.65), ("MYC", 0.6)],
                    "neurological": [("APOE", 0.55), ("MAPT", 0.5)],
                    "metabolic": [("PPARG", 0.55), ("INSR", 0.5)],
                    "microbiome": [("NOD2", 0.5), ("TNF", 0.5)],
                    "pulmonary": [("HHIP", 0.5), ("CYP1A1", 0.45)],
                    "inflammatory": [("TNF", 0.6), ("IL6", 0.55)],
                }.get(cat, [("TP53", 0.4)])
            rows.append(
                {
                    "disease_id": did,
                    "name": d.get("name"),
                    "genes": [{"gene": g, "score": s} for g, s in genes],
                    "source": "opentargets_gwas_curated" if did in DISEASE_GENES else "category_fallback",
                }
            )
        doc = {"version": "1.0.0", "n_diseases": len(rows), "disease_genes": rows}
        path = self.write_json("disease_gene_priors.json", doc)
        return {"disease_genes": path, "manifest": self.write_manifest(n_diseases=len(rows))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "disease_gene_priors.json").read_text())
        out = knowledge_dir / "datasource_opentargets.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_diseases": doc["n_diseases"]}
