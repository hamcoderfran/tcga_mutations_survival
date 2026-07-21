"""Priority 11 — Blood metabolome proxy (UKB-like) for systemic VOC precursors."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# Public-knowledge blood metabolite ↔ breath VOC links (not UKB individual data)
BLOOD_BREATH_LINKS = [
    {"blood_metabolite": "3-hydroxybutyrate", "breath_voc": "acetone", "diseases": ["type_2_diabetes", "obesity"], "direction": "up"},
    {"blood_metabolite": "glucose", "breath_voc": "acetone", "diseases": ["type_2_diabetes"], "direction": "up"},
    {"blood_metabolite": "TMAO", "breath_voc": "trimethylamine", "diseases": ["chronic_kidney_disease", "heart_disease"], "direction": "up"},
    {"blood_metabolite": "indoxyl_sulfate", "breath_voc": "indole", "diseases": ["chronic_kidney_disease", "gut_dysbiosis"], "direction": "up"},
    {"blood_metabolite": "p-cresol_sulfate", "breath_voc": "phenol", "diseases": ["chronic_kidney_disease", "gut_dysbiosis"], "direction": "up"},
    {"blood_metabolite": "lactate", "breath_voc": "acetaldehyde", "diseases": ["lung_adenocarcinoma", "sepsis"], "direction": "up"},
    {"blood_metabolite": "ammonia_blood", "breath_voc": "ammonia", "diseases": ["chronic_liver_disease"], "direction": "up"},
    {"blood_metabolite": "isoprene_blood", "breath_voc": "isoprene", "diseases": ["major_depressive_disorder", "type_2_diabetes"], "direction": "mod"},
]


class BloodProxySource(DataSource):
    priority = 11
    key = "blood_proxy"
    title = "Blood metabolome → breath VOC proxy links"
    description = (
        "UK Biobank-style systemic metabolite bridges (public pathway links; "
        "no restricted UKB microdata)"
    )

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        doc = {
            "version": "1.0.0",
            "note": "Individual-level UKB requires application; this encodes public blood↔breath bridges",
            "links": BLOOD_BREATH_LINKS,
        }
        path = self.write_json("blood_breath_links.json", doc)
        return {"links": path, "manifest": self.write_manifest(n_links=len(BLOOD_BREATH_LINKS))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "blood_breath_links.json").read_text())
        out = knowledge_dir / "datasource_blood_proxy.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_links": len(doc["links"])}
