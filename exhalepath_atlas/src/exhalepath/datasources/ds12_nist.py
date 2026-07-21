"""Priority 12 — NIST/GC-MS identification aids (RI / spectral metadata)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# Kovats RI (nonpolar) approximate literature values for identification QC
VOC_RI = {
    "acetone": 500,
    "ethanol": 500,
    "isoprene": 520,
    "pentane": 500,
    "hexane": 600,
    "benzene": 650,
    "toluene": 760,
    "ethylbenzene": 850,
    "xylene": 860,
    "styrene": 880,
    "hexanal": 800,
    "heptanal": 900,
    "octanal": 1000,
    "nonanal": 1100,
    "decanal": 1200,
    "benzaldehyde": 960,
    "limonene": 1030,
    "indole": 1290,
    "phenol": 980,
    "2_butanone": 580,
    "isopropanol": 500,
    "propanol": 550,
    "furan": 500,
    "trimethylamine": 450,
}


class NISTSource(DataSource):
    priority = 12
    key = "nist"
    title = "NIST/GC-MS identification metadata"
    description = "Retention-index aids for mapping public GC-MS peak tables → VOC IDs"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        rows = [
            {
                "voc_id": vid,
                "ri_nonpolar_approx": ri,
                "id_method": "kovats_ri_literature",
                "note": "Full NIST spectral library is licensed; RI used for peak mapping QC",
            }
            for vid, ri in VOC_RI.items()
        ]
        doc = {"version": "1.0.0", "n_vocs": len(rows), "retention_indices": rows}
        path = self.write_json("gcms_retention_indices.json", doc)
        return {"ri": path, "manifest": self.write_manifest(n_vocs=len(rows))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "gcms_retention_indices.json").read_text())
        out = knowledge_dir / "datasource_nist.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_vocs": doc["n_vocs"]}
