"""Priority 10 — BindingDB / pharmacology complement to ChEMBL."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# High-value target pharmacology notes for VOC-pathway enzymes
TARGET_PHARM = [
    {"gene": "HMGCR", "ligand_class": "statins", "effect_on_voc": {"isoprene": "down"}, "note": "Mevalonate block lowers isoprene"},
    {"gene": "CYP2E1", "ligand_class": "ethanol/solvent substrates", "effect_on_voc": {"acetaldehyde": "up"}, "note": "Inducible ethanol oxidation"},
    {"gene": "IDO1", "ligand_class": "IDO inhibitors", "effect_on_voc": {"indole": "down"}, "note": "Tryptophan–kynurenine / indole axis"},
    {"gene": "PTGS2", "ligand_class": "NSAIDs/COX2i", "effect_on_voc": {"pentane": "down", "hexanal": "mod"}, "note": "Inflammatory peroxidation coupling"},
    {"gene": "ALOX15", "ligand_class": "LOX inhibitors", "effect_on_voc": {"hexanal": "down"}, "note": "PUFA peroxidation aldehydes"},
    {"gene": "FMO3", "ligand_class": "TMA substrates", "effect_on_voc": {"trimethylamine": "mod"}, "note": "TMA→TMAO vs breath TMA"},
    {"gene": "KRAS", "ligand_class": "KRAS G12C inhibitors", "effect_on_voc": {"acetaldehyde": "down", "ethanol": "down"}, "note": "Tumor Warburg coupling"},
]


class BindingDBSource(DataSource):
    priority = 10
    key = "bindingdb"
    title = "BindingDB / pharmacology target notes"
    description = "Ligand-class effects on VOC-pathway enzymes (ChEMBL complement)"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        doc = {
            "version": "1.0.0",
            "ref": "https://www.bindingdb.org/",
            "targets": TARGET_PHARM,
            "chembl_bridge": "Use exhalepath harvest-chembl for million-scale pChEMBL rows",
        }
        path = self.write_json("binding_pharmacology.json", doc)
        return {"pharmacology": path, "manifest": self.write_manifest(n_targets=len(TARGET_PHARM))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "binding_pharmacology.json").read_text())
        out = knowledge_dir / "datasource_bindingdb.json"
        out.write_text(json.dumps(doc, indent=2))
        return {"path": str(out), "n_targets": len(doc["targets"])}
