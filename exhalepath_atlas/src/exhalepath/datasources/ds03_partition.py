"""Priority 3 — Blood:air partition / Henry’s-law λ tables for Farhi transport."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .base import DataSource

# Literature / Poulin–Krishnan–style λ blood:air for ExhalePath VOCs
# Values are experimental or QSPR-typical; used to overwrite physio defaults.
LAMBDA_TABLE = {
    "acetone": 340.0,
    "acetaldehyde": 190.0,
    "ethanol": 1800.0,
    "methanol": 2500.0,
    "isoprene": 0.75,
    "pentane": 0.4,
    "ethane": 0.2,
    "butane": 0.5,
    "hexane": 0.6,
    "octane": 0.8,
    "hexanal": 50.0,
    "heptanal": 60.0,
    "nonanal": 80.0,
    "decanal": 90.0,
    "octanal": 70.0,
    "ammonia": 1000.0,
    "dms": 10.0,
    "dmts": 15.0,
    "dimethyl_disulfide": 12.0,
    "hydrogen_sulfide": 5.0,
    "indole": 120.0,
    "phenol": 200.0,
    "toluene": 15.0,
    "benzene": 8.0,
    "limonene": 5.0,
    "formaldehyde": 3000.0,
    "2_butanone": 120.0,
    "propanol": 900.0,
    "isopropanol": 800.0,
    "styrene": 20.0,
    "ethylbenzene": 18.0,
    "xylene": 20.0,
    "cyclohexane": 1.2,
    "methyl_acetate": 100.0,
    "ethyl_acetate": 120.0,
    "acetonitrile": 300.0,
    "furan": 5.0,
    "2_pentanone": 90.0,
    "3_methylbutanal": 40.0,
    "benzaldehyde": 80.0,
    "undecane": 1.5,
    "dodecane": 2.0,
    "carbon_disulfide": 8.0,
    "methyl_mercaptan": 6.0,
    "trimethylamine": 50.0,
    "dimethyl_amine": 60.0,
    "pyrrole": 30.0,
    "allyl_methyl_sulfide": 10.0,
    "propionaldehyde": 35.0,
    "crotonaldehyde": 40.0,
}


class PartitionLambdaSource(DataSource):
    priority = 3
    key = "partition"
    title = "Blood:air partition coefficients (λ)"
    description = "Farhi alveolar-release λ table (experimental + QSPR-typical)"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        rows = []
        for voc_id, lam in LAMBDA_TABLE.items():
            rows.append(
                {
                    "voc_id": voc_id,
                    "lambda_blood_air": lam,
                    "log10_lambda": float(np.log10(max(lam, 1e-9))),
                    "source": "curated_literature_qspr",
                    "ref": "DOI:10.1088/1752-7155/10/1/017103",
                }
            )
        doc = {
            "version": "1.0.0",
            "n_vocs": len(rows),
            "description": self.description,
            "partitions": rows,
        }
        path = self.write_json("blood_air_partition.json", doc)
        return {"partitions": path, "manifest": self.write_manifest(n_vocs=len(rows))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "blood_air_partition.json").read_text())
        out = knowledge_dir / "datasource_partition.json"
        out.write_text(json.dumps(doc, indent=2))
        phys_path = knowledge_dir / "physio_constants.json"
        if phys_path.exists():
            phys = json.loads(phys_path.read_text())
            lam = phys.setdefault("blood_air_partition_lambda", {})
            n = 0
            for row in doc["partitions"]:
                lam[row["voc_id"]] = float(row["lambda_blood_air"])
                n += 1
            phys["partition_source"] = "datasource_priority3_lambda_table"
            phys_path.write_text(json.dumps(phys, indent=2))
            return {"path": str(out), "n_lambda_updated": n}
        return {"path": str(out)}
