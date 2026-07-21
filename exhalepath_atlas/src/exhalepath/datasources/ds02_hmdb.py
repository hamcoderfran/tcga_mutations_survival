"""Priority 2 — HMDB + Human Breathomics-style VOC annotations."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

from .base import DataSource

# Seed HMDB IDs for core breath VOCs (publicly known)
VOC_HMDB = {
    "acetone": "HMDB0001659",
    "isoprene": "HMDB0000071",
    "ethanol": "HMDB0000108",
    "methanol": "HMDB0001875",
    "ammonia": "HMDB0000051",
    "acetaldehyde": "HMDB0000990",
    "hexanal": "HMDB0005994",
    "pentane": "HMDB0031645",
    "indole": "HMDB0000738",
    "phenol": "HMDB0000228",
    "trimethylamine": "HMDB0000906",
    "hydrogen_sulfide": "HMDB0000606",
    "limonene": "HMDB0000536",
    "toluene": "HMDB0003405",
    "benzene": "HMDB0001505",
    "formaldehyde": "HMDB0001426",
    "2_butanone": "HMDB0000474",
    "isopropanol": "HMDB0000863",
    "propanol": "HMDB0000820",
    "benzaldehyde": "HMDB0000615",
}


class HMDBSource(DataSource):
    priority = 2
    key = "hmdb"
    title = "HMDB + breath VOC annotation layer"
    description = "Metabolite identity, biofluid presence, and breath-relevant flags"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        rows = []
        for voc_id, hmdb_id in VOC_HMDB.items():
            row: dict[str, Any] = {
                "voc_id": voc_id,
                "hmdb_id": hmdb_id,
                "hmdb_url": f"https://hmdb.ca/metabolites/{hmdb_id}",
                "breath_reported": True,
                "source": "curated_hmdb_map",
            }
            if not offline:
                # HMDB XML endpoint (may be rate-limited)
                try:
                    r = requests.get(
                        f"https://hmdb.ca/metabolites/{hmdb_id}.xml",
                        timeout=30,
                        headers={"User-Agent": "ExhalePathAtlas/1.0"},
                    )
                    row["http_status"] = r.status_code
                    if r.ok:
                        text = r.text
                        row["in_blood"] = "Blood" in text or "blood" in text
                        row["in_urine"] = "Urine" in text or "urine" in text
                        row["in_saliva"] = "Saliva" in text or "saliva" in text
                        # crude name pull
                        if "<name>" in text:
                            row["hmdb_name"] = text.split("<name>", 1)[1].split("</name>", 1)[0]
                except Exception as e:  # noqa: BLE001
                    row["error"] = str(e)
            rows.append(row)

        # Breathomics DB style catalog (literature-reported breath compounds)
        hbdb = {
            "version": "1.0.0",
            "description": "Breath-reported VOC flags aligned to ExhalePath catalog (HBDB-style)",
            "n_vocs": len(rows),
            "vocs": rows,
        }
        path = self.write_json("hmdb_voc_annotations.json", hbdb)
        man = self.write_manifest(n_vocs=len(rows))
        return {"annotations": path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "hmdb_voc_annotations.json").read_text())
        out = knowledge_dir / "datasource_hmdb.json"
        out.write_text(json.dumps(doc, indent=2))
        # also enrich voc_catalog breath flags if present
        voc_path = knowledge_dir / "voc_catalog.json"
        if voc_path.exists():
            voc_doc = json.loads(voc_path.read_text())
            by_id = {v["voc_id"]: v for v in doc["vocs"]}
            for v in voc_doc.get("vocs") or []:
                ann = by_id.get(v["voc_id"])
                if not ann:
                    continue
                v["hmdb_id"] = ann.get("hmdb_id")
                v["breath_reported"] = True
                for k in ("in_blood", "in_urine", "in_saliva"):
                    if k in ann:
                        v[k] = ann[k]
            voc_path.write_text(json.dumps(voc_doc, indent=2))
        return {"path": str(out), "n_vocs": doc["n_vocs"]}
