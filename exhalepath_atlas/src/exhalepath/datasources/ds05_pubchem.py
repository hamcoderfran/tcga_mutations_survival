"""Priority 5 — PubChem / CompTox-style physicochemical properties."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import requests

from ..ingest.public_breath import VOC_PUBCHEM_CID
from .base import DataSource


class PubChemSource(DataSource):
    priority = 5
    key = "pubchem"
    title = "PubChem physicochemical properties"
    description = "MW, XLogP, TPSA, H-bond counts for VOC catalog via PUG REST"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        rows = []
        for voc_id, cid in VOC_PUBCHEM_CID.items():
            row: dict[str, Any] = {"voc_id": voc_id, "cid": cid}
            if offline:
                row.update({"xlogp": None, "mw": None, "tpsa": None, "offline": True})
                rows.append(row)
                continue
            try:
                url = (
                    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
                    "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,Complexity/JSON"
                )
                r = requests.get(url, timeout=45, headers={"User-Agent": "ExhalePathAtlas/1.0"})
                if r.ok:
                    props = (r.json().get("PropertyTable") or {}).get("Properties") or [{}]
                    p0 = props[0]
                    row.update(
                        {
                            "mw": p0.get("MolecularWeight"),
                            "xlogp": p0.get("XLogP"),
                            "tpsa": p0.get("TPSA"),
                            "hbd": p0.get("HBondDonorCount"),
                            "hba": p0.get("HBondAcceptorCount"),
                            "complexity": p0.get("Complexity"),
                        }
                    )
                else:
                    row["http_status"] = r.status_code
                time.sleep(0.08)
            except Exception as e:  # noqa: BLE001
                row["error"] = str(e)
            rows.append(row)
        doc = {"version": "1.0.0", "n_vocs": len(rows), "properties": rows}
        path = self.write_json("pubchem_voc_properties.json", doc)
        return {"properties": path, "manifest": self.write_manifest(n_vocs=len(rows))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "pubchem_voc_properties.json").read_text())
        out = knowledge_dir / "datasource_pubchem.json"
        out.write_text(json.dumps(doc, indent=2))
        # merge into chembl_pathway_priors voc_physchem-compatible block
        chembl_path = knowledge_dir / "chembl_pathway_priors.json"
        chembl = json.loads(chembl_path.read_text()) if chembl_path.exists() else {"pathways": {}, "voc_physchem": {}}
        vp = chembl.setdefault("voc_physchem", {})
        n = 0
        for row in doc["properties"]:
            if row.get("mw") is None and row.get("xlogp") is None:
                continue
            cur = vp.setdefault(row["voc_id"], {})
            if row.get("mw") is not None:
                cur["full_mwt"] = float(row["mw"])
            if row.get("xlogp") is not None:
                cur["alogp"] = float(row["xlogp"])
            if row.get("tpsa") is not None:
                cur["psa"] = float(row["tpsa"])
            cur["pubchem_cid"] = row.get("cid")
            n += 1
        chembl_path.write_text(json.dumps(chembl, indent=2))
        return {"path": str(out), "n_vocs_enriched": n}
