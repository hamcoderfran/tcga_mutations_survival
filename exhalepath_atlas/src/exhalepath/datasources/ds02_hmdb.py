"""Priority 2 — HMDB annotations (PubChem-bridged when hmdb.ca is blocked)."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any

import requests

from ..ingest.public_breath import VOC_PUBCHEM_CID
from .base import DataSource

UA = {"User-Agent": "ExhalePathAtlas/1.0"}

# Seed HMDB IDs for core breath VOCs (publicly known) + expanded panel
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
    "heptanal": "HMDB0000517",
    "nonanal": "HMDB0005852",
    "decanal": "HMDB0011623",
    "dms": "HMDB0003237",
    "dimethyl_disulfide": "HMDB0005879",
    "ethane": "HMDB0003275",
    "butane": "HMDB0003223",
    "hexane": "HMDB0029599",
    "octane": "HMDB0031646",
    "styrene": "HMDB0003412",
    "ethylbenzene": "HMDB0003403",
    "xylene": "HMDB0003401",
    "cyclohexane": "HMDB0029598",
    "methyl_acetate": "HMDB0031647",
    "ethyl_acetate": "HMDB0031234",
    "acetonitrile": "HMDB0031296",
    "furan": "HMDB0031231",
    "2_pentanone": "HMDB0001865",
    "3_methylbutanal": "HMDB0006479",
    "octanal": "HMDB0005853",
    "undecane": "HMDB0031597",
    "dodecane": "HMDB0031598",
    "carbon_disulfide": "HMDB0003685",
    "methyl_mercaptan": "HMDB0003264",
    "dimethyl_amine": "HMDB0000087",
    "pyrrole": "HMDB0000305",
    "propionaldehyde": "HMDB0003366",
    "crotonaldehyde": "HMDB0031235",
    "dmts": "HMDB0031236",
    "allyl_methyl_sulfide": "HMDB0031237",
}


def _hmdb_from_synonyms(synonyms: list[str]) -> str | None:
    for s in synonyms:
        m = re.match(r"^(HMDB\d+)$", str(s).strip(), re.I)
        if m:
            hid = m.group(1).upper()
            # normalize to HMDB0000000 style when short
            if len(hid) < 11:
                num = re.sub(r"\D", "", hid)
                hid = f"HMDB{int(num):07d}"
            return hid
    return None


class HMDBSource(DataSource):
    priority = 2
    key = "hmdb"
    title = "HMDB + breath VOC annotation layer"
    description = "Metabolite identity, biofluid presence, and PubChem-bridged HMDB xrefs"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        prev = self.out_dir / "hmdb_voc_annotations.json"
        if offline and prev.exists():
            doc = json.loads(prev.read_text())
            path = self.write_json("hmdb_voc_annotations.json", doc)
            man = self.write_manifest(n_vocs=doc.get("n_vocs"), offline=True, reused=True)
            return {"annotations": path, "manifest": man}

        # Ensure every catalog VOC is considered
        catalog_path = self.root / "data" / "knowledge" / "voc_catalog.json"
        voc_ids = list(VOC_HMDB.keys())
        if catalog_path.exists():
            for v in json.loads(catalog_path.read_text()).get("vocs") or []:
                vid = v.get("voc_id")
                if vid and vid not in VOC_HMDB:
                    voc_ids.append(vid)
                elif vid and v.get("hmdb_id"):
                    VOC_HMDB.setdefault(vid, v["hmdb_id"])

        rows = []
        for voc_id in sorted(set(voc_ids)):
            hmdb_id = VOC_HMDB.get(voc_id)
            row: dict[str, Any] = {
                "voc_id": voc_id,
                "hmdb_id": hmdb_id,
                "hmdb_url": f"https://hmdb.ca/metabolites/{hmdb_id}" if hmdb_id else None,
                "breath_reported": True,
                "pubchem_cid": VOC_PUBCHEM_CID.get(voc_id),
                "source": "curated_hmdb_map",
            }
            if offline:
                rows.append(row)
                continue

            # Prefer live HMDB XML; fall back to PubChem synonyms for HMDB xref
            if hmdb_id:
                try:
                    r = requests.get(
                        f"https://hmdb.ca/metabolites/{hmdb_id}.xml",
                        timeout=30,
                        headers=UA,
                    )
                    row["http_status"] = r.status_code
                    if r.ok:
                        text = r.text
                        row["in_blood"] = "Blood" in text or "blood" in text
                        row["in_urine"] = "Urine" in text or "urine" in text
                        row["in_saliva"] = "Saliva" in text or "saliva" in text
                        if "<name>" in text:
                            row["hmdb_name"] = text.split("<name>", 1)[1].split("</name>", 1)[0]
                        row["source"] = "hmdb_xml"
                except Exception as e:  # noqa: BLE001
                    row["hmdb_xml_error"] = str(e)

            cid = row.get("pubchem_cid")
            if cid and (not row.get("hmdb_name") or row.get("http_status") not in (200,)):
                try:
                    r = requests.get(
                        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON",
                        timeout=45,
                        headers=UA,
                    )
                    if r.ok:
                        info = (r.json().get("InformationList") or {}).get("Information") or [{}]
                        syns = info[0].get("Synonym") or []
                        row["n_pubchem_synonyms"] = len(syns)
                        found = _hmdb_from_synonyms(syns)
                        if found:
                            row["hmdb_id_pubchem"] = found
                            if not row.get("hmdb_id"):
                                row["hmdb_id"] = found
                                row["hmdb_url"] = f"https://hmdb.ca/metabolites/{found}"
                            row["source"] = "pubchem_synonym_bridge"
                        # biofluid hints from synonym text
                        blob = " ".join(syns[:80]).lower()
                        row.setdefault("in_blood", "blood" in blob)
                        row.setdefault("in_urine", "urine" in blob)
                        row.setdefault("in_saliva", "saliva" in blob)
                    # properties for identity hardening
                    pr = requests.get(
                        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
                        "MolecularFormula,MolecularWeight,InChIKey,CanonicalSMILES/JSON",
                        timeout=45,
                        headers=UA,
                    )
                    if pr.ok:
                        props = (pr.json().get("PropertyTable") or {}).get("Properties") or [{}]
                        p0 = props[0]
                        row["formula"] = p0.get("MolecularFormula")
                        row["mw"] = p0.get("MolecularWeight")
                        row["inchikey"] = p0.get("InChIKey")
                        row["smiles"] = p0.get("CanonicalSMILES") or p0.get("ConnectivitySMILES")
                    time.sleep(0.08)
                except Exception as e:  # noqa: BLE001
                    row["pubchem_error"] = str(e)
            rows.append(row)

        doc = {
            "version": "2.0.0",
            "description": "HMDB annotations for ExhalePath VOC panel (XML + PubChem bridge)",
            "n_vocs": len(rows),
            "vocs": rows,
        }
        path = self.write_json("hmdb_voc_annotations.json", doc)
        man = self.write_manifest(n_vocs=len(rows))
        return {"annotations": path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "hmdb_voc_annotations.json").read_text())
        out = knowledge_dir / "datasource_hmdb.json"
        out.write_text(json.dumps(doc, indent=2))
        voc_path = knowledge_dir / "voc_catalog.json"
        n = 0
        if voc_path.exists():
            voc_doc = json.loads(voc_path.read_text())
            by_id = {v["voc_id"]: v for v in doc["vocs"]}
            for v in voc_doc.get("vocs") or []:
                ann = by_id.get(v["voc_id"])
                if not ann:
                    continue
                hid = ann.get("hmdb_id") or ann.get("hmdb_id_pubchem")
                if hid:
                    v["hmdb_id"] = hid
                v["breath_reported"] = True
                for k in ("in_blood", "in_urine", "in_saliva", "inchikey", "formula", "mw"):
                    if k in ann and ann[k] is not None:
                        v[k] = ann[k]
                n += 1
            voc_path.write_text(json.dumps(voc_doc, indent=2))
        return {"path": str(out), "n_vocs": doc["n_vocs"], "n_catalog_enriched": n}
