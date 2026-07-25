"""Priority 2 — HMDB annotations via Wishart fix mirror (+ PubChem bridge fallback).

``hmdb.ca`` is Cloudflare-blocked in many automated environments. The Wishart lab
mirror ``hmdbfix.wishartlab.com`` serves the same metabolite XML / bulk zip.
"""

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

HMDB_MIRROR = "https://hmdbfix.wishartlab.com"
HMDB_MIRROR_ZIP = f"{HMDB_MIRROR}/system/downloads/current/hmdb_metabolites.zip"
HMDB_CANONICAL = "https://hmdb.ca"

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
            if len(hid) < 11:
                num = re.sub(r"\D", "", hid)
                hid = f"HMDB{int(num):07d}"
            return hid
    return None


def _fetch_metabolite_xml(hmdb_id: str, timeout: int = 30) -> tuple[int, str]:
    """Try Wishart mirror first, then canonical hmdb.ca."""
    urls = [
        f"{HMDB_MIRROR}/metabolites/{hmdb_id}.xml",
        f"{HMDB_CANONICAL}/metabolites/{hmdb_id}.xml",
    ]
    last_status = 0
    last_text = ""
    for url in urls:
        try:
            r = requests.get(url, timeout=timeout, headers=UA)
            last_status = r.status_code
            last_text = r.text if r.ok else ""
            if r.ok and "<metabolite>" in r.text:
                return r.status_code, r.text
        except Exception:  # noqa: BLE001
            continue
    return last_status, last_text


class HMDBSource(DataSource):
    priority = 2
    key = "hmdb"
    title = "HMDB + breath VOC annotation layer"
    description = (
        "Metabolite identity and biofluid presence via Wishart HMDB mirror "
        "(hmdbfix.wishartlab.com); PubChem bridge as fallback"
    )

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        panel_path = self.out_dir / "hmdb_panel_enrichment.json"
        breath_path = self.out_dir / "hmdb_breath_metabolites.json"

        # Prefer committed / pre-extracted Wishart bulk parse when present
        if panel_path.exists():
            panel_doc = json.loads(panel_path.read_text())
            by_voc = panel_doc.get("by_voc") or {}
            if by_voc:
                rows = []
                for voc_id in sorted(set(list(VOC_HMDB) + list(by_voc))):
                    rec = by_voc.get(voc_id) or {}
                    hmdb_id = rec.get("hmdb_id") or VOC_HMDB.get(voc_id)
                    rows.append(
                        {
                            "voc_id": voc_id,
                            "hmdb_id": hmdb_id,
                            "hmdb_url": (
                                f"{HMDB_CANONICAL}/metabolites/{hmdb_id}"
                                if hmdb_id
                                else None
                            ),
                            "hmdb_mirror_url": (
                                f"{HMDB_MIRROR}/metabolites/{hmdb_id}.xml"
                                if hmdb_id
                                else None
                            ),
                            "hmdb_name": rec.get("name"),
                            "pubchem_cid": rec.get("pubchem_cid")
                            or VOC_PUBCHEM_CID.get(voc_id),
                            "formula": rec.get("formula"),
                            "mw": rec.get("average_mass"),
                            "inchikey": rec.get("inchikey"),
                            "smiles": rec.get("smiles"),
                            "cas": rec.get("cas"),
                            "in_breath": bool(rec.get("in_breath")),
                            "in_blood": bool(rec.get("in_blood")),
                            "in_urine": bool(rec.get("in_urine")),
                            "in_saliva": bool(rec.get("in_saliva")),
                            "biospecimens": rec.get("biospecimens") or [],
                            "breath_reported": True,
                            "source": "hmdbfix_wishart_bulk_xml",
                        }
                    )
                breath_n = None
                if breath_path.exists():
                    breath_n = json.loads(breath_path.read_text()).get("n_breath")
                doc = {
                    "version": "3.0.0",
                    "description": (
                        "HMDB annotations for ExhalePath VOC panel from Wishart "
                        "fix mirror bulk XML (hmdb.ca Cloudflare-blocked)."
                    ),
                    "mirror": HMDB_MIRROR,
                    "mirror_zip": HMDB_MIRROR_ZIP,
                    "license": "CC-BY-NC-4.0",
                    "n_vocs": len(rows),
                    "n_breath_metabolites_in_hmdb": breath_n,
                    "vocs": rows,
                }
                path = self.write_json("hmdb_voc_annotations.json", doc)
                out = {
                    "annotations": path,
                    "manifest": self.write_manifest(
                        n_vocs=len(rows),
                        n_breath_metabolites=breath_n,
                        mirror=HMDB_MIRROR,
                        source="hmdbfix_wishart_bulk_xml",
                        offline=offline,
                    ),
                }
                if breath_path.exists():
                    out["breath_metabolites"] = breath_path
                if panel_path.exists():
                    out["panel_enrichment"] = panel_path
                return out

        # Online per-VOC XML via mirror (no bulk zip required)
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
                "hmdb_url": f"{HMDB_CANONICAL}/metabolites/{hmdb_id}" if hmdb_id else None,
                "hmdb_mirror_url": (
                    f"{HMDB_MIRROR}/metabolites/{hmdb_id}.xml" if hmdb_id else None
                ),
                "breath_reported": True,
                "pubchem_cid": VOC_PUBCHEM_CID.get(voc_id),
                "source": "curated_hmdb_map",
            }
            if offline:
                rows.append(row)
                continue

            if hmdb_id:
                status, text = _fetch_metabolite_xml(hmdb_id)
                row["http_status"] = status
                if status == 200 and text:
                    row["in_blood"] = ">Blood<" in text or ">blood<" in text
                    row["in_urine"] = ">Urine<" in text or ">urine<" in text
                    row["in_saliva"] = ">Saliva<" in text or ">saliva<" in text
                    row["in_breath"] = ">Breath<" in text or ">breath<" in text
                    if "<name>" in text:
                        row["hmdb_name"] = text.split("<name>", 1)[1].split("</name>", 1)[0]
                    m = re.search(r"<inchikey>([^<]+)</inchikey>", text)
                    if m:
                        row["inchikey"] = m.group(1)
                    m = re.search(r"<chemical_formula>([^<]+)</chemical_formula>", text)
                    if m:
                        row["formula"] = m.group(1)
                    row["source"] = "hmdbfix_wishart_xml"

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
                                row["hmdb_url"] = f"{HMDB_CANONICAL}/metabolites/{found}"
                            row["source"] = "pubchem_synonym_bridge"
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
                        row["smiles"] = p0.get("CanonicalSMILES") or p0.get(
                            "ConnectivitySMILES"
                        )
                    time.sleep(0.08)
                except Exception as e:  # noqa: BLE001
                    row["pubchem_error"] = str(e)
            rows.append(row)

        doc = {
            "version": "3.0.0",
            "description": "HMDB annotations for ExhalePath VOC panel (Wishart mirror + PubChem bridge)",
            "mirror": HMDB_MIRROR,
            "mirror_zip": HMDB_MIRROR_ZIP,
            "license": "CC-BY-NC-4.0",
            "n_vocs": len(rows),
            "vocs": rows,
        }
        path = self.write_json("hmdb_voc_annotations.json", doc)
        man = self.write_manifest(n_vocs=len(rows), mirror=HMDB_MIRROR)
        return {"annotations": path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "hmdb_voc_annotations.json").read_text())
        breath_path = self.out_dir / "hmdb_breath_metabolites.json"
        breath = json.loads(breath_path.read_text()) if breath_path.exists() else None
        fused = {
            **doc,
            "breath_metabolites_summary": (
                {
                    "n_breath": breath.get("n_breath"),
                    "n_metabolites_scanned": breath.get("n_metabolites_scanned"),
                    "source": breath.get("source"),
                    "path": str(breath_path),
                }
                if breath
                else None
            ),
        }
        out = knowledge_dir / "datasource_hmdb.json"
        out.write_text(json.dumps(fused, indent=2))
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
                for k in (
                    "in_blood",
                    "in_urine",
                    "in_saliva",
                    "in_breath",
                    "inchikey",
                    "formula",
                    "mw",
                    "smiles",
                    "cas",
                ):
                    if k in ann and ann[k] is not None:
                        v[k] = ann[k]
                n += 1
            voc_path.write_text(json.dumps(voc_doc, indent=2))
        return {
            "path": str(out),
            "n_vocs": doc["n_vocs"],
            "n_catalog_enriched": n,
            "n_breath_metabolites": (breath or {}).get("n_breath"),
        }
