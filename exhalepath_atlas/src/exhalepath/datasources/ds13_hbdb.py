"""Priority 13 — Human Breathomics / EPA VOLATILOME full breath compound catalog.

HBDB (hbdb.cmdm.tw) is Cloudflare-blocked and the Zenodo SQL snapshot is access-
restricted. We integrate the open EPA CompTox VOLATILOME list (777 breath VOCs
from de Lacy Costello et al. + EPA measurements) — the same literature backbone
HBDB was built on — and fuse disease↔compound associations from atlas priors.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from ..ingest.secure_fetch import secure_fetch
from .base import DataSource

VOLATILOME_BLOB = (
    "https://clowder.edap-cluster.com/api/files/6165e731e4b0b85abf3ae0f5/blob"
)
VOLATILOME_DOI = "10.23645/epacomptox.9762410.v2"
HBDB_PAPER = "https://doi.org/10.1093/database/baz139"
DE_LACY = "https://doi.org/10.1088/1752-7155/8/1/014001"

# Literature-style disease → VOC name associations (HBDB disease view proxies)
HBDB_DISEASE_VOC_NAMES: dict[str, list[str]] = {
    "asthma": ["pentane", "ethane", "isoprene", "acetone", "nitric oxide", "hexanal"],
    "copd": ["pentane", "ethane", "hexanal", "heptanal", "acetone", "isoprene"],
    "lung_adenocarcinoma": [
        "hexanal",
        "heptanal",
        "nonanal",
        "2-butanone",
        "acetaldehyde",
        "isoprene",
    ],
    "type_2_diabetes": ["acetone", "isopropanol", "ethanol", "ethanol", "acetaldehyde"],
    "heart_failure": ["acetone", "pentane", "isoprene"],
    "cystic_fibrosis": ["hydrogen sulfide", "dimethyl sulfide", "2-pentanone", "acetone"],
    "liver_cirrhosis": ["limonene", "methanol", "2-butanone", "dimethyl sulfide"],
    "chronic_kidney_disease": ["trimethylamine", "ammonia", "dimethylamine", "isoprene"],
    "ulcerative_colitis": ["indole", "phenol", "hydrogen sulfide", "dimethyl disulfide"],
    "parkinson_disease": ["pentane", "hexane", "styrene", "ethylbenzene"],
    "alzheimer_disease": ["ethane", "pentane", "acetone"],
    "malaria": ["terpenes", "isoprene", "acetone"],
    "pneumonia_bacterial": ["hexanal", "heptanal", "acetone", "isoprene"],
    "epilepsy": ["acetone", "isoprene", "dimethyl sulfide"],
}


def _norm(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


class HBDBSource(DataSource):
    priority = 13
    key = "hbdb"
    title = "HBDB / EPA VOLATILOME breath compound catalog"
    description = (
        "Full open breath VOC compound list (VOLATILOME, 777 chemicals) + "
        "HBDB-style disease associations fused into atlas knowledge"
    )

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        compounds_path = self.out_dir / "volatilome_compounds.json"
        compounds: list[dict[str, Any]] = []

        if offline and compounds_path.exists():
            compounds = json.loads(compounds_path.read_text()).get("compounds") or []
        elif not offline:
            xls_path = self.out_dir / "Volatilome_subset_Sept_2_2019.xls"
            # Allowlisted quarantine→promote fetch; never execute blob contents.
            secure_fetch(
                VOLATILOME_BLOB,
                dest=xls_path,
                quarantine_dir=self.out_dir / ".quarantine",
                max_bytes=30 * 1024 * 1024,
                timeout=180,
            )
            df = pd.read_excel(xls_path, engine="xlrd")
            for _, row in df.iterrows():
                compounds.append(
                    {
                        "dtxsid": row.get("DTXSID"),
                        "name": row.get("PREFERRED_NAME"),
                        "casrn": row.get("CASRN"),
                        "inchikey": row.get("INCHIKEY"),
                        "iupac_name": row.get("IUPAC_NAME"),
                        "smiles": row.get("SMILES"),
                        "formula": row.get("MOLECULAR_FORMULA"),
                        "average_mass": row.get("AVERAGE_MASS"),
                        "source": "EPA_VOLATILOME",
                    }
                )
        else:
            # Minimal offline seed from atlas catalog CAS/names
            voc_path = self.root / "data" / "knowledge" / "voc_catalog.json"
            if voc_path.exists():
                for v in json.loads(voc_path.read_text()).get("vocs") or []:
                    compounds.append(
                        {
                            "name": v.get("name"),
                            "casrn": v.get("cas"),
                            "formula": v.get("formula"),
                            "atlas_voc_id": v.get("voc_id"),
                            "source": "atlas_offline_seed",
                        }
                    )

        # Map to atlas VOC panel via CAS / normalized name
        voc_path = self.root / "data" / "knowledge" / "voc_catalog.json"
        by_cas: dict[str, str] = {}
        by_name: dict[str, str] = {}
        if voc_path.exists():
            for v in json.loads(voc_path.read_text()).get("vocs") or []:
                if v.get("cas"):
                    by_cas[_norm(v["cas"])] = v["voc_id"]
                by_name[_norm(v.get("name") or "")] = v["voc_id"]
                by_name[_norm(v.get("voc_id") or "").replace("_", " ")] = v["voc_id"]

        n_mapped = 0
        for c in compounds:
            cas = _norm(c.get("casrn") or "")
            name = _norm(c.get("name") or "")
            vid = by_cas.get(cas) or by_name.get(name)
            if vid:
                c["atlas_voc_id"] = vid
                n_mapped += 1

        compound_doc = {
            "version": "1.0.0",
            "description": (
                "EPA CompTox VOLATILOME breath chemicals (de Lacy Costello backbone "
                "used by HBDB). HBDB website/SQL not redistributable here."
            ),
            "references": {
                "volatilome_doi": VOLATILOME_DOI,
                "de_lacy_costello": DE_LACY,
                "hbdb_paper": HBDB_PAPER,
                "hbdb_site": "https://hbdb.cmdm.tw/",
                "hbdb_access_note": (
                    "Live HBDB HTML blocked (Cloudflare); Zenodo hbdb2.sql restricted. "
                    "VOLATILOME is the open compound-layer substitute."
                ),
            },
            "n_compounds": len(compounds),
            "n_mapped_to_atlas": n_mapped,
            "compounds": compounds,
        }
        c_path = self.write_json("volatilome_compounds.json", compound_doc)

        # Disease associations: curated HBDB-style names + atlas voc_log2fc_prior
        disease_assoc: list[dict[str, Any]] = []
        priors_path = self.root / "data" / "knowledge" / "disease_voc_priors.json"
        if priors_path.exists():
            for d in json.loads(priors_path.read_text()).get("diseases") or []:
                did = d.get("disease_id")
                voc_prior = d.get("voc_log2fc_prior") or {}
                names = list(HBDB_DISEASE_VOC_NAMES.get(did or "", []))
                disease_assoc.append(
                    {
                        "disease_id": did,
                        "disease_name": d.get("name"),
                        "hbdb_style_voc_names": names,
                        "atlas_voc_ids": sorted(voc_prior.keys()),
                        "n_atlas_voc_priors": len(voc_prior),
                        "source": "atlas_priors+hbdb_literature_proxy",
                    }
                )
        # Ensure curated-only diseases also appear
        have = {a["disease_id"] for a in disease_assoc}
        for did, names in HBDB_DISEASE_VOC_NAMES.items():
            if did in have:
                continue
            disease_assoc.append(
                {
                    "disease_id": did,
                    "hbdb_style_voc_names": names,
                    "atlas_voc_ids": [],
                    "source": "hbdb_literature_proxy",
                }
            )

        d_path = self.write_json(
            "hbdb_disease_associations.json",
            {
                "version": "1.0.0",
                "n_diseases": len(disease_assoc),
                "diseases": disease_assoc,
            },
        )
        man = self.write_manifest(
            n_compounds=len(compounds),
            n_mapped_to_atlas=n_mapped,
            n_diseases=len(disease_assoc),
        )
        return {"compounds": c_path, "disease_associations": d_path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        compounds = json.loads((self.out_dir / "volatilome_compounds.json").read_text())
        diseases = json.loads((self.out_dir / "hbdb_disease_associations.json").read_text())
        # Slim knowledge fragment (full compound list kept under datasources/hbdb)
        mapped = [c for c in compounds.get("compounds") or [] if c.get("atlas_voc_id")]
        fused = {
            "version": "1.0.0",
            "n_compounds": compounds.get("n_compounds"),
            "n_mapped_to_atlas": compounds.get("n_mapped_to_atlas"),
            "references": compounds.get("references"),
            "atlas_mapped_compounds": mapped,
            "disease_associations": diseases.get("diseases") or [],
            "full_compound_table": str(self.out_dir / "volatilome_compounds.json"),
        }
        out = knowledge_dir / "datasource_hbdb.json"
        out.write_text(json.dumps(fused, indent=2))

        voc_path = knowledge_dir / "voc_catalog.json"
        n = 0
        if voc_path.exists():
            by_voc = {c["atlas_voc_id"]: c for c in mapped}
            voc_doc = json.loads(voc_path.read_text())
            for v in voc_doc.get("vocs") or []:
                hit = by_voc.get(v["voc_id"])
                if not hit:
                    continue
                v["breath_volatilome"] = True
                v["dtxsid"] = hit.get("dtxsid")
                if hit.get("inchikey"):
                    v["inchikey"] = hit["inchikey"]
                if hit.get("smiles"):
                    v["smiles"] = hit["smiles"]
                n += 1
            voc_path.write_text(json.dumps(voc_doc, indent=2))
        return {
            "path": str(out),
            "n_compounds": fused["n_compounds"],
            "n_atlas_enriched": n,
            "n_diseases": len(fused["disease_associations"]),
        }
