"""Priority 13 — HBDB disease↔VOC associations + EPA VOLATILOME catalog.

Live HBDB (https://hbdb.cmdm.tw/?disease_page=1&tab=disease) is Cloudflare-blocked
in this environment. The open Zenodo dump ``hbdb2_wo_sentences.sql``
(10.5281/zenodo.14958797) provides all 60 diseases and compound links; committed
JSON under ``data/datasources/hbdb/`` is the offline source of truth.

EPA CompTox VOLATILOME remains the open full breath-compound catalog layer.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from ..ingest.secure_fetch import secure_fetch
from .base import DataSource
from .hbdb_mapping import build_atlas_mapped_doc

VOLATILOME_BLOB = (
    "https://clowder.edap-cluster.com/api/files/6165e731e4b0b85abf3ae0f5/blob"
)
VOLATILOME_DOI = "10.23645/epacomptox.9762410.v2"
HBDB_PAPER = "https://doi.org/10.1093/database/baz139"
DE_LACY = "https://doi.org/10.1088/1752-7155/8/1/014001"
ZENODO_SQL = (
    "https://zenodo.org/records/14958797/files/hbdb2_wo_sentences.sql?download=1"
)
ZENODO_DOI = "10.5281/zenodo.14958797"

# Legacy literature-style proxies (used only when Zenodo extract is absent)
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
    "type_2_diabetes": ["acetone", "isopropanol", "ethanol", "acetaldehyde"],
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


def load_hbdb_atlas_associations(root: Path) -> list[dict[str, Any]]:
    """Return atlas-mapped HBDB associations (Zenodo extract preferred)."""
    mapped = root / "data" / "datasources" / "hbdb" / "hbdb_atlas_mapped.json"
    if mapped.exists():
        doc = json.loads(mapped.read_text())
        return list(doc.get("atlas_associations") or [])
    return []


class HBDBSource(DataSource):
    priority = 13
    key = "hbdb"
    title = "HBDB disease–VOC associations + EPA VOLATILOME catalog"
    description = (
        "All 60 HBDB diseases with compound links (Zenodo SQL) plus the open "
        "EPA VOLATILOME breath VOC catalog (777 chemicals)"
    )

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        compounds_path = self.out_dir / "volatilome_compounds.json"
        compounds: list[dict[str, Any]] = []

        if offline and compounds_path.exists():
            compounds = json.loads(compounds_path.read_text()).get("compounds") or []
        elif not offline:
            xls_path = self.out_dir / "Volatilome_subset_Sept_2_2019.xls"
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

        voc_path = self.root / "data" / "knowledge" / "voc_catalog.json"
        by_cas: dict[str, str] = {}
        by_name: dict[str, str] = {}
        atlas_vocs: list[dict[str, Any]] = []
        if voc_path.exists():
            atlas_vocs = list(json.loads(voc_path.read_text()).get("vocs") or [])
            for v in atlas_vocs:
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
            "version": "1.1.0",
            "description": (
                "EPA CompTox VOLATILOME breath chemicals (de Lacy Costello backbone). "
                "HBDB disease associations are in hbdb_disease_vocs_60.json."
            ),
            "references": {
                "volatilome_doi": VOLATILOME_DOI,
                "de_lacy_costello": DE_LACY,
                "hbdb_paper": HBDB_PAPER,
                "hbdb_site": "https://hbdb.cmdm.tw/?disease_page=1&tab=disease",
                "zenodo_doi": ZENODO_DOI,
                "zenodo_sql": ZENODO_SQL,
                "hbdb_access_note": (
                    "Live HBDB HTML blocked (Cloudflare). Open Zenodo dump "
                    "hbdb2_wo_sentences.sql provides all 60 disease↔compound links."
                ),
            },
            "n_compounds": len(compounds),
            "n_mapped_to_atlas": n_mapped,
            "compounds": compounds,
        }
        c_path = self.write_json("volatilome_compounds.json", compound_doc)

        # Prefer committed Zenodo extract (all 60 diseases). Optional online refresh
        # when extract is missing: fetch SQL + run scripts/extract_hbdb_zenodo.py.
        disease_vocs_path = self.out_dir / "hbdb_disease_vocs_60.json"
        compounds_hbdb_path = self.out_dir / "hbdb_compounds.json"
        if not offline and not disease_vocs_path.exists():
            self._try_refresh_from_zenodo(disease_vocs_path, compounds_hbdb_path)

        disease_assoc: list[dict[str, Any]] = []
        mapped_doc: dict[str, Any] | None = None
        if disease_vocs_path.exists() and compounds_hbdb_path.exists() and atlas_vocs:
            disease_vocs = json.loads(disease_vocs_path.read_text())
            compounds_hbdb = json.loads(compounds_hbdb_path.read_text())
            mapped_doc = build_atlas_mapped_doc(disease_vocs, compounds_hbdb, atlas_vocs)
            self.write_json("hbdb_atlas_mapped.json", mapped_doc)
            for a in mapped_doc.get("atlas_associations") or []:
                disease_assoc.append(
                    {
                        "disease_id": a["disease_id"],
                        "hbdb_style_voc_names": list(a.get("hbdb_style_voc_names") or []),
                        "atlas_voc_ids": sorted((a.get("vocs") or {}).keys()),
                        "n_atlas_voc_priors": len(a.get("vocs") or {}),
                        "n_hbdb_panel_vocs": a.get("n_atlas_vocs"),
                        "source": "hbdb_zenodo_sql",
                        "voc_evidence": a.get("vocs") or {},
                    }
                )
            # Inventory rows for HBDB diseases without panel overlap
            have = {a["disease_id"] for a in disease_assoc}
            for d in mapped_doc.get("diseases") or []:
                did = d.get("atlas_disease_id")
                if not did or did in have:
                    continue
                disease_assoc.append(
                    {
                        "disease_id": did,
                        "hbdb_disease_id": d.get("hbdb_disease_id"),
                        "hbdb_disease_name": d.get("name"),
                        "hbdb_style_voc_names": [],
                        "atlas_voc_ids": [],
                        "n_hbdb_compounds": d.get("n_compounds"),
                        "source": "hbdb_zenodo_sql_inventory",
                    }
                )

        if not disease_assoc:
            # Fallback: legacy proxy + atlas priors
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
                "version": "2.1.0",
                "n_diseases": len(disease_assoc),
                "source": (
                    "hbdb_zenodo_sql"
                    if mapped_doc
                    else "hbdb_literature_proxy"
                ),
                "n_hbdb_diseases_raw": (mapped_doc or {}).get("n_hbdb_diseases"),
                "n_disease_compound_links": (mapped_doc or {}).get(
                    "n_disease_compound_links"
                ),
                "diseases": disease_assoc,
            },
        )
        man = self.write_manifest(
            n_compounds=len(compounds),
            n_mapped_to_atlas=n_mapped,
            n_diseases=len(disease_assoc),
            n_hbdb_diseases=(mapped_doc or {}).get("n_hbdb_diseases"),
            n_hbdb_links=(mapped_doc or {}).get("n_disease_compound_links"),
            zenodo_doi=ZENODO_DOI,
        )
        out = {"compounds": c_path, "disease_associations": d_path, "manifest": man}
        if disease_vocs_path.exists():
            out["disease_vocs_60"] = disease_vocs_path
        if (self.out_dir / "hbdb_atlas_mapped.json").exists():
            out["atlas_mapped"] = self.out_dir / "hbdb_atlas_mapped.json"
        return out

    def _try_refresh_from_zenodo(
        self, disease_vocs_path: Path, compounds_hbdb_path: Path
    ) -> None:
        import importlib.util

        sql_path = self.out_dir / ".quarantine" / "hbdb2_wo_sentences.sql"
        try:
            secure_fetch(
                ZENODO_SQL,
                dest=sql_path,
                quarantine_dir=self.out_dir / ".quarantine",
                max_bytes=400 * 1024 * 1024,
                timeout=300,
            )
        except Exception:
            return
        script = self.root / "scripts" / "extract_hbdb_zenodo.py"
        if not script.exists() or not sql_path.exists():
            return
        try:
            spec = importlib.util.spec_from_file_location("extract_hbdb_zenodo", script)
            if spec is None or spec.loader is None:
                return
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            extracted = mod.extract_from_sql(sql_path)
            disease_vocs_path.write_text(json.dumps(extracted["disease_vocs"], indent=2))
            compounds_hbdb_path.write_text(json.dumps(extracted["compounds"], indent=2))
        except Exception:
            return

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        compounds = json.loads((self.out_dir / "volatilome_compounds.json").read_text())
        diseases = json.loads(
            (self.out_dir / "hbdb_disease_associations.json").read_text()
        )
        mapped_path = self.out_dir / "hbdb_atlas_mapped.json"
        mapped = json.loads(mapped_path.read_text()) if mapped_path.exists() else {}
        raw_path = self.out_dir / "hbdb_disease_vocs_60.json"
        raw_meta = {}
        if raw_path.exists():
            raw = json.loads(raw_path.read_text())
            raw_meta = {
                "n_hbdb_diseases": raw.get("n_diseases"),
                "n_disease_compound_links": raw.get("n_disease_compound_links"),
                "n_unique_compounds_linked": raw.get("n_unique_compounds_linked"),
                "diseases_with_zero_compounds": raw.get("diseases_with_zero_compounds"),
                "zenodo_doi": raw.get("zenodo_doi"),
            }

        mapped_compounds = [
            c for c in compounds.get("compounds") or [] if c.get("atlas_voc_id")
        ]
        fused = {
            "version": "2.1.0",
            "n_compounds": compounds.get("n_compounds"),
            "n_mapped_to_atlas": compounds.get("n_mapped_to_atlas"),
            "references": compounds.get("references"),
            "hbdb_extract": raw_meta,
            "atlas_mapped_summary": {
                "n_hbdb_diseases": mapped.get("n_hbdb_diseases"),
                "n_hbdb_diseases_mapped_to_atlas": mapped.get(
                    "n_hbdb_diseases_mapped_to_atlas"
                ),
                "n_disease_compound_links": mapped.get("n_disease_compound_links"),
                "n_links_mapped_to_atlas_panel": mapped.get(
                    "n_links_mapped_to_atlas_panel"
                ),
                "n_atlas_diseases_with_panel_vocs": mapped.get(
                    "n_atlas_diseases_with_panel_vocs"
                ),
            },
            "atlas_mapped_compounds": mapped_compounds,
            "disease_associations": diseases.get("diseases") or [],
            "full_compound_table": str(self.out_dir / "volatilome_compounds.json"),
            "hbdb_disease_vocs_60": str(raw_path) if raw_path.exists() else None,
            "hbdb_atlas_mapped": str(mapped_path) if mapped_path.exists() else None,
        }
        out = knowledge_dir / "datasource_hbdb.json"
        out.write_text(json.dumps(fused, indent=2))

        voc_path = knowledge_dir / "voc_catalog.json"
        n = 0
        if voc_path.exists():
            by_voc = {c["atlas_voc_id"]: c for c in mapped_compounds}
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
            "n_hbdb_links": raw_meta.get("n_disease_compound_links"),
        }
