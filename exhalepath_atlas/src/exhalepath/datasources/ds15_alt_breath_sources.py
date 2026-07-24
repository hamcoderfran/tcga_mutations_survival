"""Priority 15 — Audit + harvest open alternatives to Owlstone VOC Atlas / live HBDB.

Covers the publicly discussed breath-data options:
- Clinical Breathomics (Figshare / Scientific Data) — already in public_breath/
- Zenodo open breathomics adjuncts (common-78 features, Tedlar variability, HBDB eval)
- HMDB bulk (Cloudflare-blocked; seed annotations remain in ds02)
- PhysioNet respiratory waveform catalogs (modality ≠ VOC chemistry)
- Breathomix BreathBase / Shirley Human Breath Atlas / Owlstone VOC Atlas (gated)

Does not scrape Cloudflare-blocked sites or violate Owlstone license terms.
"""

from __future__ import annotations

import csv
import io
import json
import re
import zipfile
from pathlib import Path
from typing import Any

import requests

from .base import DataSource

UA = {"User-Agent": "ExhalePathAtlas/1.0 (research; open-data harvest)"}

ZENODO_FILES = {
    "common78_features.csv": (
        "https://zenodo.org/api/records/21238613/files/common78_features.csv/content",
        "10.5281/zenodo.21238613",
    ),
    "table3_metrics.csv": (
        "https://zenodo.org/api/records/21238613/files/table3_metrics.csv/content",
        "10.5281/zenodo.21238613",
    ),
    "Tedlar_temp_raw.xlsx": (
        "https://zenodo.org/api/records/21135284/files/Tedlar(temp)_raw.xlsx/content",
        "10.5281/zenodo.21135284",
    ),
    "hbdb_eval_dataset.zip": (
        "https://zenodo.org/api/records/14958797/files/eval_dataset.zip/content",
        "10.5281/zenodo.14958797",
    ),
}

FIGSHARE_CBD = {
    "article_doi": "10.6084/m9.figshare.23522490",
    "paper_doi": "10.1038/s41597-024-03052-2",
    "title": "Clinical Breathomics Dataset (asthma / COPD / bronchiectasis GC-MS)",
}


def _get(url: str, *, timeout: int = 120) -> bytes:
    r = requests.get(url, headers=UA, timeout=timeout)
    r.raise_for_status()
    return r.content


def _physionet_respir_catalog() -> list[dict[str, str]]:
    html = _get("https://physionet.org/about/database/", timeout=60).decode(
        "utf-8", "replace"
    )
    out: list[dict[str, str]] = []
    seen: set[str] = set()
    for m in re.finditer(
        r'href="(/content/([a-z0-9-]+)/?)"[^>]*>([^<]{3,160})', html, re.I
    ):
        slug, title = m.group(2), re.sub(r"\s+", " ", m.group(3)).strip()
        blob = f"{slug} {title}".lower()
        if not re.search(r"breath|respir|ventil|cpap|lung|airway|spirom", blob):
            continue
        if slug in seen:
            continue
        seen.add(slug)
        out.append(
            {
                "slug": slug,
                "title": title,
                "url": f"https://physionet.org/content/{slug}/",
            }
        )
    return out


def _index_eval_zip(zip_bytes: bytes) -> dict[str, Any]:
    compounds: dict[str, int] = {}
    n_json = 0
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        for name in zf.namelist():
            if "MACOSX" in name or name.endswith("/"):
                continue
            if not name.endswith(".json"):
                continue
            n_json += 1
            # eval_dataset/evaluation/<compound>/...
            parts = Path(name).parts
            if len(parts) >= 3 and parts[0] == "eval_dataset" and parts[1] == "evaluation":
                compounds[parts[2]] = compounds.get(parts[2], 0) + 1
    return {
        "n_json_records": n_json,
        "compounds": [
            {"name": k, "n_records": v} for k, v in sorted(compounds.items())
        ],
    }


def _parse_common78(csv_bytes: bytes) -> list[dict[str, Any]]:
    text = csv_bytes.decode("utf-8", "replace")
    rows = list(csv.DictReader(io.StringIO(text)))
    out = []
    for r in rows:
        cid = r.get("pubchem_CID") or r.get("pubchem_cid")
        name = r.get("IUPAC Name") or r.get("IUPAC name") or r.get("name")
        try:
            cid_i = int(float(cid)) if cid not in (None, "") else None
        except ValueError:
            cid_i = None
        out.append({"pubchem_cid": cid_i, "iupac_name": (name or "").strip()})
    return out


def build_access_matrix(
    *,
    cbd_present: bool,
    hmdb_annotations_present: bool,
    zenodo_files: list[str],
    physionet_n: int,
) -> dict[str, Any]:
    """Status of each source from the common 'human breath database' shortlist."""
    return {
        "version": "1.0.0",
        "description": (
            "Access audit for bulk VOC / breath databases suggested as alternatives "
            "to Owlstone VOC Atlas en-masse download."
        ),
        "sources": [
            {
                "key": "hbdb",
                "name": "Human Breathomics Database",
                "url": "https://hbdb.cmdm.tw/?disease_page=1&tab=disease",
                "modality": "VOC / disease associations",
                "status": "partial_open",
                "access": (
                    "Live site Cloudflare-blocked for automated scrape. "
                    "Open Zenodo SQL dump hbdb2_wo_sentences.sql "
                    "(10.5281/zenodo.14958797) is the bulk path; see HBDB scrape PR."
                ),
                "bulk_download": "zenodo_sql",
                "in_atlas": True,
                "notes": "Google-style BeautifulSoup scrape of hbdb.cmdm.tw fails with 403 here.",
            },
            {
                "key": "clinical_breathomics_figshare",
                "name": "Clinical Breathomics Dataset (Scientific Data / Figshare)",
                "url": "https://doi.org/10.6084/m9.figshare.23522490",
                "paper": FIGSHARE_CBD["paper_doi"],
                "modality": "GC-MS peak tables + clinical metadata",
                "status": "open_downloaded",
                "access": "Anonymous Figshare file download (peak CSVs + metadata xlsx).",
                "bulk_download": "figshare",
                "in_atlas": cbd_present,
                "local_dir": "data/public_breath/",
                "harvest_cmd": "voc harvest-public-breath",
            },
            {
                "key": "owlstone_voc_atlas",
                "name": "Breath Biopsy VOC Atlas (Owlstone)",
                "url": "https://www.owlstonemedical.com/",
                "modality": "VOC catalog / OMNI feature tables",
                "status": "gated_license_restricted",
                "access": (
                    "Registration required. License bars AI training / fine-tuning / "
                    "validation — do not ingest into ExhalePath ML calibrator."
                ),
                "bulk_download": "manual_registration_only",
                "in_atlas": False,
            },
            {
                "key": "physionet",
                "name": "PhysioNet respiratory / vital-sign databases",
                "url": "https://physionet.org/about/database/",
                "modality": "waveforms / ventilation (not VOC chemistry)",
                "status": "open_cataloged",
                "access": "Public database listing; waveform dumps are large and out of VOC scope.",
                "bulk_download": "physionet_per_project",
                "in_atlas": False,
                "n_respir_related_cataloged": physionet_n,
                "notes": "Useful for breathing-pattern models, not breath VOC panels.",
            },
            {
                "key": "breathomix_breathbase",
                "name": "BreathBase® Data (Breathomix)",
                "url": "https://www.breathomix.com/breathbase-data/",
                "modality": "eNose breath profiles + clinical data",
                "status": "gated_contact",
                "access": "No public CSV/ZIP; access via vendor contact / partnership.",
                "bulk_download": None,
                "in_atlas": False,
            },
            {
                "key": "hmdb",
                "name": "Human Metabolome Database",
                "url": "https://hmdb.ca/downloads",
                "modality": "metabolite structures / biofluids (incl. breath subset)",
                "status": "cloudflare_blocked_bulk",
                "access": (
                    "hmdb.ca downloads return Cloudflare 403 here. Atlas keeps "
                    "PubChem-bridged VOC↔HMDB seed annotations (ds02)."
                ),
                "bulk_download": "blocked_here",
                "in_atlas": hmdb_annotations_present,
                "local_artifact": "data/datasources/hmdb/hmdb_voc_annotations.json",
            },
            {
                "key": "shirley_human_breath_atlas",
                "name": "Human Breath Atlas (Shirley Diagnostics)",
                "url": "https://www.shirleydiagnostics.com/",
                "modality": "end-tidal VOC reference (vendor)",
                "status": "gated_marketing",
                "access": "No public bulk export found on the marketing site.",
                "bulk_download": None,
                "in_atlas": False,
            },
            {
                "key": "zenodo_breathomics_adjuncts",
                "name": "Zenodo open breathomics adjunct datasets",
                "url": "https://zenodo.org/",
                "modality": "feature lists / sampling variability / HBDB eval texts",
                "status": "open_downloaded",
                "access": "Anonymous Zenodo file API.",
                "bulk_download": "zenodo",
                "in_atlas": True,
                "files": zenodo_files,
            },
        ],
    }


class AltBreathSources(DataSource):
    priority = 15
    key = "alt_breath_sources"
    title = "Open alternatives to VOC Atlas / live HBDB (access audit + Zenodo harvest)"
    description = (
        "Catalog gated vs open breath VOC sources; download Zenodo adjuncts; "
        "confirm Figshare clinical breathomics presence"
    )

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        raw_dir = self.out_dir / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)

        zenodo_meta: list[dict[str, Any]] = []
        common78: list[dict[str, Any]] = []
        eval_index: dict[str, Any] = {}

        if offline:
            # Reuse committed artifacts
            for name in ZENODO_FILES:
                p = raw_dir / name
                if p.exists():
                    zenodo_meta.append({"file": name, "path": str(p), "offline": True})
            c78 = self.out_dir / "common78_features.json"
            if c78.exists():
                common78 = json.loads(c78.read_text()).get("features") or []
            ev = self.out_dir / "hbdb_eval_index.json"
            if ev.exists():
                eval_index = json.loads(ev.read_text())
            physio = json.loads(
                (self.out_dir / "physionet_respir_catalog.json").read_text()
            ) if (self.out_dir / "physionet_respir_catalog.json").exists() else {"databases": []}
        else:
            for name, (url, doi) in ZENODO_FILES.items():
                dest = raw_dir / name
                try:
                    blob = _get(url, timeout=180)
                    dest.write_bytes(blob)
                    zenodo_meta.append(
                        {
                            "file": name,
                            "doi": doi,
                            "url": url,
                            "n_bytes": len(blob),
                            "path": str(dest),
                        }
                    )
                    if name == "common78_features.csv":
                        common78 = _parse_common78(blob)
                    if name == "hbdb_eval_dataset.zip":
                        eval_index = _index_eval_zip(blob)
                        eval_index["doi"] = doi
                except Exception as e:  # noqa: BLE001 — record per-file failures
                    zenodo_meta.append(
                        {"file": name, "doi": doi, "url": url, "error": str(e)}
                    )
            try:
                physio_list = _physionet_respir_catalog()
            except Exception as e:  # noqa: BLE001
                physio_list = []
                physio = {"databases": [], "error": str(e)}
            else:
                physio = {"n_databases": len(physio_list), "databases": physio_list}

        cbd_dir = self.root / "data" / "public_breath"
        cbd_present = (cbd_dir / "Asthma_peaktable_ver3.csv").exists()
        hmdb_present = (
            self.root / "data" / "datasources" / "hmdb" / "hmdb_voc_annotations.json"
        ).exists()

        matrix = build_access_matrix(
            cbd_present=cbd_present,
            hmdb_annotations_present=hmdb_present,
            zenodo_files=[m.get("file") for m in zenodo_meta if not m.get("error")],
            physionet_n=int(physio.get("n_databases") or len(physio.get("databases") or [])),
        )
        matrix_path = self.write_json("ACCESS_MATRIX.json", matrix)
        physio_path = self.write_json("physionet_respir_catalog.json", physio)
        c78_path = self.write_json(
            "common78_features.json",
            {
                "version": "1.0.0",
                "doi": "10.5281/zenodo.21238613",
                "description": (
                    "Common-78 breathomics feature list (PubChem CID + IUPAC) from "
                    "Zenodo aggregation-stability study."
                ),
                "n_features": len(common78),
                "features": common78,
            },
        )
        eval_path = self.write_json(
            "hbdb_eval_index.json",
            {
                "version": "1.0.0",
                "doi": "10.5281/zenodo.14958797",
                "description": (
                    "Index of HBDB manually curated evaluation JSON snippets "
                    "(full zip kept under raw/ when harvested online)."
                ),
                **eval_index,
            },
        )
        zenodo_path = self.write_json(
            "zenodo_harvest.json",
            {"version": "1.0.0", "files": zenodo_meta},
        )

        open_n = sum(1 for s in matrix["sources"] if s["status"].startswith("open"))
        gated_n = sum(
            1
            for s in matrix["sources"]
            if "gated" in s["status"] or "blocked" in s["status"]
        )
        man = self.write_manifest(
            n_sources=len(matrix["sources"]),
            n_open_or_cataloged=open_n,
            n_gated_or_blocked=gated_n,
            n_zenodo_files=len([m for m in zenodo_meta if not m.get("error")]),
            n_common78=len(common78),
            n_physionet_respir=physio.get("n_databases")
            or len(physio.get("databases") or []),
            clinical_breathomics_present=cbd_present,
            offline=offline,
        )
        return {
            "access_matrix": matrix_path,
            "physionet": physio_path,
            "common78": c78_path,
            "hbdb_eval_index": eval_path,
            "zenodo": zenodo_path,
            "manifest": man,
        }

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        matrix = json.loads((self.out_dir / "ACCESS_MATRIX.json").read_text())
        c78 = json.loads((self.out_dir / "common78_features.json").read_text())
        physio = json.loads((self.out_dir / "physionet_respir_catalog.json").read_text())
        eval_idx = json.loads((self.out_dir / "hbdb_eval_index.json").read_text())

        # Map common-78 CIDs onto atlas panel when possible
        from ..ingest.public_breath import VOC_PUBCHEM_CID

        cid_to_voc = {int(v): k for k, v in VOC_PUBCHEM_CID.items()}
        mapped = []
        for f in c78.get("features") or []:
            cid = f.get("pubchem_cid")
            voc = cid_to_voc.get(int(cid)) if cid is not None else None
            mapped.append({**f, "atlas_voc_id": voc})

        fused = {
            "version": "1.0.0",
            "access_matrix": matrix,
            "common78": {
                "n_features": c78.get("n_features"),
                "n_mapped_to_atlas_panel": sum(1 for m in mapped if m.get("atlas_voc_id")),
                "features": mapped,
                "doi": c78.get("doi"),
            },
            "hbdb_eval_index": {
                "n_json_records": eval_idx.get("n_json_records"),
                "compounds": eval_idx.get("compounds"),
                "doi": eval_idx.get("doi"),
            },
            "physionet_respir": {
                "n_databases": physio.get("n_databases")
                or len(physio.get("databases") or []),
                "databases": physio.get("databases") or [],
                "note": "Waveform / ventilation modality — not fused into VOC priors.",
            },
        }
        out = knowledge_dir / "datasource_alt_breath_sources.json"
        out.write_text(json.dumps(fused, indent=2))
        return {
            "path": str(out),
            "n_sources": len(matrix.get("sources") or []),
            "n_common78": c78.get("n_features"),
            "n_common78_mapped": fused["common78"]["n_mapped_to_atlas_panel"],
            "n_physionet": fused["physionet_respir"]["n_databases"],
        }
