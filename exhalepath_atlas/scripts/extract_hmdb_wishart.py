#!/usr/bin/env python3
"""Extract breath + atlas-panel metabolites from HMDB via Wishart fix mirror.

``hmdb.ca`` is Cloudflare-blocked in many environments. The Wishart lab mirror
``hmdbfix.wishartlab.com`` serves the same ``hmdb_metabolites.zip`` openly.

Usage:
  python scripts/extract_hmdb_wishart.py /tmp/hmdb_metabolites.zip
  python scripts/extract_hmdb_wishart.py --download
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import iterparse

import requests

MIRROR_ZIP = (
    "https://hmdbfix.wishartlab.com/system/downloads/current/hmdb_metabolites.zip"
)
MIRROR_METABOLITE = "https://hmdbfix.wishartlab.com/metabolites/{hmdb_id}.xml"


def local(tag: str) -> str:
    return tag.split("}", 1)[-1] if "}" in tag else tag


def text(el: Any) -> str:
    return (el.text or "").strip() if el is not None else ""


def child_text(el: Any, name: str) -> str:
    for c in list(el):
        if local(c.tag) == name:
            return text(c)
    return ""


def extract(zip_path: Path, *, voc_hmdb: dict[str, str], voc_cid: dict[str, int]) -> dict[str, Any]:
    wanted_cid = {int(v): k for k, v in voc_cid.items()}
    voc_by_hmdb = {v.upper(): k for k, v in voc_hmdb.items()}

    breath: list[dict[str, Any]] = []
    panel: dict[str, dict[str, Any]] = {}
    biofluid_counts: Counter[str] = Counter()
    n_metab = 0

    with zipfile.ZipFile(zip_path) as zf, zf.open("hmdb_metabolites.xml") as fh:
        for _event, elem in iterparse(fh, events=("end",)):
            if local(elem.tag) != "metabolite":
                continue
            n_metab += 1
            acc = child_text(elem, "accession").upper()
            name = child_text(elem, "name")
            pc = None
            for c in elem.iter():
                if local(c.tag) == "pubchem_compound_id" and text(c):
                    try:
                        pc = int(text(c))
                    except ValueError:
                        pc = None
                    break
            bios = sorted(
                {
                    text(c)
                    for c in elem.iter()
                    if local(c.tag) == "biospecimen" and text(c)
                }
            )
            for b in bios:
                biofluid_counts[b] += 1
            formula = child_text(elem, "chemical_formula") or None
            inchikey = child_text(elem, "inchikey") or None
            smiles = child_text(elem, "smiles") or None
            cas = child_text(elem, "cas_registry_number") or None
            avg_mass = child_text(elem, "average_molecular_weight")
            desc = child_text(elem, "description")
            if len(desc) > 400:
                desc = desc[:400] + "…"
            syns: list[str] = []
            for c in elem.iter():
                if local(c.tag) == "synonym" and text(c):
                    syns.append(text(c))
                    if len(syns) >= 15:
                        break

            voc = voc_by_hmdb.get(acc)
            if voc is None and pc in wanted_cid:
                voc = wanted_cid[pc]

            rec = {
                "hmdb_id": acc,
                "name": name,
                "pubchem_cid": pc,
                "formula": formula,
                "inchikey": inchikey,
                "smiles": smiles,
                "cas": cas,
                "average_mass": float(avg_mass) if avg_mass else None,
                "biospecimens": bios,
                "synonyms": syns,
                "description": desc or None,
                "in_breath": any(b.lower() == "breath" for b in bios),
                "in_blood": any(b.lower() == "blood" for b in bios),
                "in_urine": any(b.lower() == "urine" for b in bios),
                "in_saliva": any(b.lower() == "saliva" for b in bios),
            }
            if voc is not None:
                rec["atlas_voc_id"] = voc
                panel[voc] = rec
            if rec["in_breath"]:
                breath.append(
                    {
                        "hmdb_id": acc,
                        "name": name,
                        "pubchem_cid": pc,
                        "formula": formula,
                        "inchikey": inchikey,
                        "cas": cas,
                        "biospecimens": bios,
                        "atlas_voc_id": voc,
                    }
                )
            elem.clear()
            if n_metab % 50000 == 0:
                print(f"… {n_metab} breath={len(breath)} panel={len(panel)}", flush=True)

    return {
        "n_metabolites_scanned": n_metab,
        "biofluid_counts": dict(biofluid_counts.most_common()),
        "breath": breath,
        "panel": panel,
        "missing_panel_vocs": sorted(set(voc_hmdb) - set(panel)),
    }


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("zip_path", nargs="?", type=Path, help="Local hmdb_metabolites.zip")
    p.add_argument("--download", action="store_true", help="Download zip from Wishart mirror")
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/datasources/hmdb"),
        help="Output directory",
    )
    args = p.parse_args(argv)

    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "src"))
    from exhalepath.datasources.ds02_hmdb import VOC_HMDB
    from exhalepath.ingest.public_breath import VOC_PUBCHEM_CID

    zip_path = args.zip_path
    if args.download or zip_path is None:
        zip_path = Path("/tmp/hmdb_metabolites.zip")
        print(f"Downloading {MIRROR_ZIP} → {zip_path}", flush=True)
        r = requests.get(MIRROR_ZIP, timeout=600, headers={"User-Agent": "ExhalePathAtlas/1.0"})
        r.raise_for_status()
        zip_path.write_bytes(r.content)
        print(f"saved {zip_path.stat().st_size} bytes", flush=True)

    extracted = extract(zip_path, voc_hmdb=VOC_HMDB, voc_cid=VOC_PUBCHEM_CID)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    breath_doc = {
        "version": "1.0.0",
        "source": "hmdbfix.wishartlab.com",
        "source_url": MIRROR_ZIP,
        "mirror_metabolite_xml": MIRROR_METABOLITE,
        "hmdb_site": "https://hmdb.ca/",
        "access_note": (
            "Live hmdb.ca is Cloudflare-blocked in many automated environments. "
            "Bulk XML pulled from the Wishart lab fix mirror hmdbfix.wishartlab.com."
        ),
        "license": "CC-BY-NC-4.0",
        "citation": "Wishart DS et al. HMDB 5.0. Nucleic Acids Res. 2022;50(D1):D622-D631.",
        "n_metabolites_scanned": extracted["n_metabolites_scanned"],
        "n_breath": len(extracted["breath"]),
        "biofluid_counts": extracted["biofluid_counts"],
        "metabolites": extracted["breath"],
    }
    panel_doc = {
        "version": "1.0.0",
        "source": "hmdbfix.wishartlab.com",
        "source_url": MIRROR_ZIP,
        "license": "CC-BY-NC-4.0",
        "n_panel_vocs": len(VOC_HMDB),
        "n_enriched": len(extracted["panel"]),
        "missing_vocs": extracted["missing_panel_vocs"],
        "by_voc": extracted["panel"],
    }
    (args.out_dir / "hmdb_breath_metabolites.json").write_text(
        json.dumps(breath_doc, indent=2) + "\n"
    )
    (args.out_dir / "hmdb_panel_enrichment.json").write_text(
        json.dumps(panel_doc, indent=2) + "\n"
    )
    print(
        f"wrote breath={breath_doc['n_breath']} panel={panel_doc['n_enriched']} "
        f"scanned={breath_doc['n_metabolites_scanned']} → {args.out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
