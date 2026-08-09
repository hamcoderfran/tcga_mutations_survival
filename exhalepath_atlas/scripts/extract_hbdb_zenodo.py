#!/usr/bin/env python3
"""Extract HBDB 60-disease VOC associations from Zenodo SQL dump.

Source: https://zenodo.org/records/14958797 (hbdb2_wo_sentences.sql)
Live site https://hbdb.cmdm.tw/ is often Cloudflare-blocked; this dump is open.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


def parse_all_inserts(sql: str, table: str) -> list[dict[str, Any]]:
    marker = f"INSERT INTO `{table}` ("
    rows: list[dict[str, Any]] = []
    pos = 0
    while True:
        start = sql.find(marker, pos)
        if start < 0:
            break
        paren = sql.find("(", start)
        col_end = sql.find(")", paren)
        cols = [c.strip().strip("`") for c in sql[paren + 1 : col_end].split(",")]
        v = sql.find("VALUES", col_end)
        i = v + len("VALUES")
        n = len(sql)
        while i < n:
            while i < n and sql[i] in " \n\r\t,":
                i += 1
            if i >= n or sql[i] == ";":
                if i < n and sql[i] == ";":
                    i += 1
                break
            if sql[i] != "(":
                raise ValueError(f"bad token at {i}: {sql[i : i + 60]!r}")
            i += 1
            vals: list[Any] = []
            while True:
                while i < n and sql[i] in " \n\r\t":
                    i += 1
                if sql[i] == ")":
                    i += 1
                    break
                if sql[i] == "'":
                    i += 1
                    buf: list[str] = []
                    while i < n:
                        ch = sql[i]
                        if ch == "\\" and i + 1 < n:
                            nxt = sql[i + 1]
                            mapping = {
                                "n": "\n",
                                "r": "\r",
                                "t": "\t",
                                "0": "\0",
                                "b": "\b",
                                "Z": "\x1a",
                                "\\": "\\",
                                "'": "'",
                                '"': '"',
                            }
                            buf.append(mapping.get(nxt, nxt))
                            i += 2
                            continue
                        if ch == "'":
                            if i + 1 < n and sql[i + 1] == "'":
                                buf.append("'")
                                i += 2
                                continue
                            i += 1
                            break
                        buf.append(ch)
                        i += 1
                    vals.append("".join(buf))
                else:
                    j = i
                    while j < n and sql[j] not in ",)":
                        j += 1
                    token = sql[i:j].strip()
                    vals.append(None if token.upper() == "NULL" else token)
                    i = j
                while i < n and sql[i] in " \n\r\t":
                    i += 1
                if i < n and sql[i] == ",":
                    i += 1
                    continue
                if i < n and sql[i] == ")":
                    i += 1
                    break
            rows.append(dict(zip(cols, vals)))
        pos = i
    return rows


def extract_from_sql(sql_path: Path) -> dict[str, Any]:
    sql = sql_path.read_text(errors="replace")
    diseases = parse_all_inserts(sql, "diseases")
    rels = parse_all_inserts(sql, "compound_disease_relationships")
    compounds = parse_all_inserts(sql, "compounds")

    comp_by_id = {int(c["id"]): c for c in compounds}
    by_disease: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for r in rels:
        cid = int(r["compound_id"])
        did = int(r["disease_id"])
        c = comp_by_id.get(cid)
        if not c:
            continue
        aliases = [
            a.strip()
            for a in (c.get("alias") or "").replace("<br>", "\n").split("\n")
            if a.strip()
        ]
        by_disease[did].append(
            {
                "compound_id": cid,
                "name": c.get("name"),
                "pubchem_cid": c.get("pubchemcid") or None,
                "formula": c.get("formula") or None,
                "kegg": c.get("kegg") or None,
                "references_count": int(c.get("references_count") or 0),
                "aliases": aliases[:12],
                "evidence": "relationship",
            }
        )

    for did in list(by_disease):
        seen: set[int] = set()
        uniq: list[dict[str, Any]] = []
        for x in sorted(
            by_disease[did],
            key=lambda x: (-x["references_count"], (x["name"] or "").lower()),
        ):
            if x["compound_id"] in seen:
                continue
            seen.add(x["compound_id"])
            uniq.append(x)
        by_disease[did] = uniq

    out_diseases = []
    for d in sorted(diseases, key=lambda x: int(x["id"])):
        did = int(d["id"])
        comps = by_disease.get(did, [])
        out_diseases.append(
            {
                "hbdb_disease_id": did,
                "name": d.get("name"),
                "mesh_id": d.get("meshid") or None,
                "mesh_name": d.get("meshname") or None,
                "location": d.get("location"),
                "references_count": int(d.get("references_count") or 0),
                "n_compounds": len(comps),
                "compounds": comps,
            }
        )

    disease_vocs = {
        "version": "2.1.0",
        "source": "zenodo_hbdb2_wo_sentences",
        "zenodo_doi": "10.5281/zenodo.14958797",
        "zenodo_record": "14958797",
        "zenodo_file": "hbdb2_wo_sentences.sql",
        "hbdb_site": "https://hbdb.cmdm.tw/?disease_page=1&tab=disease",
        "hbdb_paper": "https://doi.org/10.1093/database/baz139",
        "access_note": (
            "Live hbdb.cmdm.tw is Cloudflare-blocked in many environments. "
            "Associations extracted from the open Zenodo SQL dump "
            "hbdb2_wo_sentences.sql (public HBDB content)."
        ),
        "n_diseases": len(out_diseases),
        "n_disease_compound_links": sum(d["n_compounds"] for d in out_diseases),
        "n_unique_compounds_linked": len(
            {c["compound_id"] for d in out_diseases for c in d["compounds"]}
        ),
        "n_compounds_catalog": len(compounds),
        "n_relationship_rows": len(rels),
        "diseases_with_zero_compounds": [
            d["hbdb_disease_id"] for d in out_diseases if d["n_compounds"] == 0
        ],
        "diseases": out_diseases,
    }

    cat = []
    for c in compounds:
        cat.append(
            {
                "compound_id": int(c["id"]),
                "name": c.get("name"),
                "pubchem_cid": c.get("pubchemcid") or None,
                "formula": c.get("formula") or None,
                "kegg": c.get("kegg") or None,
                "references_count": int(c.get("references_count") or 0),
                "aliases": [
                    a.strip()
                    for a in (c.get("alias") or "").replace("<br>", "\n").split("\n")
                    if a.strip()
                ][:20],
            }
        )
    compounds_doc = {
        "version": "2.1.0",
        "source": "zenodo_hbdb2_wo_sentences",
        "zenodo_doi": "10.5281/zenodo.14958797",
        "n_compounds": len(cat),
        "compounds": cat,
    }
    return {"disease_vocs": disease_vocs, "compounds": compounds_doc}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sql_path", type=Path, help="Path to hbdb2_wo_sentences.sql")
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("data/datasources/hbdb"),
        help="Output directory for JSON artifacts",
    )
    args = p.parse_args(argv)
    extracted = extract_from_sql(args.sql_path)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "hbdb_disease_vocs_60.json").write_text(
        json.dumps(extracted["disease_vocs"], indent=2) + "\n"
    )
    (args.out_dir / "hbdb_compounds.json").write_text(
        json.dumps(extracted["compounds"], indent=2) + "\n"
    )
    # Rebuild atlas mapping if voc catalog present
    root = Path(__file__).resolve().parents[1]
    voc_path = root / "data" / "knowledge" / "voc_catalog.json"
    if voc_path.exists():
        sys.path.insert(0, str(root / "src"))
        from exhalepath.datasources.hbdb_mapping import build_atlas_mapped_doc

        atlas_vocs = json.loads(voc_path.read_text()).get("vocs") or []
        mapped = build_atlas_mapped_doc(
            extracted["disease_vocs"], extracted["compounds"], atlas_vocs
        )
        (args.out_dir / "hbdb_atlas_mapped.json").write_text(
            json.dumps(mapped, indent=2) + "\n"
        )
        print(
            "wrote",
            args.out_dir,
            "diseases",
            mapped["n_hbdb_diseases"],
            "links",
            mapped["n_disease_compound_links"],
            "panel links",
            mapped["n_links_mapped_to_atlas_panel"],
        )
    else:
        print("wrote disease_vocs + compounds (no voc_catalog for mapping)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
