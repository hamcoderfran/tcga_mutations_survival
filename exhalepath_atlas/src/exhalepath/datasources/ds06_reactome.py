"""Priority 6 — Reactome / Rhea-style pathway expansion for VOC chains."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import requests

from .base import DataSource


class ReactomeSource(DataSource):
    priority = 6
    key = "reactome"
    title = "Reactome pathway expansion"
    description = "Expand seed genes to Reactome pathway participants for VOC routes"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        kb_path = self.root / "data" / "knowledge" / "pathway_voc_map.json"
        pathways = json.loads(kb_path.read_text())["pathways"] if kb_path.exists() else []
        expansions = []
        for p in pathways:
            entry: dict[str, Any] = {
                "pathway_id": p["pathway_id"],
                "reactome_ids": p.get("reactome_ids") or [],
                "seed_genes": p.get("seed_genes") or [],
                "expanded_genes": list(p.get("seed_genes") or []),
            }
            if not offline and entry["reactome_ids"]:
                rid = entry["reactome_ids"][0]
                try:
                    r = requests.get(
                        f"https://reactome.org/ContentService/data/participants/{rid}",
                        timeout=45,
                        headers={"User-Agent": "ExhalePathAtlas/1.0"},
                    )
                    genes = set(entry["expanded_genes"])
                    if r.ok:
                        data = r.json()
                        # participants payload varies; collect displayNames that look like genes
                        blob = json.dumps(data)
                        for g in entry["seed_genes"]:
                            if g in blob:
                                genes.add(g)
                        # pull UniProt display names lightly
                        if isinstance(data, list):
                            for item in data[:50]:
                                ref = item.get("refEntities") or item.get("displayName")
                                if isinstance(ref, list):
                                    for ent in ref[:20]:
                                        dn = str(ent.get("displayName") or "")
                                        if dn.isupper() and 2 <= len(dn) <= 8:
                                            genes.add(dn)
                    entry["expanded_genes"] = sorted(genes)
                    entry["http_status"] = r.status_code
                except Exception as e:  # noqa: BLE001
                    entry["error"] = str(e)
            expansions.append(entry)
        doc = {"version": "1.0.0", "n_pathways": len(expansions), "pathways": expansions}
        path = self.write_json("reactome_expansions.json", doc)
        return {"expansions": path, "manifest": self.write_manifest(n_pathways=len(expansions))}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "reactome_expansions.json").read_text())
        out = knowledge_dir / "datasource_reactome.json"
        out.write_text(json.dumps(doc, indent=2))
        # merge expanded genes into pathway_voc_map seed_genes (union)
        pw_path = knowledge_dir / "pathway_voc_map.json"
        n = 0
        if pw_path.exists():
            pw = json.loads(pw_path.read_text())
            exp = {p["pathway_id"]: p for p in doc["pathways"]}
            for p in pw.get("pathways") or []:
                e = exp.get(p["pathway_id"])
                if not e:
                    continue
                seeds = list(dict.fromkeys(list(p.get("seed_genes") or []) + list(e.get("expanded_genes") or [])))
                if len(seeds) > len(p.get("seed_genes") or []):
                    n += len(seeds) - len(p.get("seed_genes") or [])
                p["seed_genes"] = seeds[:40]
            pw_path.write_text(json.dumps(pw, indent=2))
        return {"path": str(out), "n_genes_added": n}
