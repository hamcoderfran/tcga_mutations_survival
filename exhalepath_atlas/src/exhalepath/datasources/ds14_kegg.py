"""Priority 14 — KEGG VOC compound → reaction → pathway → enzyme maps."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote

import requests

from .base import DataSource

UA = {"User-Agent": "ExhalePathAtlas/1.0 (academic research; contact: exhalepath)"}
KEGG = "https://rest.kegg.jp"
SLEEP = 0.35  # stay under ~3 requests/sec

# High-confidence KEGG compound IDs for catalog VOCs (skip ambiguous find)
VOC_KEGG_SEED = {
    "acetone": "C00207",
    "isoprene": "C16521",
    "ethanol": "C00469",
    "methanol": "C00132",
    "ammonia": "C00014",
    "acetaldehyde": "C00084",
    "formaldehyde": "C00067",
    "hydrogen_sulfide": "C00283",
    "phenol": "C00146",
    "indole": "C00463",
    "toluene": "C01455",
    "benzene": "C01407",
    "trimethylamine": "C00565",
    "dimethyl_amine": "C00543",
    "propanol": "C00483",
    "isopropanol": "C01845",
    "2_butanone": "C02845",
    "benzaldehyde": "C00261",
    "limonene": "C06099",
    "dms": "C00580",
    "methyl_mercaptan": "C00409",
    "carbon_disulfide": "C18904",
    "pyrrole": "C05085",
    "furan": "C14248",
    "acetonitrile": "C01548",
    "propionaldehyde": "C00466",
    "hexanal": "C02233",
    "pentane": "C13388",
}


def _parse_link(text: str) -> list[str]:
    out = []
    for line in text.strip().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out.append(parts[1].split(":", 1)[-1] if ":" in parts[1] else parts[1])
    return out


def _parse_find_compound(text: str) -> str | None:
    for line in text.strip().splitlines():
        # cpd:C00207\tAcetone; ...
        m = re.match(r"cpd:(C\d+)\t", line)
        if m:
            return m.group(1)
    return None


def _parse_get_compound(text: str) -> dict[str, Any]:
    name = None
    formula = None
    pathways: list[dict[str, str]] = []
    enzymes: list[str] = []
    reactions: list[str] = []
    for line in text.splitlines():
        if line.startswith("NAME"):
            name = line.split(None, 1)[1].rstrip(";").strip() if len(line.split(None, 1)) > 1 else None
        elif line.startswith("FORMULA"):
            formula = line.split(None, 1)[1].strip()
        elif line.startswith("PATHWAY"):
            rest = line.split(None, 1)[1]
            pid = rest.split(None, 1)[0]
            pname = rest.split(None, 1)[1] if len(rest.split(None, 1)) > 1 else pid
            pathways.append({"kegg_pathway_id": pid, "name": pname})
        elif line.startswith("            map") or (
            line.startswith("            ") and pathways and line.strip().startswith("map")
        ):
            rest = line.strip()
            pid = rest.split(None, 1)[0]
            pname = rest.split(None, 1)[1] if len(rest.split(None, 1)) > 1 else pid
            pathways.append({"kegg_pathway_id": pid, "name": pname})
        elif line.startswith("ENZYME"):
            enzymes.extend(re.findall(r"\d+\.\d+\.\d+\.\d+", line))
            enzymes.extend(re.findall(r"\d+\.\d+\.\d+\.-", line))
        elif line.startswith("            ") and "ENZYME" in text[: text.find(line) + 1]:
            enzymes.extend(re.findall(r"\d+\.\d+\.\d+\.\d+", line))
        elif line.startswith("REACTION"):
            reactions.extend(re.findall(r"R\d+", line))
        elif line.startswith("            ") and reactions is not None:
            reactions.extend(re.findall(r"R\d+", line))
    # Dedupe
    seen_p = set()
    uniq_pw = []
    for p in pathways:
        if p["kegg_pathway_id"] in seen_p:
            continue
        seen_p.add(p["kegg_pathway_id"])
        uniq_pw.append(p)
    return {
        "name": name,
        "formula": formula,
        "pathways": uniq_pw,
        "enzymes": list(dict.fromkeys(enzymes)),
        "reactions": list(dict.fromkeys(reactions)),
    }


class KEGGSource(DataSource):
    priority = 14
    key = "kegg"
    title = "KEGG VOC reaction / pathway / enzyme maps"
    description = "Map ExhalePath VOCs to KEGG compounds, reactions, human pathways, and EC numbers"

    def _get(self, path: str) -> str:
        r = requests.get(f"{KEGG}/{path}", timeout=45, headers=UA)
        time.sleep(SLEEP)
        r.raise_for_status()
        return r.text

    def _resolve_compound(self, voc_id: str, display_name: str) -> str | None:
        if voc_id in VOC_KEGG_SEED:
            return VOC_KEGG_SEED[voc_id]
        # try find by name
        q = display_name.split("(")[0].strip()
        try:
            text = self._get(f"find/compound/{quote(q)}")
            return _parse_find_compound(text)
        except Exception:  # noqa: BLE001
            return None

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        cache = self.out_dir / "kegg_voc_maps.json"
        if offline and cache.exists():
            doc = json.loads(cache.read_text())
            man = self.write_manifest(n_vocs=doc.get("n_mapped"), offline=True)
            return {"maps": cache, "manifest": man}

        voc_path = self.root / "data" / "knowledge" / "voc_catalog.json"
        vocs = []
        if voc_path.exists():
            vocs = json.loads(voc_path.read_text()).get("vocs") or []
        if not vocs:
            vocs = [{"voc_id": k, "name": k} for k in VOC_KEGG_SEED]

        rows: list[dict[str, Any]] = []
        if offline:
            # Seed-only offline stub
            for vid, cid in VOC_KEGG_SEED.items():
                rows.append(
                    {
                        "voc_id": vid,
                        "kegg_compound_id": cid,
                        "pathways": [],
                        "enzymes": [],
                        "reactions": [],
                        "offline": True,
                    }
                )
        else:
            for v in vocs:
                vid = v["voc_id"]
                name = v.get("name") or vid
                row: dict[str, Any] = {"voc_id": vid, "name": name}
                try:
                    cid = self._resolve_compound(vid, name)
                    row["kegg_compound_id"] = cid
                    if not cid:
                        row["status"] = "unmapped"
                        rows.append(row)
                        continue
                    text = self._get(f"get/cpd:{cid}")
                    parsed = _parse_get_compound(text)
                    row.update(parsed)
                    # Extra link endpoints for completeness
                    try:
                        rxn = _parse_link(self._get(f"link/reaction/cpd:{cid}"))
                        if rxn:
                            row["reactions"] = list(
                                dict.fromkeys(list(row.get("reactions") or []) + rxn)
                            )
                    except Exception:  # noqa: BLE001
                        pass
                    try:
                        enz = _parse_link(self._get(f"link/enzyme/cpd:{cid}"))
                        if enz:
                            row["enzymes"] = list(
                                dict.fromkeys(list(row.get("enzymes") or []) + enz)
                            )
                    except Exception:  # noqa: BLE001
                        pass
                    row["status"] = "mapped"
                    row["kegg_url"] = f"https://www.kegg.jp/entry/{cid}"
                except Exception as e:  # noqa: BLE001
                    row["error"] = str(e)
                    row["status"] = "error"
                rows.append(row)

        doc = {
            "version": "1.0.0",
            "source": "KEGG REST API",
            "license_note": "KEGG academic use; cite Kanehisa Laboratories",
            "n_vocs": len(rows),
            "n_mapped": sum(1 for r in rows if r.get("kegg_compound_id")),
            "vocs": rows,
        }
        path = self.write_json("kegg_voc_maps.json", doc)
        man = self.write_manifest(n_vocs=doc["n_vocs"], n_mapped=doc["n_mapped"])
        return {"maps": path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "kegg_voc_maps.json").read_text())
        out = knowledge_dir / "datasource_kegg.json"
        out.write_text(json.dumps(doc, indent=2))

        # Attach KEGG IDs onto voc_catalog
        voc_path = knowledge_dir / "voc_catalog.json"
        n = 0
        if voc_path.exists():
            by_id = {r["voc_id"]: r for r in doc.get("vocs") or []}
            voc_doc = json.loads(voc_path.read_text())
            for v in voc_doc.get("vocs") or []:
                hit = by_id.get(v["voc_id"])
                if not hit or not hit.get("kegg_compound_id"):
                    continue
                v["kegg_compound_id"] = hit["kegg_compound_id"]
                v["kegg_pathways"] = [
                    p.get("kegg_pathway_id") for p in (hit.get("pathways") or [])[:12]
                ]
                v["kegg_enzymes"] = (hit.get("enzymes") or [])[:20]
                n += 1
            voc_path.write_text(json.dumps(voc_doc, indent=2))

        # Enrich pathway_voc_map with KEGG pathway crossrefs when names overlap
        pw_path = knowledge_dir / "pathway_voc_map.json"
        n_pw = 0
        if pw_path.exists():
            pw = json.loads(pw_path.read_text())
            # Build voc → kegg pathways
            voc_kegg_pw: dict[str, list[dict[str, str]]] = {}
            for r in doc.get("vocs") or []:
                if r.get("pathways"):
                    voc_kegg_pw[r["voc_id"]] = r["pathways"]
            for p in pw.get("pathways") or []:
                linked = []
                voc_ids = list((p.get("voc_effects") or {}).keys())
                for vid in voc_ids:
                    for kp in voc_kegg_pw.get(vid) or []:
                        linked.append(kp)
                if not linked:
                    # match by pathway name tokens against KEGG pathway names
                    pname = _norm_simple(p.get("name") or p.get("pathway_id") or "")
                    for r in doc.get("vocs") or []:
                        for kp in r.get("pathways") or []:
                            kn = _norm_simple(kp.get("name") or "")
                            if pname and (pname in kn or kn in pname):
                                linked.append(kp)
                if linked:
                    # unique
                    seen = set()
                    uniq = []
                    for kp in linked:
                        kid = kp.get("kegg_pathway_id")
                        if kid in seen:
                            continue
                        seen.add(kid)
                        uniq.append(kp)
                    p["kegg_pathways"] = uniq[:15]
                    n_pw += 1
            pw_path.write_text(json.dumps(pw, indent=2))

        return {
            "path": str(out),
            "n_mapped": doc.get("n_mapped"),
            "n_voc_catalog_enriched": n,
            "n_pathways_enriched": n_pw,
        }


def _norm_simple(s: str) -> str:
    s = s.lower().replace("_", " ")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()
