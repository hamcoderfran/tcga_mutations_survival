"""Map HBDB disease/compound records onto ExhalePath atlas IDs.

HBDB live site (hbdb.cmdm.tw) is Cloudflare-blocked in many environments.
Canonical associations ship from the open Zenodo dump
``hbdb2_wo_sentences.sql`` (10.5281/zenodo.14958797).
"""

from __future__ import annotations

import re
from typing import Any

# HBDB disease_id → atlas disease_id (None = inventory-only, no prior fuse)
HBDB_DISEASE_TO_ATLAS: dict[int, str | None] = {
    1: "alcohol_use",
    2: None,  # Acute Mountain Sickness
    3: "ild",
    4: "ild",
    5: "asthma",
    6: "breast_invasive_carcinoma",
    7: "chronic_bronchitis",  # bronchiectasis → nearest airway atlas id
    8: "chronic_bronchitis",
    9: "lung_adenocarcinoma",
    10: "asthma",
    11: None,
    12: None,
    13: "cystic_fibrosis",
    14: "atopic_dermatitis",
    15: "type1_diabetes",
    16: None,
    17: "atopic_dermatitis",
    18: None,
    19: "gerd",
    20: None,
    21: "heart_failure",
    22: "chronic_kidney_disease",
    23: "ild",
    24: "pulmonary_hypertension",
    25: None,
    26: None,
    27: "chronic_kidney_disease",
    28: "copd",
    29: "lung_adenocarcinoma",
    30: "sleep_apnea",
    31: "ild",
    32: None,
    33: "pneumonia_bacterial",
    35: None,
    36: "ards",
    37: "asthma",
    38: None,
    39: "asthma",
    40: "ild",
    41: "scleroderma",
    42: "ild",
    43: "sleep_apnea",
    44: "tobacco_use",
    45: "asthma",
    46: None,
    47: "tuberculosis",
    48: None,
    49: "inflammatory_bowel_disease",
    50: "vasculitis",
    51: None,
    52: None,
    53: "asthma",
    54: "ild",
    55: "pneumonia_bacterial",
    56: "sleep_apnea",
    57: "copd",
    58: "ild",
    59: "asthma",
    60: "asthma",
    61: None,
}

# Extra name keys (normalized) → atlas voc_id for HBDB IUPAC / alias forms
_VOC_EXTRA_NAMES: dict[str, list[str]] = {
    "isoprene": ["2 methylbuta 1 3 diene", "2 methyl 1 3 butadiene"],
    "2_butanone": ["butan 2 one", "methyl ethyl ketone", "mek"],
    "dms": ["methylsulfanylmethane", "dimethyl sulfide", "dimethyl sulphide"],
    "limonene": ["1 methyl 4 prop 1 en 2 ylcyclohexene", "d limonene"],
    "carbon_disulfide": ["methanedithione", "carbon disulphide"],
    "trimethylamine": ["n n dimethylmethanamine", "trimethyl amine"],
    "methyl_mercaptan": ["methanethiol"],
    "allyl_methyl_sulfide": ["3 methylsulfanylprop 1 ene"],
    "pyrrole": ["1h pyrrole"],
    "crotonaldehyde": ["e but 2 enal", "but 2 enal"],
    "propanol": ["propan 1 ol", "1 propanol", "n propanol"],
    "propionaldehyde": ["propanal"],
    "hexanal": ["hexaldehyde", "n hexanal"],
    "xylene": [
        "p xylene",
        "m xylene",
        "o xylene",
        "1 2 dimethylbenzene",
        "1 3 dimethylbenzene",
        "1 4 dimethylbenzene",
    ],
    "dimethyl_disulfide": ["methyldisulfanylmethane"],
    "hydrogen_sulfide": ["hydrogen sulphide", "sulfane"],
    "dmts": ["dimethyl trisulfide"],
    "dimethyl_amine": ["n methylmethanamine", "dimethylamine"],
    "indole": ["1h indole"],
}


def norm_name(s: str) -> str:
    s = str(s or "").lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def build_compound_id_to_voc(
    hbdb_compounds: list[dict[str, Any]],
    atlas_vocs: list[dict[str, Any]],
) -> dict[int, str]:
    """Map HBDB compound_id → atlas voc_id via names/aliases."""
    name_to_cids: dict[str, set[int]] = {}
    for c in hbdb_compounds:
        cid = int(c["compound_id"] if "compound_id" in c else c["id"])
        names = [c.get("name") or ""] + list(c.get("aliases") or [])
        for n in names:
            k = norm_name(n)
            if not k:
                continue
            name_to_cids.setdefault(k, set()).add(cid)

    cid_to_voc: dict[int, str] = {}
    for v in atlas_vocs:
        voc_id = v["voc_id"]
        keys = [norm_name(v.get("name") or ""), norm_name(voc_id.replace("_", " "))]
        for a in v.get("aliases") or []:
            keys.append(norm_name(a))
        for a in _VOC_EXTRA_NAMES.get(voc_id, []):
            keys.append(norm_name(a))
        for k in keys:
            if not k:
                continue
            for cid in name_to_cids.get(k, ()):
                cid_to_voc[cid] = voc_id
    return cid_to_voc


def build_atlas_mapped_doc(
    disease_vocs: dict[str, Any],
    compounds_doc: dict[str, Any],
    atlas_vocs: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build ``hbdb_atlas_mapped.json`` content from raw HBDB extracts."""
    cid_to_voc = build_compound_id_to_voc(
        list(compounds_doc.get("compounds") or []),
        atlas_vocs,
    )
    atlas_edges: dict[str, dict[str, list[dict[str, Any]]]] = {}
    mapped_diseases: list[dict[str, Any]] = []
    n_total = 0
    n_mapped = 0

    for d in disease_vocs.get("diseases") or []:
        hid = int(d["hbdb_disease_id"])
        atlas = HBDB_DISEASE_TO_ATLAS.get(hid)
        comps_out: list[dict[str, Any]] = []
        for c in d.get("compounds") or []:
            n_total += 1
            cid = int(c["compound_id"])
            vid = cid_to_voc.get(cid)
            comps_out.append(
                {
                    "compound_id": cid,
                    "name": c.get("name"),
                    "pubchem_cid": c.get("pubchem_cid"),
                    "formula": c.get("formula"),
                    "references_count": c.get("references_count"),
                    "evidence": c.get("evidence"),
                    "atlas_voc_id": vid,
                }
            )
            if not vid:
                continue
            n_mapped += 1
            if not atlas:
                continue
            atlas_edges.setdefault(atlas, {}).setdefault(vid, []).append(
                {
                    "hbdb_disease_id": hid,
                    "hbdb_disease_name": d.get("name"),
                    "compound_id": cid,
                    "compound_name": c.get("name"),
                    "pubchem_cid": c.get("pubchem_cid"),
                    "references_count": c.get("references_count"),
                }
            )
        mapped_diseases.append(
            {
                "hbdb_disease_id": hid,
                "name": d.get("name"),
                "mesh_id": d.get("mesh_id"),
                "mesh_name": d.get("mesh_name"),
                "location": d.get("location"),
                "references_count": d.get("references_count"),
                "atlas_disease_id": atlas,
                "n_compounds": d.get("n_compounds") or len(comps_out),
                "n_mapped_to_atlas_panel": sum(1 for x in comps_out if x["atlas_voc_id"]),
                "compounds": comps_out,
            }
        )

    atlas_assoc: list[dict[str, Any]] = []
    for did in sorted(atlas_edges):
        voc_map: dict[str, Any] = {}
        for voc, sources in atlas_edges[did].items():
            sources_sorted = sorted(
                sources,
                key=lambda s: (-(s.get("references_count") or 0), s["compound_id"]),
            )
            voc_map[voc] = {
                "n_hbdb_links": len(sources_sorted),
                "sources": sources_sorted[:8],
            }
        atlas_assoc.append(
            {
                "disease_id": did,
                "n_atlas_vocs": len(voc_map),
                "vocs": voc_map,
                "hbdb_style_voc_names": sorted(voc_map.keys()),
                "source": "hbdb_zenodo_sql",
            }
        )

    return {
        "version": "2.1.0",
        "description": (
            "HBDB all-60-disease compound associations from Zenodo SQL, mapped onto "
            "ExhalePath atlas disease_id and 50-VOC prediction panel where possible."
        ),
        "source": disease_vocs.get("source") or "zenodo_hbdb2_wo_sentences",
        "zenodo_doi": disease_vocs.get("zenodo_doi") or "10.5281/zenodo.14958797",
        "zenodo_file": disease_vocs.get("zenodo_file") or "hbdb2_wo_sentences.sql",
        "hbdb_site": "https://hbdb.cmdm.tw/?disease_page=1&tab=disease",
        "hbdb_paper": disease_vocs.get("hbdb_paper")
        or "https://doi.org/10.1093/database/baz139",
        "access_note": disease_vocs.get("access_note"),
        "n_hbdb_diseases": len(mapped_diseases),
        "n_hbdb_diseases_mapped_to_atlas": sum(
            1 for d in mapped_diseases if d["atlas_disease_id"]
        ),
        "n_disease_compound_links": n_total,
        "n_links_mapped_to_atlas_panel": n_mapped,
        "n_atlas_diseases_with_panel_vocs": len(atlas_assoc),
        "n_hbdb_compound_ids_on_panel": len(cid_to_voc),
        "disease_id_map": {str(k): v for k, v in sorted(HBDB_DISEASE_TO_ATLAS.items())},
        "compound_id_to_atlas_voc": {str(k): v for k, v in sorted(cid_to_voc.items())},
        "diseases": mapped_diseases,
        "atlas_associations": atlas_assoc,
    }
