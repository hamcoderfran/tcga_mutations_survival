"""HBDB 60-disease Zenodo extract + atlas mapping."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.datasources.ds13_hbdb import HBDBSource, load_hbdb_atlas_associations
from exhalepath.datasources.hbdb_mapping import HBDB_DISEASE_TO_ATLAS, build_atlas_mapped_doc
from exhalepath.knowledge.external_evidence import collect_external_evidence

ROOT = Path(__file__).resolve().parents[1]
HBDB = ROOT / "data" / "datasources" / "hbdb"


def test_hbdb_sixty_disease_extract_present():
    raw = json.loads((HBDB / "hbdb_disease_vocs_60.json").read_text())
    assert raw["n_diseases"] == 60
    assert raw["n_disease_compound_links"] >= 700
    assert len(raw["diseases"]) == 60
    # Eleven diseases have no compound_disease_relationships rows in the dump
    assert len(raw.get("diseases_with_zero_compounds") or []) == 11
    asthma = next(d for d in raw["diseases"] if d["hbdb_disease_id"] == 5)
    assert asthma["n_compounds"] >= 150
    copd = next(d for d in raw["diseases"] if d["hbdb_disease_id"] == 57)
    assert copd["n_compounds"] >= 100


def test_hbdb_atlas_mapping_rebuilds():
    raw = json.loads((HBDB / "hbdb_disease_vocs_60.json").read_text())
    compounds = json.loads((HBDB / "hbdb_compounds.json").read_text())
    vocs = json.loads((ROOT / "data" / "knowledge" / "voc_catalog.json").read_text())["vocs"]
    mapped = build_atlas_mapped_doc(raw, compounds, vocs)
    assert mapped["n_hbdb_diseases"] == 60
    assert mapped["n_links_mapped_to_atlas_panel"] >= 30
    assert mapped["n_atlas_diseases_with_panel_vocs"] >= 6
    assert "asthma" in {a["disease_id"] for a in mapped["atlas_associations"]}
    assert "copd" in {a["disease_id"] for a in mapped["atlas_associations"]}
    # Mapping table covers every HBDB disease id present in extract
    ids = {d["hbdb_disease_id"] for d in raw["diseases"]}
    assert ids <= set(HBDB_DISEASE_TO_ATLAS)


def test_hbdb_harvest_offline_uses_zenodo_extract():
    src = HBDBSource(ROOT)
    paths = src.harvest(offline=True)
    assert paths["disease_associations"].exists()
    doc = json.loads(paths["disease_associations"].read_text())
    assert doc.get("source") == "hbdb_zenodo_sql"
    assert (doc.get("n_disease_compound_links") or 0) >= 700
    info = src.fuse(ROOT / "data" / "knowledge")
    assert (info.get("n_hbdb_links") or 0) >= 700
    fused = json.loads((ROOT / "data" / "knowledge" / "datasource_hbdb.json").read_text())
    assert fused.get("hbdb_extract", {}).get("n_hbdb_diseases") == 60


def test_external_evidence_uses_hbdb_sql_not_only_proxy():
    assoc = load_hbdb_atlas_associations(ROOT)
    assert assoc
    vocs = {
        v["voc_id"]
        for v in json.loads((ROOT / "data" / "knowledge" / "voc_catalog.json").read_text())[
            "vocs"
        ]
    }
    ids = {
        d["disease_id"]
        for d in json.loads((ROOT / "data" / "knowledge" / "disease_voc_priors.json").read_text())[
            "diseases"
        ]
    }
    packed = collect_external_evidence(
        vocs,
        ids,
        include_lit_bench=False,
        include_public_breath=False,
        include_panels=False,
        include_hbdb=True,
    )
    assert packed["stats"]["n_hbdb_edges"] >= 20
    # Sources should prefer zenodo sql tag
    sources = {
        s
        for did in packed["evidence"].values()
        for slot in did.values()
        for s in slot.get("sources") or []
    }
    assert "hbdb_zenodo_sql" in sources
