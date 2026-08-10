"""HMDB via Wishart fix mirror (hmdb.ca Cloudflare bypass)."""

from __future__ import annotations

import json
from pathlib import Path

import requests

from exhalepath.datasources.ds02_hmdb import HMDB_MIRROR, HMDBSource

ROOT = Path(__file__).resolve().parents[1]
HMDB = ROOT / "data" / "datasources" / "hmdb"


def test_wishart_mirror_serves_metabolite_xml():
    try:
        r = requests.get(
            f"{HMDB_MIRROR}/metabolites/HMDB0001659.xml",
            timeout=45,
            headers={"User-Agent": "ExhalePathAtlas/1.0"},
        )
    except requests.RequestException as exc:
        import pytest

        pytest.skip(f"Wishart mirror unreachable: {exc}")
    if r.status_code != 200:
        import pytest

        pytest.skip(f"Wishart mirror status={r.status_code}")
    assert "<metabolite>" in r.text
    assert "Acetone" in r.text
    assert ">Breath<" in r.text


def test_hmdb_breath_extract_committed():
    breath = json.loads((HMDB / "hmdb_breath_metabolites.json").read_text())
    assert breath["n_breath"] == 60
    assert breath["n_metabolites_scanned"] >= 200_000
    assert breath["source"] == "hmdbfix.wishartlab.com"
    ids = {m["hmdb_id"] for m in breath["metabolites"]}
    assert "HMDB0001659" in ids  # acetone


def test_hmdb_panel_enrichment_complete():
    panel = json.loads((HMDB / "hmdb_panel_enrichment.json").read_text())
    assert panel["n_enriched"] >= 50
    assert "butyric_acid" in panel["by_voc"]
    assert panel["by_voc"]["butyric_acid"]["hmdb_id"] == "HMDB0000039"
    assert panel["by_voc"]["butyric_acid"]["in_breath"] is True
    assert not panel.get("missing_vocs")
    assert panel["by_voc"]["acetone"]["in_breath"] is True


def test_hmdb_harvest_uses_wishart_bulk():
    src = HMDBSource(ROOT)
    paths = src.harvest(offline=True)
    doc = json.loads(paths["annotations"].read_text())
    assert doc["n_vocs"] >= 50
    assert doc.get("mirror") == HMDB_MIRROR
    assert any(v.get("source") == "hmdbfix_wishart_bulk_xml" for v in doc["vocs"])
    info = src.fuse(ROOT / "data" / "knowledge")
    assert info.get("n_breath_metabolites") == 60
