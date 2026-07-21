from pathlib import Path

from exhalepath.ingest.public_breath import (
    build_public_breath_benchmark,
    map_compound_to_voc,
)
from exhalepath.knowledge.loader import KnowledgeBase


def test_map_benzaldehyde_and_hexanal():
    kb = KnowledgeBase()
    assert map_compound_to_voc(pubchem_cid=240, iupac_name="benzaldehyde", kb=kb) == "benzaldehyde"
    assert map_compound_to_voc(pubchem_cid=6184, iupac_name="hexanal", kb=kb) == "hexanal"


def test_build_public_breath_benchmark(tmp_path: Path):
    # Uses already-downloaded Figshare files if present; otherwise downloads
    paths = build_public_breath_benchmark(out_dir=tmp_path / "pb")
    assert paths["benchmark"].exists()
    import json

    bench = json.loads(paths["benchmark"].read_text())
    assert len(bench["cases"]) >= 4
    assert any(c["source"].startswith("scientific_data") for c in bench["cases"])
