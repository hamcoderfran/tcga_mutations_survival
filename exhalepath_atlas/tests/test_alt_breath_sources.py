"""Open alternatives to Owlstone / live HBDB — access audit + Zenodo harvest."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.datasources.ds15_alt_breath_sources import AltBreathSources

ROOT = Path(__file__).resolve().parents[1]


def test_alt_breath_harvest_online_or_offline():
    src = AltBreathSources(ROOT)
    # Prefer online once so artifacts exist; fall back offline in CI without net.
    try:
        paths = src.harvest(offline=False)
    except Exception:
        paths = src.harvest(offline=True)
    assert paths["access_matrix"].exists()
    matrix = json.loads(paths["access_matrix"].read_text())
    keys = {s["key"] for s in matrix["sources"]}
    assert "clinical_breathomics_figshare" in keys
    assert "owlstone_voc_atlas" in keys
    assert "breathomix_breathbase" in keys
    assert "hmdb" in keys
    assert "physionet" in keys
    owl = next(s for s in matrix["sources"] if s["key"] == "owlstone_voc_atlas")
    assert owl["in_atlas"] is False
    assert "license" in owl["access"].lower() or "Registration" in owl["access"]

    info = src.fuse(ROOT / "data" / "knowledge")
    fused = json.loads(Path(info["path"]).read_text())
    assert len(fused["access_matrix"]["sources"]) >= 7
    # Clinical breathomics should already be present in this repo
    cbd = next(
        s for s in matrix["sources"] if s["key"] == "clinical_breathomics_figshare"
    )
    assert cbd["in_atlas"] is True


def test_alt_breath_offline_reuses_artifacts():
    src = AltBreathSources(ROOT)
    if not (src.out_dir / "ACCESS_MATRIX.json").exists():
        src.harvest(offline=False)
    paths = src.harvest(offline=True)
    assert paths["manifest"].exists()
    man = json.loads(paths["manifest"].read_text())
    assert man.get("offline") is True
