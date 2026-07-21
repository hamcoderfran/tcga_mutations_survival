from pathlib import Path

from exhalepath.config import DATA_DIR
from exhalepath.physio.census_fractions import census_cell_state_fractions_for_disease


def test_census_fractions_available_for_alzheimer_if_harvested():
    frac_path = DATA_DIR / "census" / "disease_exhalepath_state_fractions.csv"
    if not frac_path.exists():
        return  # harvest not present in CI clone without data
    fracs = census_cell_state_fractions_for_disease("alzheimer_disease")
    assert fracs
    assert abs(sum(fracs.values()) - 1.0) < 1e-5


def test_census_manifest_records_scale():
    manifest = DATA_DIR / "census" / "census_manifest.json"
    if not manifest.exists():
        return
    import json

    m = json.loads(manifest.read_text())
    assert m["n_primary_cells"] > 10_000_000
    assert m["n_normal_cells"] > 1_000_000
    assert m["n_tissues"] >= 20
