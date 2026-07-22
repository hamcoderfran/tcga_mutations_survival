"""Tests for the 50-disease vision evaluation suite."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.eval.vision import evaluate_vision, _fuzzy_hit, _load_profiles


def test_vision_profiles_exist_and_count_50():
    doc = _load_profiles()
    assert doc["n_profiles"] == 50
    assert len(doc["profiles"]) == 50
    grades = {p["grade"] for p in doc["profiles"]}
    assert grades <= {"A", "B", "C", "D"}
    assert "A" in {p["grade"] for p in doc["profiles"]}
    ids = [p["id"] for p in doc["profiles"]]
    assert len(ids) == len(set(ids))


def test_fuzzy_cell_match():
    assert _fuzzy_hit("Ketogenic hepatocyte", ["Ketogenic hepatocyte", "Lipolytic adipocyte"])
    assert _fuzzy_hit("airway epithelium", ["Airway epithelium", "eosinophil"])
    assert not _fuzzy_hit("plasmodium", ["Ketogenic hepatocyte"])


def test_evaluate_vision_smoke(tmp_path: Path):
    report = evaluate_vision(out_dir=tmp_path, limit=3, top_k=15, mode="hybrid")
    assert report["overall"]["n_profiles"] == 3
    assert report["overall"]["vision_fidelity_pct"] is not None
    assert (tmp_path / "vision_eval.json").exists()
    assert (tmp_path / "VISION_EVAL.md").exists()
    raw = json.loads((tmp_path / "vision_eval.json").read_text())
    assert len(raw["cases"]) == 3
    for c in raw["cases"]:
        assert "scores" in c
        assert c["scores"]["composite"] is not None
