"""Tests for the 16-model Great Disease Prediction Stack."""

from __future__ import annotations

from pathlib import Path

from exhalepath.great_stack import GreatDiseaseStack, run_great_stack
from exhalepath.great_stack.data import great_stack_dir, model_registry
from exhalepath.great_stack.report import save_stack_report
from exhalepath.great_stack.types import StackQuery


def test_registry_has_16_models():
    reg = model_registry()
    assert len(reg.get("models") or []) >= 16
    assert (great_stack_dir() / "primekg_lite_graph.json").exists()
    assert (great_stack_dir() / "opera_adme_panel.json").exists()


def test_stack_runs_and_fuses(tmp_path: Path):
    result = run_great_stack(
        "depression",
        location="brain",
        comorbidities=["obesity"],
        age=24,
        sex="male",
        top_n=12,
    )
    assert result.disease_id
    assert len(result.model_outputs) >= 16
    ok = [m for m in result.model_outputs if m.status == "ok"]
    assert len(ok) >= 10
    assert result.fused_vocs
    assert result.summary["n_models_ok"] >= 10
    paths = save_stack_report(result, tmp_path / "stack")
    assert paths["html"].is_file()
    assert paths["dashboard"].is_file()
    assert paths["json"].is_file()


def test_stack_novel_disease():
    stack = GreatDiseaseStack()
    result = stack.predict(
        StackQuery(
            disease="Maple syrup urine disease",
            location="systemic",
            genes=["BCKDHA", "BCKDHB"],
            top_n=10,
        )
    )
    assert result.fused_vocs
    ids = {m.model_id for m in result.model_outputs}
    assert "humangem_flux" in ids
    assert "primekg_graph" in ids
    assert "agora_microbiome" in ids
