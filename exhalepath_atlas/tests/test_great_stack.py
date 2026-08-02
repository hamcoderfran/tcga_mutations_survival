"""Tests for the 20-model Great Disease Prediction Stack."""

from __future__ import annotations

from pathlib import Path

from exhalepath.great_stack import GreatDiseaseStack, run_great_stack
from exhalepath.great_stack.data import great_stack_dir, model_registry
from exhalepath.great_stack.fusion import detect_zero_shot
from exhalepath.great_stack.report import save_stack_report
from exhalepath.great_stack.types import StackQuery


def test_registry_has_20_models():
    reg = model_registry()
    assert len(reg.get("models") or []) >= 20
    assert (great_stack_dir() / "primekg_lite_graph.json").exists()
    assert (great_stack_dir() / "opera_adme_panel.json").exists()
    assert (great_stack_dir() / "fusion_weights_calibrated.json").exists()


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
    assert len(result.model_outputs) >= 20
    ok = [m for m in result.model_outputs if m.status == "ok"]
    assert len(ok) >= 12
    assert result.fused_vocs
    assert result.summary["n_models_ok"] >= 12
    assert result.summary.get("zero_shot_mode") is False
    assert result.fused_vocs[0].epistemic_std >= 0.0
    assert "anti_dilution" in (result.summary.get("fusion_method") or "")
    paths = save_stack_report(result, tmp_path / "stack")
    assert paths["html"].is_file()
    assert paths["dashboard"].is_file()
    assert paths["json"].is_file()


def test_stack_novel_disease_zero_shot():
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
    assert detect_zero_shot(result.model_outputs, result.query) is True
    ids = {m.model_id: m for m in result.model_outputs}
    assert "humangem_flux" in ids
    assert "primekg_graph" in ids
    assert "agora_microbiome" in ids
    assert "zero_shot_evidence" in ids
    assert "phenotype_mondo" in ids
    assert "meta_ensemble" in ids
    assert len(ids["zero_shot_mechanism"].voc_signals) >= 3
    assert len(ids["zero_shot_evidence"].voc_signals) >= 3
    pred = {v.voc_id: v.fused_log2fc for v in result.fused_vocs}
    assert pred.get("acetone", 0) > 0 or pred.get("2_butanone", 0) > 0
