"""Tests for the five irreplaceable VOC OS diligence pillars."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.eval.leaderboard import (
    build_locked_truth_leaderboard,
    write_locked_truth_leaderboard,
)
from exhalepath.eval.panel_decision import decide_panel, write_panel_decision
from exhalepath.eval.partner_diligence import run_partner_diligence_pipeline
from exhalepath.gcms.claim_ledger import (
    build_claim_ledger,
    export_claim_ledger,
    filter_claim_ledger,
)
from exhalepath.knowledge.voc_alias_registry import (
    audit_column_names,
    build_voc_alias_registry,
    map_name_to_voc,
    write_voc_alias_registry,
)


def test_locked_truth_leaderboard_has_nested_rows(tmp_path: Path):
    report = write_locked_truth_leaderboard(out_dir=tmp_path)
    assert report["schema_version"] == "LockedTruthLeaderboard-1.0"
    assert report["n_rows"] >= 1
    assert (tmp_path / "LOCKED_TRUTH_LEADERBOARD.json").exists()
    assert (tmp_path / "LOCKED_TRUTH_LEADERBOARD.md").exists()
    nested_rows = [r for r in report["rows"] if r.get("nested_auroc") is not None]
    assert nested_rows
    assert report["clinical_claim"] is False if "clinical_claim" in report else True
    for r in report["rows"]:
        assert r.get("clinical_claim") is False


def test_partner_diligence_pipeline(tmp_path: Path):
    report = run_partner_diligence_pipeline(out_dir=tmp_path / "pd")
    assert report["n_subjects"] >= 5
    assert report["split_sha256"]
    assert (tmp_path / "pd" / "PARTNER_DILIGENCE.json").exists()
    assert report.get("alias_audit") is not None
    assert report["clinical_claim"] is False


def test_panel_decision_enum(tmp_path: Path):
    report = write_panel_decision(disease_id="malaria", out_dir=tmp_path)
    assert report["decision"] in {"proceed_research", "hold", "stop"}
    assert report["clinical_claim"] is False
    assert "NOT clinical clearance" in report["honesty"]
    assert (tmp_path / "PANEL_DECISION_malaria.md").exists()

    thin = decide_panel(disease_id="anxiety_disorder")
    # thin MH conditions may hold/stop depending on ledger contents
    assert thin["decision"] in {"proceed_research", "hold", "stop"}


def test_claim_ledger_quantified_only_filter(tmp_path: Path):
    full = build_claim_ledger(disease_ids=["major_depressive_disorder"])
    assert full["n_claims"] >= 1
    filtered = filter_claim_ledger(full, quantified_only=True)
    assert filtered["n_claims"] >= 1
    assert all(c["evidence_grade"] == "quantified" for c in filtered["claims"])
    path = export_claim_ledger(
        tmp_path / "CLAIM_LEDGER.json",
        disease_ids=["major_depressive_disorder"],
        quantified_only=True,
    )
    doc = json.loads(path.read_text())
    assert doc["filter"]["quantified_only"] is True
    assert path.with_suffix(".md").exists()


def test_voc_alias_registry_maps_magdeburg_and_aliases(tmp_path: Path):
    reg = write_voc_alias_registry(out_dir=tmp_path)
    assert reg["n_aliases"] >= 20
    assert map_name_to_voc("acetone", reg) == "acetone"
    assert map_name_to_voc("m/z 60", reg) == "trimethylamine"
    assert map_name_to_voc("n-butylamine", reg) == "butylamine"
    audit = audit_column_names(
        ["sample_id", "acetone", "isoprene", "totally_fake_voc_xyz", "label"],
        registry=reg,
    )
    assert audit["n_mapped"] >= 2
    assert audit["n_unmapped"] >= 1
    assert 0.0 < float(audit["unmapped_rate"]) < 1.0
    assert (tmp_path / "VOC_ALIAS_REGISTRY.json").exists()
    assert (tmp_path / "CONFOUNDER_REGISTRY.json").exists()
    conf = json.loads((tmp_path / "CONFOUNDER_REGISTRY.json").read_text())
    assert "smoking" in conf["strata"]


def test_build_leaderboard_idempotent_shape():
    a = build_locked_truth_leaderboard()
    b = build_voc_alias_registry()
    assert a["n_rows"] == build_locked_truth_leaderboard()["n_rows"]
    assert b["n_aliases"] == build_voc_alias_registry()["n_aliases"]
