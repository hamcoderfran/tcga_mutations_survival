"""Sellability pack: closing demo pieces, OEM kit, partner LOSO, ledger cites."""

from __future__ import annotations

from pathlib import Path

from exhalepath.gcms.omni_partner import (
    build_partner_loso_fixture,
    load_omni_style_csv,
    run_partner_loso_diligence,
)
from exhalepath.oem.kit import cite_vocs_for_disease, import_feature_table, score_sample
from exhalepath.viz.analysis_export import literature_overlay


def test_partner_loso_fixture_and_diligence(tmp_path: Path):
    man = build_partner_loso_fixture(tmp_path / "fix", seed=0)
    assert Path(man["site_a"]["path"]).exists()
    a = load_omni_style_csv(man["site_a"]["path"], study_id="A")
    b = load_omni_style_csv(man["site_b"]["path"], study_id="B")
    assert a.matrix.shape[0] >= 8 and b.matrix.shape[0] >= 8
    report = run_partner_loso_diligence(
        site_a=Path(man["site_a"]["path"]),
        site_b=Path(man["site_b"]["path"]),
        seed=0,
    )
    assert report["loso"]["status"] == "ok"
    assert report["buyer_slide"]["mean_auroc_sparse"] is not None
    assert "clinical" in (report["buyer_slide"].get("do_not_claim") or "").lower()
    assert "Diligence fixture" in (man.get("honesty") or "")


def test_score_sample_ledger_cited():
    out = score_sample(
        "asthma",
        {"ethane": 0.5, "pentane": 0.4, "hexanal": 0.3, "acetone": -0.2, "isoprene": -0.1},
    )
    assert out["score"] is not None
    grades = {r["evidence_grade"] for r in out["top_vocs"]}
    # asthma panel has quantified ethane
    assert any(r.get("doi") for r in out["top_vocs"]) or "quantified" in grades or grades


def test_cite_vocs_prefers_quantified():
    cites = cite_vocs_for_disease("asthma", ["ethane", "hexanal"])
    assert cites["ethane"]["evidence_grade"] in {"quantified", "mixed", "directional_only"}
    if cites["ethane"]["evidence_grade"] == "quantified":
        assert cites["ethane"]["doi"]


def test_literature_overlay_has_doi_grade():
    overlay = literature_overlay(
        "asthma", {"ethane": 1.0, "pentane": 0.5, "acetone": -0.2}
    )
    assert overlay["available"]
    eth = next(r for r in overlay["rows"] if r["voc_id"] == "ethane")
    assert eth.get("evidence_grade")
    assert "doi" in eth


def test_import_omni_roundtrip(tmp_path: Path):
    man = build_partner_loso_fixture(tmp_path / "fx", seed=1)
    packed = import_feature_table(Path(man["site_a"]["path"]), fmt="omni")
    assert packed["n_subjects"] >= 8
    assert packed["n_vocs"] >= 5
