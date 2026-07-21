"""Creutzfeldt–Jakob disease atlas registration + clinical profile smoke tests."""

from __future__ import annotations

from exhalepath.knowledge.loader import KnowledgeBase
from exhalepath.nl.rules import parse_rules


def test_cjd_resolves():
    kb = KnowledgeBase()
    d = kb.resolve_disease("CJD")
    assert d["disease_id"] == "creutzfeldt_jakob"
    assert d["default_site"] == "brain"
    assert d["category"] == "neurological"
    assert "lipid_peroxidation" in d["pathway_bias"]
    assert d["voc_log2fc_prior"].get("hexanal", 0) > 0


def test_cjd_nl_parse():
    slots = parse_rules("62yo woman with rapidly progressive Creutzfeldt-Jakob disease")
    assert slots.disease
    low = slots.disease.lower()
    assert "jakob" in low or "creutzfeldt" in low
    assert slots.location == "brain"
    assert slots.age_years == 62
    assert slots.sex == "female"


def test_cjd_profile_tempo():
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    from scripts.cjd_clinical_profile import run_profile

    profile = run_profile(out_dir=root / "runs/cjd_clinical_profile")
    tempo = profile["summary"]["tempo"]
    assert tempo["incubating_verdict"] == "subtle_or_absent"
    assert tempo["early_clinical_verdict"] in {"meaningful_hint", "obvious"}
    assert profile["summary"]["first_meaningful_hint_phase"] in {
        "prodromal",
        "early_clinical",
        "fulminant",
    }
    m = tempo["mean_abs_log2fc"]
    assert m["fulminant"] >= m["incubating"]
    assert m["prodromal"] > m["incubating"]
