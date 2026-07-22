"""Regression tests for vision-hardening (tissues, fuzzy match, alias collisions)."""

from __future__ import annotations

from exhalepath.body.tissues import resolve_location, default_body_map
from exhalepath.eval.vision import _fuzzy_hit
from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache


def test_head_neck_does_not_resolve_to_head():
    default_body_map.cache_clear()
    loc = resolve_location("head_neck")
    assert loc["matched"] is True
    assert loc["tissue_id"] == "head_neck"


def test_new_anatomic_sites_resolve():
    default_body_map.cache_clear()
    for q, expect in [
        ("pharynx", "pharynx"),
        ("joint", "joint"),
        ("pelvis", "pelvis"),
        ("thyroid", "thyroid"),
        ("oral", "oral"),
    ]:
        loc = resolve_location(q)
        assert loc["matched"] is True, q
        assert loc["tissue_id"] == expect, (q, loc)


def test_fuzzy_cell_no_false_positive_on_cell_token():
    # Shared filler "cell" must not match unrelated atlas states
    assert not _fuzzy_hit(
        "endothelial cell",
        ["Oxidative-stress / lipid-peroxidizing cell", "oxidative_stress_cell"],
    )
    assert _fuzzy_hit(
        "Oxidative endothelial cell",
        ["Oxidative endothelial cell", "endothelial_oxidative"],
    )
    assert _fuzzy_hit(
        "Ketogenic hepatocyte",
        ["Ketogenic hepatocyte", "hepatocyte_ketogenic"],
    )


def test_influenza_no_longer_steals_pneumonia_alias():
    clear_knowledge_cache()
    kb = KnowledgeBase()
    flu = kb.diseases["influenza"]
    assert "pneumonia" not in [a.lower() for a in (flu.get("aliases") or [])]
    assert kb.resolve_disease("pneumonia")["disease_id"] == "pneumonia_bacterial"
    assert kb.resolve_disease("influenza")["disease_id"] == "influenza"


def test_heart_failure_alias_not_on_heart_disease():
    clear_knowledge_cache()
    kb = KnowledgeBase()
    hd = kb.diseases["heart_disease"]
    assert "heart failure" not in [a.lower() for a in (hd.get("aliases") or [])]
    assert kb.resolve_disease("heart failure")["disease_id"] == "heart_failure"


def test_new_cell_states_present():
    clear_knowledge_cache()
    kb = KnowledgeBase()
    ids = {c["state_id"] for c in kb.cell_states}
    for sid in [
        "airway_epithelium_stressed",
        "macrophage_activated",
        "erythrocyte_oxidative",
        "cardiomyocyte_stressed",
        "endothelial_oxidative",
        "synovial_inflammatory",
        "thyroid_follicular_hypermetabolic",
        "colonocyte_inflamed",
    ]:
        assert sid in ids, sid
