"""Regression tests for adversarial stress-suite hardenings."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.body.tissues import default_body_map, resolve_location
from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge
from exhalepath.schemas import DiseaseQuery


@pytest.fixture(autouse=True)
def _clear_caches():
    clear_knowledge_cache()
    default_body_map.cache_clear()
    yield
    clear_knowledge_cache()
    default_body_map.cache_clear()


def test_punctuation_location_unmatched():
    loc = resolve_location("???")
    assert loc["matched"] is False
    assert "custom" in loc["tissue_id"]


def test_gene_symbols_not_diseases():
    kb = default_knowledge()
    for gene in ("BRCA1", "KRAS", "TP53", "EGFR"):
        d = kb.resolve_disease(gene)
        assert d.get("_unresolved") is True, gene
    # Exact TCGA project alias still maps
    assert kb.resolve_disease("BRCA")["disease_id"] == "breast_invasive_carcinoma"


def test_comorbidity_weight_clamped_not_crash():
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    r = eng.predict(
        disease="depression",
        location="brain",
        comorbidities=["obesity"],
        comorbidity_weight=3.5,
        top_n=5,
        explain=False,
    )
    assert r.disease_id
    r2 = eng.predict(
        disease="depression",
        location="brain",
        comorbidities=["obesity"],
        comorbidity_weight=-1.0,
        top_n=5,
        explain=False,
    )
    assert r2.disease_id
    with pytest.raises(ValidationError):
        DiseaseQuery(disease="x", comorbidity_weight=3.5)


def test_schizophrenia_suppresses_magdeburg_panel():
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    r = eng.predict(disease="schizophrenia", location="brain", top_n=40, explain=False)
    by = {p.voc_id: p for p in r.result.bundle.predictions}
    assert by["acetone"].fold_change < 1.0
    assert by["isoprene"].fold_change < 1.0
    assert by["pentane"].fold_change > 1.0
    assert by["ethane"].fold_change > 1.0


def test_autism_mondo_not_heart_failure():
    kb = default_knowledge()
    autism = kb.resolve_disease("autism")
    hf = kb.resolve_disease("heart failure")
    assert autism["disease_id"] == "autism_spectrum_disorder"
    assert hf["disease_id"] == "heart_failure"
    assert autism.get("mondo_id") != hf.get("mondo_id")
