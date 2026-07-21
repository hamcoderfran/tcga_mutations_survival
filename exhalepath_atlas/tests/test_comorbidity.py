"""Tests for comorbidity fusion and clinical eval wiring."""

from pathlib import Path

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.knowledge.comorbidity import merge_disease_with_comorbidities
from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache


ROOT = Path(__file__).resolve().parents[1]


def test_merge_comorbidities_boosts_voc_priors():
    kb = KnowledgeBase(knowledge_dir=ROOT / "data" / "knowledge")
    primary = kb.resolve_disease("depression")
    obesity = kb.resolve_disease("obesity")
    fused = merge_disease_with_comorbidities(primary, [obesity], weight=0.65)
    assert "obesity" in fused["_comorbid_ids"]
    # Obesity acetone prior should increase fused acetone prior
    assert fused["voc_log2fc_prior"].get("acetone", 0) > primary["voc_log2fc_prior"].get(
        "acetone", 0
    )
    assert fused["pathway_bias"].get("ketone_body_metabolism", 1.0) > 1.0


def test_biomarker_comorbidity_changes_ranking():
    clear_knowledge_cache()
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    base = engine.predict(
        "depression",
        location="brain",
        age_years=24,
        sex="male",
        top_n=20,
        explain=False,
    )
    comorbid = engine.predict(
        "depression",
        location="brain",
        age_years=24,
        sex="male",
        comorbidities=["obesity"],
        top_n=20,
        explain=False,
    )
    assert comorbid.result.bundle.metadata.get("comorbidity_ids") == ["obesity"]
    base_acetone = next(p for p in base.result.bundle.predictions if p.voc_id == "acetone")
    com_acetone = next(
        p for p in comorbid.result.bundle.predictions if p.voc_id == "acetone"
    )
    # Obesity comorbidity should increase acetone Δppb vs depression alone
    assert com_acetone.delta_ppb >= base_acetone.delta_ppb - 1e-6


def test_clinical_benchmark_file_or_harvestable():
    # Prefer committed knowledge file; otherwise harvest builds it
    path = ROOT / "data" / "knowledge" / "comorbidity_clinical_benchmarks.json"
    if not path.exists():
        from exhalepath.ingest.clinical_comorbidity import (
            build_comorbidity_clinical_benchmarks,
        )

        path = build_comorbidity_clinical_benchmarks()
    import json

    doc = json.loads(path.read_text())
    assert len(doc["cases"]) >= 4
    ids = {c["case_id"] for c in doc["cases"]}
    assert "clinical_mdd_obesity_comorbid" in ids
    assert "clinical_sz_heart_comorbid" in ids
