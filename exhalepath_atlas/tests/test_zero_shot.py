from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache
from exhalepath.knowledge.disease_gene_index import clear_gene_index_cache
from exhalepath.eval.zero_shot_reliability import (
    evaluate_holdout_profiles,
    evaluate_leave_disease_out_priors,
)


def setup_function():
    clear_knowledge_cache()
    clear_gene_index_cache()


def test_bare_cancer_still_unresolved():
    kb = KnowledgeBase()
    d = kb.resolve_disease("cancer")
    assert d.get("_unresolved") is True
    # category cue may enrich mechanism, but must not fuzzy-bind to a specific cancer
    assert str(d["disease_id"]).startswith("custom::")


def test_nonsense_stays_near_healthy():
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    report = eng.predict("completely novel ailment QZZ-42", location="systemic", explain=False)
    assert report.disease_id.startswith("custom::")
    peak = max(abs(p.log2_fold_change) for p in report.result.bundle.predictions)
    assert peak < 0.45
    meta = report.result.bundle.metadata
    assert meta.get("zero_shot_mode") == "near_healthy_fallback"
    assert not (eng.kb.resolve_disease("completely novel ailment QZZ-42").get("voc_log2fc_prior"))


def test_msud_zero_shot_mechanism():
    kb = KnowledgeBase()
    d = kb.resolve_disease("Maple syrup urine disease")
    assert d.get("_unresolved") is True
    assert d.get("_zero_shot_mode") == "mechanism"
    assert d.get("voc_log2fc_prior") == {}
    assert "ketone_body_metabolism" in (d.get("pathway_bias") or {})
    assert "BCKDHA" in (d.get("driver_genes") or [])

    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    report = eng.predict("Maple syrup urine disease", location="systemic", explain=False)
    by = {p.voc_id: p for p in report.result.bundle.predictions}
    assert by["acetone"].log2_fold_change > 0.1


def test_user_genes_move_custom_disease():
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    report = eng.predict(
        "brand new unnamed ailment",
        location="liver",
        genes=["HMGCS2", "CPT1A"],
        explain=False,
    )
    by = {p.voc_id: p for p in report.result.bundle.predictions}
    assert by["acetone"].log2_fold_change > 0.2


def test_ontology_nn_no_voc_prior_copy():
    kb = KnowledgeBase()
    # Cholestatic phrase should pick liver metabolic mechanisms without VOC priors
    d = kb.resolve_disease("primary biliary cholangitis")
    assert d.get("_unresolved") is True
    assert d.get("voc_log2fc_prior") == {}
    assert d.get("default_site") in {"liver", "Liver"} or "liver" in str(d.get("default_site")).lower()
    assert float(d.get("_mechanism_confidence") or 0) >= 0.55


def test_holdout_and_ldo_smoke():
    hold = evaluate_holdout_profiles()
    assert hold["n_profiles"] >= 10
    assert hold["pass_rate"] >= 0.75
    ldo = evaluate_leave_disease_out_priors(max_diseases=25)
    assert not ldo.get("skipped")
    assert ldo["mean_directional_accuracy"] >= 0.5
