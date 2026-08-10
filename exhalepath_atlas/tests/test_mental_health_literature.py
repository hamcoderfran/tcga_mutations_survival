"""Mental-health literature panels, claim grades, and de-cloned priors."""

from __future__ import annotations

import math

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.data.literature_panels import (
    literature_panel_for_disease,
    load_literature_panels,
    load_thin_evidence_notes,
)
from exhalepath.gcms.claim_ledger import build_claim_ledger
from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge
from exhalepath.viz.analysis_export import literature_overlay


def test_mh_literature_panels_loaded_with_dois():
    panels = {p["disease_id"]: p for p in load_literature_panels()}
    assert "schizophrenia" in panels
    assert "major_depressive_disorder" in panels
    assert "bipolar" in panels
    scz = panels["schizophrenia"]
    assert any(r.get("doi") == "10.1080/15622975.2022.2040052" for r in scz["refs"])
    assert scz["measured_log2fc"]["trimethylamine"] < 0
    assert scz["measured_log2fc"]["butyric_acid"] < 0
    assert scz["measured_log2fc"]["pentane"] > 0
    mdd = panels["major_depressive_disorder"]
    assert mdd["voc_evidence"]["butyric_acid"]["evidence"] == "quantified"
    assert abs(mdd["measured_log2fc"]["butyric_acid"] - math.log2(116 / 169)) < 1e-3
    bd = panels["bipolar"]
    expected = math.log2(18.62 / 9.45)
    assert abs(bd["measured_log2fc"]["methyl_mercaptan"] - expected) < 1e-3
    assert bd["voc_evidence"]["methyl_mercaptan"]["evidence"] == "quantified"


def test_butyric_acid_in_catalog_and_predictions():
    clear_knowledge_cache()
    kb = default_knowledge()
    for vid in ("butyric_acid", "acetic_acid", "valeric_acid", "butylamine"):
        assert vid in kb.vocs
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    r = eng.predict(disease="schizophrenia", location="brain", top_n=60, explain=False)
    by = {p.voc_id: p for p in r.result.bundle.predictions}
    assert by["butyric_acid"].fold_change < 1.0
    assert by["butylamine"].fold_change < 1.0


def test_mdd_quantified_scfa_panel():
    panels = {p["disease_id"]: p for p in load_literature_panels()}
    mdd = panels["major_depressive_disorder"]
    assert mdd["voc_evidence"]["acetic_acid"]["evidence"] == "quantified"
    assert mdd["voc_evidence"]["valeric_acid"]["evidence"] == "quantified"
    assert abs(mdd["measured_log2fc"]["acetic_acid"] - math.log2(124 / 146)) < 1e-3
    assert abs(mdd["measured_log2fc"]["valeric_acid"] - math.log2(4 / 8)) < 1e-3


def test_overlay_reports_circularity_outside_prior():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    r = eng.predict(disease="major_depressive_disorder", location="brain", top_n=60, explain=False)
    pred = {p.voc_id: float(p.log2_fold_change) for p in r.result.bundle.predictions}
    # empty prior → outside-prior accuracy equals raw
    ov = literature_overlay("major_depressive_disorder", pred, prior_log2fc={})
    assert ov["available"]
    assert ov["n_compared_outside_prior"] == ov["n_compared"]
    # full prior overlap → outside-prior null
    prior = dict(
        (eng.kb.diseases.get("major_depressive_disorder") or {}).get("voc_log2fc_prior") or {}
    )
    ov2 = literature_overlay("major_depressive_disorder", pred, prior_log2fc=prior)
    assert "directional_accuracy_outside_prior" in ov2
    assert ov2["n_panel_vocs_also_in_prior"] >= 1
    assert any(row.get("in_atlas_prior") for row in ov2["rows"])


def test_thin_mh_conditions_documented_not_faked():
    notes = load_thin_evidence_notes()
    for did in ("anxiety", "ptsd", "adhd", "autism_spectrum_disorder"):
        assert did in notes
        assert notes[did].get("refs")
    # no fake measured panels for thin conditions
    panels = {p["disease_id"] for p in load_literature_panels()}
    assert "anxiety" not in panels
    assert "ptsd" not in panels
    assert "adhd" not in panels
    assert "autism_spectrum_disorder" not in panels


def test_mh_claim_ledger_not_atlas_prior_only():
    ledger = build_claim_ledger(
        disease_ids=["schizophrenia", "major_depressive_disorder", "bipolar"],
        include_atlas_priors=True,
    )
    by = {}
    for c in ledger["claims"]:
        by.setdefault(c["disease_id"], set()).add(c["evidence_grade"])
    for did in ("schizophrenia", "major_depressive_disorder", "bipolar"):
        grades = by[did]
        assert "atlas_prior" in grades  # still have prior-only VOCs
        assert grades & {"directional_only", "quantified", "mixed"}
    # bipolar has quantified CH3SH with DOI
    quant = [
        c
        for c in ledger["claims"]
        if c["disease_id"] == "bipolar" and c["voc_id"] == "methyl_mercaptan"
    ]
    assert quant and quant[0]["evidence_grade"] == "quantified"
    assert quant[0]["doi"] == "10.3390/jcm14062025"


def test_mh_priors_no_longer_template_clones():
    clear_knowledge_cache()
    kb = default_knowledge()
    ids = ["bipolar", "anxiety", "ptsd", "adhd"]
    vecs = {d: dict((kb.diseases[d] or {}).get("voc_log2fc_prior") or {}) for d in ids}
    # exact clone check
    assert vecs["bipolar"] != vecs["anxiety"]
    assert vecs["anxiety"] != vecs["ptsd"]
    assert vecs["ptsd"] != vecs["adhd"]
    assert kb.diseases["bipolar"].get("atlas_source") == "literature_oralchroma_ch3sh"
    assert "thin" in (kb.diseases["anxiety"].get("atlas_source") or "")
    assert "not_exhaled" in (kb.diseases["autism_spectrum_disorder"].get("atlas_source") or "")


def test_schizophrenia_and_mdd_literature_overlay_available():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    for did, loc in (("schizophrenia", "brain"), ("major_depressive_disorder", "brain"), ("bipolar", "brain")):
        r = eng.predict(disease=did, location=loc, top_n=50, explain=False)
        pred = {p.voc_id: float(p.log2_fold_change) for p in r.result.bundle.predictions}
        overlay = literature_overlay(did, pred)
        assert overlay["available"] is True
        assert overlay["n_compared"] >= 2
        assert any(ref.get("doi") for ref in overlay["refs"])


def test_panel_for_disease_helper():
    p = literature_panel_for_disease("schizophrenia")
    assert p and p["disease_id"] == "schizophrenia"
