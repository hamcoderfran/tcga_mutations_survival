"""Mental-health literature panels, claim grades, and de-cloned priors."""

from __future__ import annotations

import json
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


def test_mh_panel_masked_prior_eval_runs():
    from exhalepath.eval.mental_health_lit import evaluate_mental_health_literature

    report = evaluate_mental_health_literature(out_dir=None)
    assert report["n_panel_diseases"] == 3
    by = {c["disease_id"]: c for c in report["cases"] if not c.get("skipped")}
    assert "schizophrenia" in by
    assert "major_depressive_disorder" in by
    assert by["major_depressive_disorder"]["n_quantified_vocs"] >= 3
    # masking should remove panel VOC prior entries
    assert by["major_depressive_disorder"]["n_prior_vocs_removed"] >= 1
    assert by["major_depressive_disorder"]["raw"]["directional_accuracy"] is not None
    assert by["major_depressive_disorder"]["panel_masked_prior"]["directional_accuracy"] is not None
    # non-circular recovery via pathway_bias + disease-local negative voc_effects
    # (not VOC prior re-injection; global ketone/gut TMA maps stay positive)
    assert by["schizophrenia"]["panel_masked_prior"]["directional_accuracy"] >= 0.95
    assert by["major_depressive_disorder"]["panel_masked_prior"]["directional_accuracy"] >= 0.95
    assert report["mean_panel_masked_directional_accuracy"] >= 0.95
    assert report["mean_mechanism_backed_directional_accuracy"] is not None
    assert report["mean_mechanism_backed_directional_accuracy"] >= 0.9
    # MDD SCFAs / butylamine must be mechanism-backed (not near-floor sign luck)
    mdd_rows = {
        r["voc_id"]: r
        for r in by["major_depressive_disorder"]["panel_masked_prior"]["voc_rows"]
    }
    for voc in ("butyric_acid", "acetic_acid", "valeric_acid", "butylamine", "trimethylamine"):
        assert mdd_rows[voc]["agree"] is True
        assert mdd_rows[voc]["mechanism_backed"] is True
        assert abs(float(mdd_rows[voc]["predicted_log2fc"])) >= 0.02
    scz_rows = {r["voc_id"]: r for r in by["schizophrenia"]["panel_masked_prior"]["voc_rows"]}
    assert scz_rows["acetone"]["agree"] is True
    assert scz_rows["acetone"]["mechanism_backed"] is True
    assert abs(float(scz_rows["acetone"]["predicted_log2fc"])) >= 0.02
    for voc in ("isoprene", "carbon_disulfide", "methanol"):
        assert scz_rows[voc]["agree"] is True
        assert scz_rows[voc]["mechanism_backed"] is True
        assert abs(float(scz_rows[voc]["predicted_log2fc"])) >= 0.02
    assert mdd_rows["isoprene"]["mechanism_backed"] is True
    assert mdd_rows["ethanol"]["mechanism_backed"] is True
    assert abs(float(mdd_rows["ethanol"]["predicted_log2fc"])) >= 0.02
    # near-floor list should be empty for panel diseases after pathway wiring
    assert by["schizophrenia"]["panel_masked_prior"].get("near_floor_vocs") == []
    assert by["major_depressive_disorder"]["panel_masked_prior"].get("near_floor_vocs") == []
    # thin conditions documented
    thin_ids = {t["disease_id"] for t in report["thin_evidence_conditions"]}
    assert {"anxiety", "ptsd", "adhd", "autism_spectrum_disorder"} <= thin_ids
    # honesty: raw often perfect when circular; masked is the research metric
    assert "panel_masked_prior" in report["honesty"] or "de-circular" in report["honesty"].lower() or "masked" in report["honesty"]
    assert "mechanism" in report["honesty"].lower()


def test_mh_negative_pathway_effects_exist_for_masked_recovery():
    from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge

    clear_knowledge_cache()
    kb = default_knowledge()
    hypo = kb.pathways["brain_energy_hypometabolism"]
    choline = kb.pathways["choline_TMA_TMAO_axis"]
    scfa = kb.pathways["scfa_metabolism"]
    amine = kb.pathways["amino_acid_decarboxylation"]
    assert hypo["voc_effects"]["acetone"] < 0
    assert choline["voc_effects"]["trimethylamine"] < 0
    assert scfa["voc_effects"]["butyric_acid"] < 0
    assert amine["voc_effects"]["butylamine"] < 0
    assert kb.pathways["mevalonate_flux_suppression"]["voc_effects"]["isoprene"] < 0
    assert kb.pathways["methanol_one_carbon_suppression"]["voc_effects"]["methanol"] < 0
    assert kb.pathways["pyruvate_ethanol_axis"]["voc_effects"]["ethanol"] > 0
    # global ketone / gut fermentation acetone/TMA stay non-negative (T2D/SIBO honesty)
    assert kb.pathways["ketone_body_metabolism"]["voc_effects"]["acetone"] > 0
    assert kb.pathways["gut_microbiome_fermentation"]["voc_effects"]["trimethylamine"] > 0
    assert kb.pathways["mevalonate_cholesterol"]["voc_effects"]["isoprene"] > 0
    assert kb.diseases["schizophrenia"]["pathway_bias"].get("brain_energy_hypometabolism", 1.0) > 1.05
    assert kb.diseases["major_depressive_disorder"]["pathway_bias"].get("choline_TMA_TMAO_axis", 1.0) > 1.05
    assert kb.diseases["major_depressive_disorder"]["pathway_bias"].get("scfa_metabolism", 1.0) > 1.05
    assert kb.diseases["schizophrenia"]["pathway_bias"].get("methionine_transsulfuration", 1.0) > 1.05
    assert kb.diseases["major_depressive_disorder"]["pathway_bias"].get("pyruvate_ethanol_axis", 1.0) > 1.05


def test_negative_mh_pathways_do_not_invert_t2d_or_sibo():
    """Disease-local MH negative effects must not flip T2D acetone or SIBO TMA."""
    from exhalepath.biomarker import ExhaleBiomarkerEngine
    from exhalepath.knowledge.loader import clear_knowledge_cache

    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    t2d = eng.predict(disease="type_2_diabetes", location="liver", top_n=40, explain=False)
    sibo = eng.predict(disease="sibo", location="gut", top_n=40, explain=False)
    t2d_by = {p.voc_id: p.log2_fold_change for p in t2d.result.bundle.predictions}
    sibo_by = {p.voc_id: p.log2_fold_change for p in sibo.result.bundle.predictions}
    assert t2d_by["acetone"] > 0.5
    assert sibo_by["trimethylamine"] > 0.2


def test_bipolar_h2s_dms_not_invented_as_panels():
    panels = {p["disease_id"]: p for p in load_literature_panels()}
    bd = panels["bipolar"]
    measured = set(bd["measured_log2fc"])
    assert "methyl_mercaptan" in measured
    assert "hydrogen_sulfide" not in measured
    assert "dms" not in measured
    note = bd.get("assayed_but_unreported") or {}
    assert "hydrogen_sulfide" in (note.get("vocs") or [])
    assert "dms" in (note.get("vocs") or [])


def test_packaged_mh_readiness_matches_data_tree():
    """Avoid data/ vs src/ honesty drift for mental_health_readiness.json."""
    from exhalepath.config import DATA_DIR, PACKAGE_ROOT

    a = DATA_DIR / "real_breath" / "literature_panels" / "mental_health_readiness.json"
    b = PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "mental_health_readiness.json"
    # When both exist under the same resolve, still assert near-floor cleared
    doc = json.loads(a.read_text() if a.exists() else b.read_text())
    for did in ("schizophrenia", "major_depressive_disorder", "bipolar"):
        assert doc["ready"][did]["de_circularized"].get("near_floor_vocs") == []
        assert doc["ready"][did]["de_circularized"].get("mechanism_backed") == 1.0
    if a.exists() and b.exists() and a.resolve() != b.resolve():
        assert a.read_text() == b.read_text()


def test_gbaoui_mdd_scfa_mean_fixture_log2fc():
    from pathlib import Path
    import math

    from exhalepath.config import DATA_DIR, PACKAGE_ROOT

    paths = [
        DATA_DIR / "real_breath" / "literature_panels" / "gbaoui_mdd_scfa_mean_fixture.json",
        PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "gbaoui_mdd_scfa_mean_fixture.json",
        Path("data/real_breath/literature_panels/gbaoui_mdd_scfa_mean_fixture.json"),
    ]
    path = next(p for p in paths if p.exists())
    doc = json.loads(path.read_text())
    assert "NOT patient-level" in doc["honesty"] or "not patient" in doc["honesty"].lower()
    by = {s["label_int"]: s for s in doc["samples"]}
    mdd, hc = by[1], by[0]
    for voc in ("butyric_acid", "acetic_acid", "valeric_acid"):
        got = math.log2(mdd[voc] / hc[voc])
        assert abs(got - doc["expected_log2fc"][voc]) < 1e-3



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
