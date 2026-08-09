"""Malaria lit remap + nested sparse + fixed-sens + external catalog."""

from __future__ import annotations

from pathlib import Path

from exhalepath.eval.malaria_diagnostic import evaluate_malaria_diagnostic
from exhalepath.gcms import load_mw_patient_matrix, map_malaria_lit_voc
from exhalepath.gcms.metrics import fixed_sensitivity_metrics
from exhalepath.gcms.nested_sparse import nested_sparse_logistic
from exhalepath.gcms.locked_split import lock_split


def test_malaria_lit_remap_names():
    assert map_malaria_lit_voc("?-Pinene(93@8.46)") == "alpha_pinene"
    assert map_malaria_lit_voc("3-carene(41@10.91)") == "delta_3_carene"
    assert map_malaria_lit_voc("Cyclohexanone (43@7.34)") == "cyclohexanone"
    assert map_malaria_lit_voc("Tridecane(43@20.26)") == "tridecane"
    assert map_malaria_lit_voc("Allyl methyl sulfide") == "allyl_methyl_sulfide"


def test_st000883_malaria_lit_expands_features():
    atlas = load_mw_patient_matrix("ST000883", feature_map="atlas")
    lit = load_mw_patient_matrix("ST000883", feature_map="malaria_lit")
    assert lit.matrix.shape[1] > atlas.matrix.shape[1]
    for v in ("alpha_pinene", "delta_3_carene", "cyclohexanone", "tridecane"):
        assert v in lit.matrix.columns


def test_nested_sparse_and_fixed_sens():
    m = load_mw_patient_matrix("ST000883", feature_map="malaria_lit").log1p()
    man = lock_split(m, strategy="stratified_kfold", n_splits=5, seed=0)
    sparse = nested_sparse_logistic(m, man, max_features=6, selector="kbest", seed=0)
    assert sparse["mean_test_auroc"] is not None
    assert 0.0 <= sparse["mean_test_auroc"] <= 1.0
    assert sparse["consensus_features"]
    # fixed sens on labels vs random-ish scores still returns structure
    y = m.labels.tolist()
    s = m.matrix.mean(axis=1).tolist()
    fs = fixed_sensitivity_metrics(y, s, target_sensitivities=[0.9], n_boot=50, seed=0)
    assert fs["points"]


def test_evaluate_malaria_diagnostic(tmp_path: Path):
    report = evaluate_malaria_diagnostic(
        out_dir=tmp_path / "mal",
        feature_map="malaria_lit",
        n_splits=5,
        seed=1,
        max_features=6,
    )
    t = report["transferable_signature"]
    sk = report["fit_on_cohort_sparse_kbest"]
    assert t["nested_auroc"] is not None
    assert sk["mean_test_auroc"] is not None
    assert (tmp_path / "mal" / "MALARIA_DIAGNOSTIC_UPGRADE.md").exists()
    assert (tmp_path / "mal" / "fixed_sensitivity.json").exists()
    assert (tmp_path / "mal" / "learning_curve.json").exists()
    assert report["comparison_to_baseline_atlas_map"]["newly_mapped_vocs"]
    # sparse ceiling should be at least competitive with transferable on this tiny n
    assert sk["mean_test_auroc"] >= 0.5
