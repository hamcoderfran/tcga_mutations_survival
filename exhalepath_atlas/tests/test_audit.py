from exhalepath.eval.audit import (
    evaluate_literature_benchmarks,
    evaluate_zero_shot,
    run_full_audit,
    stress_test_pipeline,
)
from exhalepath.model.predict import ExhalePathPredictor


def test_literature_directional_accuracy_high():
    pred = ExhalePathPredictor(use_opentargets=False)
    lit = evaluate_literature_benchmarks(pred)
    assert lit["directional_accuracy"] >= 0.85
    assert lit["case_pass_rate"] >= 0.8


def test_zero_shot_metabolic_and_novel():
    pred = ExhalePathPredictor(use_opentargets=False)
    zs = evaluate_zero_shot(pred)
    assert zs["pass_rate"] >= 0.75
    assert zs["novel_vs_t2d_acetone_delta_ppb"] > 100


def test_stress_suite_passes():
    pred = ExhalePathPredictor(use_opentargets=False)
    stress = stress_test_pipeline(pred)
    assert stress["pass"], stress["failures"]


def test_full_audit_passes(tmp_path):
    report = run_full_audit(out_path=tmp_path / "audit.json")
    assert (tmp_path / "audit.json").exists()
    assert report.passed, report.to_dict()
