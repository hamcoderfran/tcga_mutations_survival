from exhalepath.eval.audit import evaluate_calibrator_holdout
from exhalepath.eval.zero_shot_hard import evaluate_zero_shot_hard


def test_calibrator_holdout_not_skipped_on_real_corpus():
    cal = evaluate_calibrator_holdout()
    assert not cal.get("skipped"), cal
    assert cal.get("pass") is True
    assert cal.get("n", 0) >= 50


def test_expanded_hard_break_suite():
    report = evaluate_zero_shot_hard()
    assert report["n_profiles"] >= 50
    assert report["pass_rate"] >= 0.85
