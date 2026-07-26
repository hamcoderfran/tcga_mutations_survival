from exhalepath.eval.zero_shot_hard import evaluate_zero_shot_hard


def test_hard_break_suite_majority_pass():
    report = evaluate_zero_shot_hard()
    assert report["n_profiles"] >= 30
    assert report["pass_rate"] >= 0.75
    assert report["pass"] is True
