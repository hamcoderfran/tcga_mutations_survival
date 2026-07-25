"""Novel diseases with empty VOC priors stay near-healthy."""

from pathlib import Path

from exhalepath.eval.novel_diseases import run_novel_disease_eval


def test_novel_diseases_near_healthy(tmp_path: Path):
    report = run_novel_disease_eval(out_dir=tmp_path / "novel")
    assert report["n_novel"] >= 8
    assert report["n_unresolved"] == report["n_novel"]
    assert report["n_near_healthy"] >= int(0.8 * report["n_novel"])
    assert report["mean_peak_abs_log2fc"] < 0.6
    assert report["acetone_sep_t2d_minus_novel_ppb"] > 50
    assert report["passed"] is True
