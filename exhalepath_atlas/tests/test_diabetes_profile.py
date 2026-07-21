"""Diabetes clinical VOC profile + literature evaluation tests."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.diabetes_clinical_profile import (  # noqa: E402
    LIT_ACETONE_MIN_FOLD,
    run_profile,
)


def test_diabetes_profile_passes_literature_gate():
    profile = run_profile(out_dir=ROOT / "runs/diabetes_clinical_profile")
    assert profile["summary"]["overall_passed"] is True
    primary = profile["summary"]["primary_literature_eval"]
    assert primary["acetone_pass"] is True
    assert primary["directional_accuracy"] >= 0.67
    poorly = next(p for p in profile["phases"] if p["id"] == "poorly_controlled")
    assert poorly["ketone_panel"]["acetone_fold"] >= LIT_ACETONE_MIN_FOLD


def test_diabetes_tempo_increases_acetone():
    profile = run_profile(out_dir=ROOT / "runs/diabetes_clinical_profile")
    tempo = profile["summary"]["tempo"]
    pre = tempo["prediabetes"]["acetone_fold"]
    controlled = tempo["controlled_t2d"]["acetone_fold"]
    poorly = tempo["poorly_controlled"]["acetone_fold"]
    assert pre is not None and controlled is not None and poorly is not None
    assert controlled >= pre
    assert poorly >= controlled * 0.95  # allow plateau at physio cap


def test_diabetes_cli():
    from typer.testing import CliRunner

    from exhalepath.cli import app

    runner = CliRunner()
    result = runner.invoke(app, ["diabetes-profile", "--age", "55", "--sex", "male"])
    assert result.exit_code == 0, result.output
    assert "PASSED" in result.output
    assert "acetone" in result.output.lower()
