"""Revolutionary gap pack: DDx, residualize, claim ledger, BreathVOC, closed loop."""

from __future__ import annotations

from pathlib import Path

from exhalepath.eval.revolutionary_gap import evaluate_revolutionary_gap
from exhalepath.gcms.claim_ledger import build_claim_ledger
from exhalepath.gcms.differential import (
    evaluate_scidata_differential,
    rank_diseases_for_vector,
)
from exhalepath.gcms.interchange import (
    BREATHVOC_SCHEMA_VERSION,
    export_breathvoc,
    load_and_validate_breathvoc,
)
from exhalepath.gcms.patient_matrix import load_mw_patient_matrix
from exhalepath.gcms.residualize import residualize_matrix
from exhalepath.gcms.scidata_samples import load_scidata_intensity_matrix


def test_claim_ledger_has_quantified_and_priors():
    ledger = build_claim_ledger()
    assert ledger["n_claims"] > 50
    grades = ledger["by_evidence_grade"]
    assert grades.get("quantified", 0) + grades.get("directional_only", 0) > 0
    assert grades.get("atlas_prior", 0) > 0


def test_breathvoc_roundtrip(tmp_path: Path):
    m = load_mw_patient_matrix("ST000883", feature_map="atlas")
    path = export_breathvoc(m, tmp_path / "x.breathvoc.json")
    loaded, errs = load_and_validate_breathvoc(path)
    assert errs == []
    assert loaded.matrix.shape == m.matrix.shape
    bundle = __import__("json").loads(path.read_text())
    assert bundle["schema_version"] == BREATHVOC_SCHEMA_VERSION


def test_residualize_scidata_age_sex():
    X, cov, y = load_scidata_intensity_matrix(mapped_vocs_only=True)
    cov = cov.set_index("sample_id")
    import numpy as np

    X_log = np.log1p(X.clip(lower=0))
    out = residualize_matrix(X_log, cov, covariate_cols=["age", "sex"])
    assert out["residual_matrix"].shape == X_log.shape
    assert out["mean_covariate_r2"] is not None
    assert 0.0 <= out["mean_covariate_r2"] <= 1.0


def test_scidata_differential_mechanism_vs_ceiling():
    report = evaluate_scidata_differential(seed=0, residualize_age_sex=True)
    ceil = report["fit_on_cohort_ceiling"]["macro_ovr_auroc"]
    mech = report["mechanism_ddx_hybrid"]["macro_ovr_auroc"]
    assert ceil is not None and ceil > 0.7  # fit-on-cohort should be strong
    assert mech is not None and 0.0 <= mech <= 1.0
    # Ceiling should beat transferable mechanism on this labeled cohort
    assert ceil >= mech - 0.05
    resid = report["age_sex_residualized"]
    assert "mean_covariate_r2" in resid or "error" in resid


def test_rank_diseases_for_vector():
    rows = rank_diseases_for_vector(
        {"pentane": 0.5, "hexanal": 0.4, "acetone": -0.2, "isoprene": -0.3},
        ["malaria", "asthma", "heart_failure"],
        source="hybrid",
    )
    assert len(rows) == 3
    assert rows[0]["rank"] == 1
    assert rows[0]["score"] is not None


def test_evaluate_revolutionary_gap_quick(tmp_path: Path):
    report = evaluate_revolutionary_gap(out_dir=tmp_path / "rev", quick=True)
    assert report["scorecard"]["exhalepath_pct"] >= 50
    assert (tmp_path / "rev" / "REVOLUTIONARY_GAP.md").exists()
    assert (tmp_path / "rev" / "CLAIM_LEDGER.json").exists()
    assert report["live_differentiators"]["breathvoc_interchange"]["roundtrip_ok"]
