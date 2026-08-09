"""CSIRO labels, LOSO, confounder PTR, stop-chasing artifacts."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from exhalepath.eval.confounder_ptr import evaluate_confounder_ptr, load_st003200_ptr_matrix
from exhalepath.eval.malaria_diagnostic import evaluate_malaria_diagnostic
from exhalepath.gcms.csiro_chmi import csiro_readiness, load_csiro_labels
from exhalepath.gcms.external_malaria import leave_one_study_out, load_external_patient_csv
from exhalepath.gcms.patient_matrix import load_mw_patient_matrix


def test_csiro_labels_bundled():
    ready = csiro_readiness()
    assert ready["labels_bundled"] is True
    assert ready["n_subjects"] == 7
    assert ready["n_labeled_baseline_peak"] == 14
    lab = load_csiro_labels()
    assert set(lab["label"].unique()) == {0, 1}
    assert (lab["label"] == 0).sum() == 7
    assert (lab["label"] == 1).sum() == 7


def test_leave_one_study_out_two_synthetic(tmp_path: Path):
    m = load_mw_patient_matrix("ST000883", feature_map="malaria_lit").log1p()
    # Stratified half-split so each pseudo-study has both classes
    pos = list(m.labels[m.labels == 1].index.astype(str))
    neg = list(m.labels[m.labels == 0].index.astype(str))
    a_ids = pos[: len(pos) // 2] + neg[: len(neg) // 2]
    b_ids = pos[len(pos) // 2 :] + neg[len(neg) // 2 :]
    assert len(a_ids) >= 8 and len(b_ids) >= 8

    def _slice(study_id, keep):
        mat = m.matrix.loc[keep]
        lab = m.labels.loc[keep]
        path_m = tmp_path / f"{study_id}_X.csv"
        path_y = tmp_path / f"{study_id}_y.csv"
        mat.to_csv(path_m)
        lab.to_frame("label").to_csv(path_y)
        return load_external_patient_csv(path_m, path_y, study_id=study_id)

    a = _slice("STUDY_A", a_ids)
    b = _slice("STUDY_B", b_ids)
    sig = {c: 1.0 for c in list(m.matrix.columns)[:6]}
    out = leave_one_study_out([a, b], signature=sig, max_features=4, seed=0)
    assert out["status"] == "ok"
    assert len(out["folds"]) == 2
    assert out["mean_auroc_sparse"] is not None or any(
        f.get("auroc_sparse") is not None for f in out["folds"]
    )


def test_evaluate_malaria_writes_stop_chasing_and_csiro(tmp_path: Path):
    report = evaluate_malaria_diagnostic(
        out_dir=tmp_path / "mal",
        feature_map="malaria_lit",
        n_splits=5,
        seed=2,
        max_features=6,
        loso=True,  # no external → skipped with CSIRO hint
    )
    assert report["csiro_chmi"]["labels_bundled"] is True
    assert report["leave_one_study_out"]["status"] == "skipped"
    assert report["story"]["transferable_nested_auroc"] is not None
    assert (tmp_path / "mal" / "STOP_CHASING_ST000883.md").exists()
    assert (tmp_path / "mal" / "csiro_chmi" / "csiro_chmi_labels.csv").exists()


def test_confounder_ptr_st003200(tmp_path: Path):
    matrix, demo = load_st003200_ptr_matrix()
    assert matrix.matrix.shape[0] >= 400
    assert "sex" in demo.columns and "smoking_status" in demo.columns
    report = evaluate_confounder_ptr(out_dir=tmp_path / "cf", include_scidata=False)
    c = report["st003200"]["confounder_proxy_auroc"]
    # Healthy PTR: smoking should leave a detectable VOC footprint
    assert c["smoking_current_vs_never"] is not None
    assert c["smoking_current_vs_never"] > 0.55
    assert (tmp_path / "cf" / "CONFOUNDER_PTR.md").exists()
