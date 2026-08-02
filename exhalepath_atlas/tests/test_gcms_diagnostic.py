"""Tests for GC-MS patient-level diagnostic research enablement."""

from __future__ import annotations

from pathlib import Path

from exhalepath.gcms import (
    disease_signature,
    load_mw_patient_matrix,
    lock_split,
    score_observed_vector,
    score_patients,
    verify_split_manifest,
)
from exhalepath.eval.patient_diagnostic import evaluate_patient_diagnostic


def test_load_st000883_patient_matrix():
    m = load_mw_patient_matrix("ST000883")
    assert m.disease_id == "malaria"
    assert m.matrix.shape[0] >= 30
    assert m.matrix.shape[1] >= 5
    assert set(m.labels.unique()) == {0, 1}
    assert int(m.labels.sum()) >= 10


def test_lock_split_integrity(tmp_path: Path):
    m = load_mw_patient_matrix("ST000883")
    path = tmp_path / "split.json"
    manifest = lock_split(m, strategy="stratified_kfold", n_splits=5, seed=42, out_path=path)
    assert path.exists()
    assert len(manifest.content_sha256) == 64
    assert verify_split_manifest(manifest, m) == []
    # no train/test leakage
    for fold in manifest.folds:
        assert not (set(fold.train_subject_ids) & set(fold.test_subject_ids))


def test_signature_scoring_runs():
    m = load_mw_patient_matrix("ST000883").log1p()
    sig = disease_signature("malaria", source="hybrid", voc_ids=list(m.matrix.columns))
    assert len(sig) >= 2
    scores = score_patients(m, sig, method="cosine")
    assert len(scores) == len(m.subject_ids)
    one = score_observed_vector(
        {"acetone": 0.3, "pentane": 0.4, "isoprene": -0.2},
        sig,
        method="cosine",
    )
    assert one["n_overlap"] >= 2
    assert one["score"] is not None


def test_evaluate_patient_diagnostic_malaria(tmp_path: Path):
    report = evaluate_patient_diagnostic(
        study_id="ST000883",
        signature_source="hybrid",
        out_dir=tmp_path / "diag",
        n_splits=5,
        seed=42,
    )
    m = report["metrics"]
    assert m["n_subjects"] >= 30
    assert m["auroc"] is not None
    assert 0.0 <= m["auroc"] <= 1.0
    assert (tmp_path / "diag" / "PATIENT_DIAGNOSTIC_REPORT.md").exists()
    assert (tmp_path / "diag" / "roc_curve.png").exists()
    assert report["nested"]["content_sha256"]
