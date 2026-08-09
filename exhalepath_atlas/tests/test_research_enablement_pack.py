"""Sci Data per-sample adapters, filters, MetaboLights export, paper pack."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from exhalepath.eval.scidata_samples import evaluate_scidata_samples
from exhalepath.gcms import (
    blank_ratio_filter,
    export_metabolights_bundle,
    export_paper_pack,
    list_scidata_cohorts,
    load_mw_patient_matrix,
    load_scidata_ovr_matrix,
)
from exhalepath.gcms.metrics import stratified_auroc
from exhalepath.gcms.preprocess import apply_feature_filter
from exhalepath.gcms.scidata_samples import load_scidata_intensity_matrix


def test_scidata_cohorts_available():
    rows = list_scidata_cohorts()
    cohorts = {r["cohort"] for r in rows if "cohort" in r}
    assert {"asthma", "copd", "bronchiectasis"} <= cohorts
    assert all(r.get("available") for r in rows if "cohort" in r)


def test_scidata_per_sample_matrix_shape():
    m = load_scidata_ovr_matrix("asthma", mapped_vocs_only=True)
    assert m.matrix.shape[0] >= 100  # asthma+copd+bronchi samples
    assert m.matrix.shape[1] >= 3
    assert set(m.labels.unique()) == {0, 1}
    assert int(m.labels.sum()) >= 40  # asthma positives
    # age/sex joined for most subjects
    ages = [s.age for s in m.samples if s.age is not None]
    sexes = [s.sex for s in m.samples if s.sex is not None]
    assert len(ages) >= 80
    assert len(sexes) >= 80
    assert m.metadata.get("smoking_available") is False


def test_blank_ratio_and_detection_filter():
    intensity, _, _ = load_scidata_intensity_matrix(mapped_vocs_only=False)
    # no blanks in Sci Data → detection filter only
    filt = blank_ratio_filter(intensity, min_detect_frac=0.5, min_ratio=2.0)
    assert filt["n_blanks"] == 0
    assert len(filt["kept_features"]) >= 1
    kept = apply_feature_filter(intensity, filt["kept_features"])
    assert kept.shape[1] == len(filt["kept_features"])

    # synthetic blank row should drop low-ratio features
    fake = intensity.copy()
    blank = fake.median(axis=0) * 10.0
    fake.loc["BLANK_01"] = blank
    filt2 = blank_ratio_filter(fake, blank_ids=["BLANK_01"], min_ratio=2.0, min_detect_frac=0.1)
    assert filt2["n_blanks"] == 1
    assert len(filt2["dropped_features"]) >= 1


def test_stratified_auroc_age_sex():
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1])
    s = np.array([0.1, 0.2, 0.15, 0.3, 0.7, 0.8, 0.6, 0.9, 0.25, 0.55, 0.35, 0.65])
    ages = [20, 22, 40, 45, 60, 62, 25, 70, 30, 55, 28, 50]
    sex = ["male", "female"] * 6
    out = stratified_auroc(y, s, {"age": ages, "sex": sex}, min_n=4)
    assert out["overall"] is not None
    assert "age" in out["strata"]
    assert "sex" in out["strata"]


def test_export_metabolights_and_paper_pack(tmp_path: Path):
    m = load_mw_patient_matrix("ST000883")
    paths = export_metabolights_bundle(m, tmp_path / "mtbls")
    assert paths["mztab"].exists()
    assert paths["investigation"].exists()
    assert "mzTab-version" in paths["mztab"].read_text()

    # minimal run dir for paper pack
    run = tmp_path / "run"
    run.mkdir()
    (run / "PATIENT_DIAGNOSTIC_REPORT.md").write_text("# report\n")
    (run / "roc_curve.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"\0" * 32)
    (run / "ST000883_split_manifest.json").write_text(
        '{"dataset_id":"ST000883","content_sha256":"abc123"}'
    )
    info = export_paper_pack(
        run,
        study_id="ST000883",
        disease_id="malaria",
        n_subjects=35,
        split_sha256="abc123",
    )
    assert Path(info["zip_path"]).exists()
    assert info["n_files"] >= 2
    assert info["split_content_sha256"] == "abc123"
    assert (run / "METHODS.md").exists()


def test_evaluate_scidata_samples_asthma(tmp_path: Path):
    report = evaluate_scidata_samples(
        positive_cohort="asthma",
        out_dir=tmp_path / "asthma",
        n_splits=3,
        seed=0,
        min_detect_frac=0.2,
        make_paper_pack=True,
    )
    m = report["metrics"]
    assert m["n_subjects"] >= 100
    assert m["auroc"] is not None
    assert 0.0 <= m["auroc"] <= 1.0
    assert "age" in (report.get("stratified_auroc") or {}).get("strata", {})
    assert (tmp_path / "asthma" / "SCIDATA_SAMPLE_EVAL.md").exists()
    assert (tmp_path / "asthma" / "stratified_auroc.json").exists()
    assert report.get("paper_pack", {}).get("zip_path")
    assert Path(report["paper_pack"]["zip_path"]).exists()
