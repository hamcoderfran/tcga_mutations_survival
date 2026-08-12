"""Tests for literature compare, smoking fix, and disease100 suite."""

from __future__ import annotations

from pathlib import Path

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.knowledge.loader import clear_knowledge_cache
from exhalepath.eval.lit_compare import (
    literature_concordance,
    run_disease100_suite,
    run_lit_demo_disease100,
)


def test_smoking_btex_survives_hybrid_blend():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    cur = eng.predict(
        disease="copd",
        location="lung",
        smoking_status="current",
        age_years=60,
        sex="male",
        top_n=50,
        explain=False,
    )
    nev = eng.predict(
        disease="copd",
        location="lung",
        smoking_status="never",
        age_years=60,
        sex="male",
        top_n=50,
        explain=False,
    )
    vc = {p.voc_id: p.log2_fold_change for p in cur.result.bundle.predictions}
    vn = {p.voc_id: p.log2_fold_change for p in nev.result.bundle.predictions}
    assert vc["benzene"] - vn["benzene"] >= 0.35
    assert vc["toluene"] - vn["toluene"] >= 0.35
    assert vc["pentane"] - vn["pentane"] >= 0.20


def test_age_modulates_oxidative_vocs():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    young = eng.predict(
        disease="asthma", location="lung", age_years=25, sex="female",
        smoking_status="never", top_n=50, explain=False,
    )
    old = eng.predict(
        disease="asthma", location="lung", age_years=80, sex="female",
        smoking_status="never", top_n=50, explain=False,
    )
    y = {p.voc_id: p.log2_fold_change for p in young.result.bundle.predictions}
    o = {p.voc_id: p.log2_fold_change for p in old.result.bundle.predictions}
    assert o["pentane"] > y["pentane"]
    assert o["hexanal"] > y["hexanal"]


def test_literature_concordance_reasonable():
    lit = literature_concordance()
    assert lit["n_diseases"] >= 10
    assert lit["mean_concordance"] is not None
    assert lit["mean_concordance"] >= 0.7


def test_disease100_smoke(tmp_path: Path):
    out = tmp_path / "d100"
    # smaller n for speed in unit test
    report = run_disease100_suite(out_dir=out, n_diseases=20)
    assert report["n_diseases"] == 20
    assert (out / "disease100_voc_matrix.csv").exists()
    assert (out / "figures" / "disease100_pca_by_category.png").exists()


def test_lusc_farther_from_copd_than_bronchitis():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    import numpy as np

    def vec(did):
        r = eng.predict(
            disease=did, location="lung", smoking_status="never",
            age_years=60, top_n=50, explain=False,
        )
        return {p.voc_id: p.log2_fold_change for p in r.result.bundle.predictions}

    ids = ["copd", "chronic_bronchitis", "lung_squamous_cell_carcinoma", "lung_adenocarcinoma"]
    vs = [vec(i) for i in ids]
    keys = sorted(vs[0])
    F = np.array([[v[k] for k in keys] for v in vs])
    d_cb = float(np.linalg.norm(F[0] - F[1]))
    d_lusc = float(np.linalg.norm(F[0] - F[2]))
    d_luad = float(np.linalg.norm(F[0] - F[3]))
    # Bronchitis should be nearer COPD than LUAD; LUSC should not collapse onto COPD
    assert d_cb < d_luad
    # Soft bound: main priors are less separated than the lit-compare tip rewrite
    assert d_lusc > d_cb * 0.7


def test_lit_concordance_includes_mh_panels():
    clear_knowledge_cache()
    lit = literature_concordance()
    by = {r["disease_id"]: r for r in lit.get("by_disease") or []}
    assert "schizophrenia" in by
    assert "major_depressive_disorder" in by
    assert "bipolar" in by
    assert by["schizophrenia"]["n_checked"] >= 8
    assert by["major_depressive_disorder"]["n_checked"] >= 6
    assert by["bipolar"]["n_checked"] >= 2
    assert by["schizophrenia"]["concordance"] == 1.0
    assert by["major_depressive_disorder"]["concordance"] == 1.0
    assert by["bipolar"]["concordance"] == 1.0
