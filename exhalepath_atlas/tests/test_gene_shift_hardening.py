"""Gene drivers must shift pathway scores / VOC profiles."""

from __future__ import annotations

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.knowledge.loader import clear_knowledge_cache


def _l2(a, b) -> float:
    m1 = {p.voc_id: p.log2_fold_change for p in a.result.bundle.predictions}
    m0 = {p.voc_id: p.log2_fold_change for p in b.result.bundle.predictions}
    return sum((m1[k] - m0[k]) ** 2 for k in m1) ** 0.5


def test_copd_gwas_genes_shift_vocs():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    with_g = eng.predict(
        disease="copd",
        location="lung",
        genes=["SERPINA1", "HHIP", "CHRNA3", "FAM13A"],
        smoking_status="current",
        top_n=50,
        explain=False,
    )
    without = eng.predict(
        disease="copd",
        location="lung",
        genes=None,
        smoking_status="current",
        top_n=50,
        explain=False,
    )
    assert _l2(with_g, without) > 0.05


def test_bronchitis_resolves_and_predicts():
    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    r = eng.predict(
        disease="bronchitis",
        location="lung",
        genes=["TNF", "IL6"],
        smoking_status="former",
        top_n=20,
        explain=False,
    )
    assert r.disease_id == "chronic_bronchitis"
    by = {p.voc_id: p for p in r.result.bundle.predictions}
    assert by["pentane"].fold_change > 1.0
