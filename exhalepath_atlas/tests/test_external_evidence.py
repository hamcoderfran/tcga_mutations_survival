"""Tests for external evidence fuse + 100-disease external validation."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.knowledge.external_evidence import (
    collect_external_evidence,
    enrich_voc_catalog_identity,
    fuse_external_into_priors,
    load_pristine_doc,
)
from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge
from exhalepath.eval.model_external_100 import external_validation


def test_collect_external_evidence_has_open_sources():
    clear_knowledge_cache()
    kb = default_knowledge()
    packed = collect_external_evidence(set(kb.vocs), set(kb.diseases))
    assert packed["stats"]["n_voc_edges"] >= 20
    assert packed["stats"]["n_diseases"] >= 10
    assert "asthma" in packed["evidence"]
    assert "ethane" in packed["evidence"]["asthma"] or "pentane" in packed["evidence"]["asthma"]


def test_fuse_dry_run_does_not_write_catalogs(tmp_path: Path):
    """dry_run must not mutate voc_catalog on disk."""
    clear_knowledge_cache()
    kb_paths = []
    from exhalepath.knowledge import external_evidence as ee

    for know in ee._know_dirs():
        p = know / "voc_catalog.json"
        if p.exists():
            kb_paths.append((p, p.read_text()))
    report = fuse_external_into_priors(dry_run=True)
    assert report["n_diseases_updated"] >= 5
    assert report["policy"]["calibrator_untouched"] is True
    assert report["policy"]["blend_from_pristine"] is True
    for p, before in kb_paths:
        assert p.read_text() == before


def test_fuse_idempotent_from_pristine():
    clear_knowledge_cache()
    a = fuse_external_into_priors(dry_run=False, include_lit_bench=True, backup=False)
    clear_knowledge_cache()
    b = fuse_external_into_priors(dry_run=False, include_lit_bench=True, backup=False)
    assert a["n_voc_prior_updates"] == b["n_voc_prior_updates"]
    # Sample a fused disease prior value is stable
    doc1 = json.loads((Path("data/knowledge/disease_voc_priors.json")).read_text())
    clear_knowledge_cache()
    fuse_external_into_priors(dry_run=False, include_lit_bench=True, backup=False)
    doc2 = json.loads((Path("data/knowledge/disease_voc_priors.json")).read_text())
    d1 = next(d for d in doc1["diseases"] if d["disease_id"] == "asthma")
    d2 = next(d for d in doc2["diseases"] if d["disease_id"] == "asthma")
    assert d1["voc_log2fc_prior"] == d2["voc_log2fc_prior"]


def test_pristine_snapshot_exists():
    doc = load_pristine_doc()
    assert len(doc.get("diseases") or []) >= 100
    assert not any(d.get("external_evidence_fused") for d in doc["diseases"])


def test_external_validation_includes_expect_diseases():
    fuse_external_into_priors(dry_run=False, backup=False)
    report = external_validation(n_diseases=100)
    assert report["n_diseases_with_external_gt"] >= 8
    assert report["mean_concordance"] is not None
    assert report["mean_concordance"] >= 0.5
    ids = {r["disease_id"] for r in report["by_disease"]}
    # Prefer-expect ordering must keep key diseases with GT
    assert "type_2_diabetes" in ids or "asthma" in ids
    assert report.get("circular_with_fuse") is True
    assert Path("data/knowledge/EXTERNAL_EVIDENCE_FUSE.json").exists()


def test_enrich_identity_dry_run_no_write():
    from exhalepath.knowledge import external_evidence as ee

    before = {}
    for know in ee._know_dirs():
        p = know / "voc_catalog.json"
        if p.exists():
            before[str(p)] = p.read_text()
    enrich_voc_catalog_identity(dry_run=True)
    for path, text in before.items():
        assert Path(path).read_text() == text
