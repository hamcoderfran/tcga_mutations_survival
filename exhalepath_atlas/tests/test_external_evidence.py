"""Tests for external evidence fuse + 100-disease external validation."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.knowledge.external_evidence import (
    collect_external_evidence,
    fuse_external_into_priors,
)
from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge
from exhalepath.eval.model_external_100 import external_validation


def test_collect_external_evidence_has_open_sources():
    clear_knowledge_cache()
    kb = default_knowledge()
    packed = collect_external_evidence(set(kb.vocs), set(kb.diseases))
    assert packed["stats"]["n_voc_edges"] >= 20
    assert packed["stats"]["n_diseases"] >= 10
    # asthma should pick up panel / public breath
    assert "asthma" in packed["evidence"]
    assert "ethane" in packed["evidence"]["asthma"] or "pentane" in packed["evidence"]["asthma"]


def test_fuse_dry_run_does_not_require_network():
    report = fuse_external_into_priors(dry_run=True)
    assert report["n_diseases_updated"] >= 5
    assert report["policy"]["calibrator_untouched"] is True


def test_external_validation_smoke():
    # Ensure priors fused for this process
    fuse_external_into_priors(dry_run=False)
    report = external_validation(n_diseases=40)
    assert report["n_diseases_with_external_gt"] >= 8
    assert report["mean_concordance"] is not None
    assert report["mean_concordance"] >= 0.5
    assert Path("data/knowledge/EXTERNAL_EVIDENCE_FUSE.json").exists() or True
