"""Real-breath corpus: no synthetic labels; priority diseases covered."""

from __future__ import annotations

import json
from pathlib import Path

from exhalepath.config import PROCESSED_DIR, DATA_DIR
from exhalepath.ingest.real_breath_corpus import (
    PRIORITY_DISEASES,
    build_real_training_corpus,
    collect_real_label_rows,
)
from exhalepath.knowledge.loader import KnowledgeBase


def test_no_synthetic_in_manifest():
    man = json.loads((PROCESSED_DIR / "corpus_manifest.json").read_text())
    assert man.get("corpus_type") == "real_breath_only"
    assert man.get("synthetic_labels") is False
    assert man.get("n_voc_targets", 0) > 50
    assert man.get("n_voc_targets", 0) < 50_000  # not the old 600k synthetic dump


def test_priority_diseases_have_labels_or_atlas():
    labels = collect_real_label_rows()
    labeled = set(labels["disease_id"].astype(str))
    kb = KnowledgeBase()
    for did in PRIORITY_DISEASES:
        d = kb.resolve_disease(did)
        assert d.get("disease_id") not in {None, "unknown"}, did
        assert did in labeled or len(d.get("voc_log2fc_prior") or {}) > 0


def test_offline_demo_blocked():
    import os
    from exhalepath.ingest.build_corpus import build_training_corpus

    os.environ.pop("VOC_ALLOW_SYNTHETIC", None)
    try:
        build_training_corpus(offline_demo=True)
        assert False, "should have raised"
    except RuntimeError as e:
        assert "Synthetic" in str(e) or "synthetic" in str(e).lower()


def test_malaria_and_ards_resolve():
    kb = KnowledgeBase()
    assert kb.resolve_disease("malaria")["disease_id"] == "malaria"
    assert kb.resolve_disease("ARDS")["disease_id"] == "ards"
    assert kb.resolve_disease("VAP")["disease_id"] == "pneumonia_bacterial"
