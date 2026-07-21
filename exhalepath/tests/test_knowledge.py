from exhalepath.knowledge.loader import KnowledgeBase


def test_disease_alias_resolution():
    kb = KnowledgeBase()
    d = kb.resolve_disease("TCGA-PAAD")
    assert d["disease_id"] == "pancreatic_adenocarcinoma"
    assert "KRAS" in str(kb.pathways["Kras_mapk_proliferation"]["seed_genes"])


def test_voc_catalog_has_ppb_baselines():
    kb = KnowledgeBase()
    assert "acetone" in kb.vocs
    assert kb.vocs["acetone"]["healthy_ppb_median"] > 0
    assert len(kb.pathways) >= 8
