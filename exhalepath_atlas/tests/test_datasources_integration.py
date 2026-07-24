from pathlib import Path

from exhalepath.datasources import DATASOURCE_REGISTRY, integrate_all_datasources
from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache


def test_registry_has_priorities_1_to_15():
    root = Path(__file__).resolve().parents[1]
    reg = DATASOURCE_REGISTRY(root)
    pris = sorted(s.priority for s in reg)
    assert pris == list(range(1, 16))


def test_integrate_all_offline(tmp_path: Path):
    # Run against the atlas tree (writes into data/knowledge — expected for integration)
    root = Path(__file__).resolve().parents[1]
    man = integrate_all_datasources(root=root, offline=True)
    assert man["n_datasources"] == 15
    assert (root / "data" / "datasources" / "INTEGRATION_MANIFEST.json").exists()
    assert (root / "data" / "knowledge" / "atlas_capability.json").exists()
    clear_knowledge_cache()
    kb = KnowledgeBase(knowledge_dir=root / "data" / "knowledge")
    assert len(kb.diseases) >= 100
    assert len(kb.datasources) >= 10
    assert "hbdb" in kb.datasources or (root / "data" / "knowledge" / "datasource_hbdb.json").exists()
    assert "kegg" in kb.datasources or (root / "data" / "knowledge" / "datasource_kegg.json").exists()
    assert (root / "data" / "knowledge" / "datasource_alt_breath_sources.json").exists()
