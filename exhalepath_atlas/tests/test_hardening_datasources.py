"""Hardening-pass checks for MW quantified pulls, HBDB/VOLATILOME, and KEGG."""

from pathlib import Path

from exhalepath.datasources.ds01_metabolomics import MetabolomicsReposSource
from exhalepath.datasources.ds13_hbdb import HBDBSource
from exhalepath.datasources.ds14_kegg import KEGGSource, VOC_KEGG_SEED
from exhalepath.knowledge.loader import clear_knowledge_cache
from exhalepath.explain.mechanisms import MechanismExplainer


ROOT = Path(__file__).resolve().parents[1]


def test_metabolomics_offline_fuse_shape():
    src = MetabolomicsReposSource(ROOT)
    src.harvest(offline=True)
    info = src.fuse(ROOT / "data" / "knowledge")
    assert info["n_studies"] >= 5
    doc = (ROOT / "data" / "knowledge" / "datasource_metabolomics.json").read_text()
    assert "ST003200" in doc or "studies" in doc


def test_hbdb_offline_or_cached():
    src = HBDBSource(ROOT)
    paths = src.harvest(offline=True)
    assert paths["compounds"].exists()
    info = src.fuse(ROOT / "data" / "knowledge")
    assert info["n_compounds"] >= 40
    hbdb = (ROOT / "data" / "knowledge" / "datasource_hbdb.json")
    assert hbdb.exists()


def test_kegg_offline_seed_maps():
    src = KEGGSource(ROOT)
    # Prefer cached online harvest if present; else seed stub
    paths = src.harvest(offline=True)
    assert paths["maps"].exists()
    info = src.fuse(ROOT / "data" / "knowledge")
    assert info["n_mapped"] >= len(VOC_KEGG_SEED) // 2 or info["n_mapped"] >= 10
    assert (ROOT / "data" / "knowledge" / "datasource_kegg.json").exists()


def test_explain_includes_kegg_or_hbdb_when_fused():
    clear_knowledge_cache()
    # Ensure fused fragments exist
    HBDBSource(ROOT).harvest(offline=True)
    HBDBSource(ROOT).fuse(ROOT / "data" / "knowledge")
    KEGGSource(ROOT).harvest(offline=True)
    KEGGSource(ROOT).fuse(ROOT / "data" / "knowledge")
    clear_knowledge_cache()
    ex = MechanismExplainer(use_opentargets=False)
    pack = ex.build_disease_pack("type 2 diabetes", location="pancreas", top_n=8, mode="hybrid")
    assert pack["n_vocs_explained"] >= 1
    # At least one mechanism should carry kegg or hbdb evidence after fuse
    has_ev = any(
        (m.get("kegg") and m["kegg"].get("kegg_compound_id"))
        or (m.get("hbdb") and (m["hbdb"].get("volatilome") or m["hbdb"].get("in_atlas_priors")))
        for m in pack.get("voc_mechanisms") or []
    )
    assert has_ev
