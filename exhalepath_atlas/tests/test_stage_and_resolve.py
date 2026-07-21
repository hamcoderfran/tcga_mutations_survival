from exhalepath.knowledge.loader import KnowledgeBase
from exhalepath.utils import stage_multiplier, stage_ordinal


def test_stage_multiplier_longest_match():
    assert stage_multiplier("I") == 1.0
    assert stage_multiplier("II") == 1.15
    assert stage_multiplier("III") == 1.35
    assert stage_multiplier("IV") == 1.55
    assert stage_multiplier("Stage IV") == 1.55
    assert stage_multiplier("IIIA") == 1.35
    assert stage_ordinal("IV") > stage_ordinal("III") > stage_ordinal("II") > stage_ordinal("I")


def test_disease_resolve_does_not_overfit_cancer_substring():
    kb = KnowledgeBase()
    assert kb.resolve_disease("cancer").get("_unresolved") is True
    assert kb.resolve_disease("lung cancer")["disease_id"] == "lung_adenocarcinoma"
    assert kb.resolve_disease("PAAD")["disease_id"] == "pancreatic_adenocarcinoma"
    assert kb.resolve_disease("cirrhosis")["disease_id"] == "chronic_liver_disease"


def test_acsl4_not_acysl4():
    kb = KnowledgeBase()
    genes = kb.pathways["lipid_peroxidation"]["seed_genes"]
    assert "ACSL4" in genes
    assert "ACYSL4" not in genes
