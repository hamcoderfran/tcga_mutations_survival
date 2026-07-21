from exhalepath.knowledge.loader import KnowledgeBase
from exhalepath.model.predict import ExhalePathPredictor
from exhalepath.schemas import DiseaseQuery, TumorContext


def test_dysbiosis_and_brain_diseases_resolve():
    kb = KnowledgeBase()
    assert kb.resolve_disease("dysbiosis")["disease_id"] == "gut_dysbiosis"
    assert kb.resolve_disease("SIBO")["disease_id"] == "sibo"
    assert kb.resolve_disease("Alzheimer's disease")["disease_id"] == "alzheimer_disease"
    assert kb.resolve_disease("Parkinson's")["disease_id"] == "parkinson_disease"
    assert kb.resolve_disease("depression")["disease_id"] == "major_depressive_disorder"
    assert kb.resolve_disease("glioblastoma")["disease_id"] == "glioblastoma"
    assert kb.resolve_disease("multiple sclerosis")["disease_id"] == "multiple_sclerosis"


def test_dysbiosis_elevates_microbial_vocs():
    pred = ExhalePathPredictor(use_opentargets=False)
    df = pred.predict(DiseaseQuery(disease="gut dysbiosis")).to_dataframe().set_index("voc_id")
    assert df.loc["indole", "fold_change"] >= 1.8
    assert df.loc["hydrogen_sulfide", "fold_change"] >= 1.6
    assert df.loc["phenol", "fold_change"] > 1.0


def test_alzheimer_and_parkinson_oxidative_signature():
    pred = ExhalePathPredictor(use_opentargets=False)
    ad = pred.predict(
        DiseaseQuery(
            disease="Alzheimer's",
            tumor=TumorContext(primary_site="brain"),
            mutated_genes=["APP", "PSEN1", "TREM2"],
        )
    ).to_dataframe().set_index("voc_id")
    pd_df = pred.predict(
        DiseaseQuery(
            disease="Parkinson's disease",
            tumor=TumorContext(primary_site="brain"),
            mutated_genes=["SNCA", "PRKN"],
        )
    ).to_dataframe().set_index("voc_id")
    assert ad.loc["hexanal", "fold_change"] >= 1.4
    assert pd_df.loc["pentane", "fold_change"] >= 1.4
    assert any(p.pathway_id == "neuroinflammation" and p.score > 0 for p in pred.predict(
        DiseaseQuery(disease="Alzheimer's", mutated_genes=["TREM2"])
    ).bundle.pathway_scores)


def test_microbiome_vocs_in_catalog():
    kb = KnowledgeBase()
    assert "indole" in kb.vocs
    assert "phenol" in kb.vocs
    assert "gut_microbiome_fermentation" in kb.pathways
    assert "neuroinflammation" in kb.pathways
