from exhalepath.model.predict import ExhalePathPredictor
from exhalepath.schemas import DiseaseQuery, TumorContext


def test_pancreatic_kras_tp53_predicts_elevated_ketones_aldehydes():
    predictor = ExhalePathPredictor(use_opentargets=False)
    result = predictor.predict(
        DiseaseQuery(
            disease="PAAD",
            tumor=TumorContext(
                stage="III",
                primary_site="pancreas",
                histology="adenocarcinoma",
                metastatic=False,
            ),
            mutated_genes=["KRAS", "TP53", "CDKN2A"],
        )
    )
    df = result.to_dataframe().set_index("voc_id")
    assert result.bundle.disease_id == "pancreatic_adenocarcinoma"
    assert df.loc["2_butanone", "fold_change"] > 1.2
    assert df.loc["hexanal", "predicted_ppb"] > df.loc["hexanal", "healthy_ppb"]
    assert any(p.score > 0 for p in result.bundle.pathway_scores)


def test_unknown_disease_still_returns_full_voc_panel():
    predictor = ExhalePathPredictor(use_opentargets=False)
    result = predictor.predict(
        DiseaseQuery(disease="obscure mitochondrial myopathy XYZ", mutated_genes=["MT-CO1"])
    )
    assert len(result.bundle.predictions) >= 50
    assert result.bundle.metadata["calibrator_loaded"] in {True, False}


def test_stage_iv_increases_burden_sensitive_vocs():
    predictor = ExhalePathPredictor(use_opentargets=False)
    early = predictor.predict(
        DiseaseQuery(
            disease="LUAD",
            tumor=TumorContext(stage="I", primary_site="lung"),
            mutated_genes=["TP53"],
        )
    )
    late = predictor.predict(
        DiseaseQuery(
            disease="LUAD",
            tumor=TumorContext(stage="IV", primary_site="lung", metastatic=True),
            mutated_genes=["TP53"],
        )
    )
    e = early.to_dataframe().set_index("voc_id")
    l = late.to_dataframe().set_index("voc_id")
    assert l.loc["hexanal", "predicted_ppb"] >= e.loc["hexanal", "predicted_ppb"]
