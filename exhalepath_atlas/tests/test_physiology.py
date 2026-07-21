from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache
from exhalepath.model.predict import ExhalePathPredictor
from exhalepath.physio.transport import farhi_alveolar_fraction
from exhalepath.schemas import DiseaseQuery, TumorContext


def setup_module():
    clear_knowledge_cache()


def test_farhi_high_lambda_less_alveolar_release():
    # Acetone (λ~340) retained more than pentane (λ~0.4)
    f_acetone = farhi_alveolar_fraction(340.0, va=5.0, q=5.5)
    f_pentane = farhi_alveolar_fraction(0.4, va=5.0, q=5.5)
    assert f_acetone < f_pentane


def test_diabetes_acetone_physiology_chain():
    pred = ExhalePathPredictor(use_opentargets=False)
    r = pred.predict(DiseaseQuery(disease="type 2 diabetes", mode="physiology"))
    df = r.to_dataframe().set_index("voc_id")
    assert df.loc["acetone", "fold_change"] >= 2.0
    trace = next(p.physiology for p in r.bundle.predictions if p.voc_id == "acetone")
    assert trace is not None
    assert trace.blood_delivery > 0
    assert "acetone_ketogenesis" in trace.contributing_chains
    assert any(s.state_id == "hepatocyte_ketogenic" for s in r.bundle.cell_states)


def test_cancer_acetaldehyde_warburg_physiology():
    pred = ExhalePathPredictor(use_opentargets=False)
    r = pred.predict(
        DiseaseQuery(
            disease="LUAD",
            mode="physiology",
            tumor=TumorContext(stage="III", primary_site="lung", histology="adenocarcinoma"),
            mutated_genes=["KRAS", "TP53", "HK2", "LDHA"],
        )
    )
    df = r.to_dataframe().set_index("voc_id")
    assert df.loc["acetaldehyde", "fold_change"] > 1.15
    trace = next(p.physiology for p in r.bundle.predictions if p.voc_id == "acetaldehyde")
    assert "acetaldehyde_warburg" in (trace.contributing_chains if trace else [])


def test_single_cell_fraction_override_increases_tumor_vocs():
    pred = ExhalePathPredictor(use_opentargets=False)
    base = pred.predict(
        DiseaseQuery(disease="PAAD", mode="physiology", mutated_genes=["KRAS"])
    )
    boosted = pred.predict(
        DiseaseQuery(
            disease="PAAD",
            mode="physiology",
            mutated_genes=["KRAS"],
            cell_state_fractions={"tumor_epithelial_warburg": 0.85},
            cell_state_activity={"tumor_epithelial_warburg": 1.5},
        )
    )
    b = base.to_dataframe().set_index("voc_id")
    u = boosted.to_dataframe().set_index("voc_id")
    assert u.loc["acetaldehyde", "predicted_ppb"] >= b.loc["acetaldehyde", "predicted_ppb"]


def test_pathway_chains_cover_major_vocs():
    kb = KnowledgeBase()
    vocs = {c["voc_id"] for c in kb.pathway_chains}
    for v in ["acetone", "acetaldehyde", "hexanal", "indole", "ammonia", "dms"]:
        assert v in vocs
