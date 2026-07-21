from exhalepath.explain.mechanisms import MechanismExplainer
from exhalepath.physio.census_fractions import resolve_us_disease_id


def test_resolve_us_disease_for_top_diseases():
    assert resolve_us_disease_id("alzheimer_disease") == "alzheimer"
    assert resolve_us_disease_id("lung_adenocarcinoma") == "cancer_lung"
    assert resolve_us_disease_id("type2_diabetes") == "type2_diabetes"


def test_mechanism_pack_has_why_cells_genes():
    explainer = MechanismExplainer(use_opentargets=False)
    pack = explainer.build_disease_pack(
        "lung adenocarcinoma",
        location="lung",
        genes=["KRAS", "TP53"],
        top_n=8,
    )
    assert pack["disease_id"] == "lung_adenocarcinoma"
    assert pack["voc_mechanisms"]
    m0 = pack["voc_mechanisms"][0]
    assert m0["why"]
    assert "pathways" in m0
    assert "cell_states" in m0
    assert pack["top_pathways"]
    # Census should resolve for lung cancer harvest
    assert pack["census"].get("us_disease_id") == "cancer_lung"
