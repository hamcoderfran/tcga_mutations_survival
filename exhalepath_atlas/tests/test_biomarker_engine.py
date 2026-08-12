from exhalepath import ExhaleBiomarkerEngine
from exhalepath.body.tissues import WholeBodyMap, resolve_location
from exhalepath.knowledge.loader import clear_knowledge_cache, default_knowledge


def test_atlas_has_core_vocs_and_100_diseases():
    clear_knowledge_cache()
    kb = default_knowledge()
    assert len(kb.vocs) >= 50
    assert "butyric_acid" in kb.vocs
    assert len(kb.diseases) >= 100
    assert len(kb.tissues) >= 50


def test_resolve_any_body_location():
    lung = resolve_location("lung")
    assert lung["matched"] is True
    assert lung["tissue_id"] == "lung"
    breast = resolve_location("left breast upper outer")
    assert breast["matched"] is True
    assert "breast" in breast["tissue_id"]
    brain = resolve_location("brain")
    assert brain["tissue_id"] == "brain"
    custom = resolve_location("xyzzy_unmatched_organ_999")
    assert custom["matched"] is False
    assert custom["name"]
    # Free-text containing a known alias still resolves (e.g. retroperitoneal → peritoneum)
    peri = resolve_location("retroperitoneal soft tissue mass")
    assert peri["matched"] is True
    assert peri["tissue_id"] == "peritoneum"


def test_biomarker_top50_any_location():
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    report = engine.predict(
        "lung adenocarcinoma",
        location="lung",
        genes=["KRAS", "TP53"],
        top_n=50,
    )
    assert report.n_vocs_modeled >= 50
    assert len(report.top_vocs) == 50
    assert report.tissue == "lung"
    # Ranked by |Δppb|
    deltas = [abs(p.delta_ppb) for p in report.ranked]
    assert deltas == sorted(deltas, reverse=True)
    # Cancer lipid / glycolysis panel should move
    by_id = {p.voc_id: p for p in report.top_vocs}
    assert abs(by_id["hexanal"].delta_ppb) > 0.5
    assert abs(by_id["acetaldehyde"].delta_ppb) > 0.5
    df = report.to_dataframe()
    assert list(df["rank"]) == list(range(1, 51))


def test_biomarker_brain_and_gut_locations():
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    ad = engine.predict("Alzheimer's disease", location="brain", top_n=50)
    assert ad.location["tissue_id"] == "brain"
    assert ad.n_vocs_modeled >= 50
    assert abs(ad.top_vocs[0].delta_ppb) > 0

    gut = engine.predict("gut dysbiosis", location="small intestine", top_n=20)
    assert gut.location["matched"] is True
    by_id = {p.voc_id: p for p in gut.top_vocs}
    # Microbiome VOCs should appear among changers
    assert any(
        abs(by_id[v].delta_ppb) > 0.05
        for v in ("indole", "phenol", "hydrogen_sulfide", "trimethylamine")
        if v in by_id
    )


def test_list_diseases_and_locations():
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    diseases = engine.list_diseases()
    locs = engine.list_locations()
    assert len(diseases) >= 100
    assert len(locs) == len(WholeBodyMap().tissues)
    assert any(d["disease_id"] == "lung_adenocarcinoma" for d in diseases)
