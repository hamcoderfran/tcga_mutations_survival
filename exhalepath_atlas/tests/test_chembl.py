from pathlib import Path

import pandas as pd

from exhalepath.ingest.chembl_harvest import (
    distill_chembl_priors,
    harvest_chembl_for_pathways,
)
from exhalepath.knowledge.loader import KnowledgeBase, clear_knowledge_cache
from exhalepath.model.features import build_feature_vector
from exhalepath.model.train_chembl import train_chembl_aux_model
from exhalepath.pathways.score import score_pathways
def test_offline_chembl_harvest_and_train(tmp_path: Path):
    paths = harvest_chembl_for_pathways(
        out_dir=tmp_path / "chembl",
        offline_demo=True,
        demo_rows=5_000,
        include_voc_properties=True,
        install_knowledge=False,
    )
    assert paths["activities"].exists()
    acts = pd.read_csv(paths["activities"])
    assert len(acts) == 5_000
    assert "pchembl_value" in acts.columns
    priors = distill_chembl_priors(acts, kb=KnowledgeBase())
    assert priors["n_source_activities"] == 5_000
    assert len(priors["pathways"]) >= 1

    result = train_chembl_aux_model(
        activities_path=paths["activities"],
        out_dir=tmp_path / "models",
        max_rows=5_000,
        update_knowledge_priors=False,
    )
    assert result["metrics"]["n_rows"] >= 100
    assert Path(result["model_path"]).exists()


def test_chembl_features_in_predict_vector(tmp_path: Path):
    paths = harvest_chembl_for_pathways(
        out_dir=tmp_path / "chembl",
        offline_demo=True,
        demo_rows=2_000,
        install_knowledge=True,
    )
    clear_knowledge_cache()
    kb = KnowledgeBase()
    assert kb.chembl_priors.get("pathways")
    disease = kb.resolve_disease("LUAD")
    scores = score_pathways(
        kb=kb, disease=disease, mutated_genes=["KRAS", "TP53"], tumor=None
    )
    feats = build_feature_vector(
        kb=kb,
        disease=disease,
        pathway_scores=scores,
        tumor=None,
        voc_id="hexanal",
    )
    assert "chembl_voc_ligandability" in feats
    assert any(k.startswith("chembl_lig_") for k in feats)
    assert paths["priors"].exists()
