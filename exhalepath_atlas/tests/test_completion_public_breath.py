from exhalepath.eval.public_breath import evaluate_public_breath
from exhalepath.knowledge.loader import clear_knowledge_cache


def test_literature_public_breath_meets_completion_bar():
    clear_knowledge_cache()
    report = evaluate_public_breath(top_k=15, mode="hybrid")
    o = report["overall"]
    assert o["literature_elevated_recall_at_k"] is not None
    assert o["literature_elevated_recall_at_k"] >= 0.85
    assert o["literature_directional_accuracy"] >= 0.90
    assert (o["scientific_data_elevated_directional_accuracy"] or 0) >= 0.60
