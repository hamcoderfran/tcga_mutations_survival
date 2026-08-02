"""Smoke tests for stack holdout + adversarial parser suite."""

from __future__ import annotations

from exhalepath.eval.stack_holdout import (
    eval_naturalistic_to_stack_smoke,
    eval_patient_template_break,
)
from exhalepath.nl import parse_patient_template


def test_multi_disease_prefers_depression():
    tpl = parse_patient_template(
        "depression with obesity and type 2 diabetes", llm="rules"
    )
    assert tpl.primary_disease
    assert "depress" in tpl.primary_disease.lower()


def test_gene_only_does_not_invent_disease():
    tpl = parse_patient_template("Patient has KRAS and TP53 mutations", llm="rules")
    assert "KRAS" in tpl.genes and "TP53" in tpl.genes
    assert not tpl.primary_disease


def test_adversarial_hard_pass_rate():
    report = eval_patient_template_break()
    assert report["hard_total"] >= 10
    assert report["hard_pass_rate"] >= 0.9


def test_naturalistic_e2e_smoke():
    report = eval_naturalistic_to_stack_smoke()
    assert report["n_cases"] == 3
    assert report["mean_models_ok"] >= 14
