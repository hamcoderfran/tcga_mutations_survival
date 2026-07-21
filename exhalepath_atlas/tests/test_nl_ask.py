"""Tests for natural-language parse + QuerySlots questionnaire path."""

from __future__ import annotations

from exhalepath.nl import parse_rules, parse_with_optional_llm, QuerySlots
from exhalepath.nl.ask import fill_slots_interactively


def test_rules_depression_obese_male():
    slots = parse_rules("24yo obese male with depression")
    assert slots.disease and "depress" in slots.disease.lower()
    assert "obesity" in [c.lower() for c in slots.comorbidities]
    assert slots.age_years == 24
    assert slots.sex == "male"
    assert slots.location == "brain"
    assert slots.parse_method == "rules"


def test_rules_luad_lll_genes():
    slots = parse_rules(
        "Stage II lung adenocarcinoma in the left lower lobe with KRAS and TP53"
    )
    assert slots.disease and "lung" in slots.disease.lower()
    assert slots.stage in {"II", "2"} or (slots.stage or "").startswith("II")
    assert slots.location and "lower" in slots.location.lower()
    assert "KRAS" in slots.genes
    assert "TP53" in slots.genes


def test_rules_schizophrenia_heart():
    slots = parse_rules("18 year old male with schizophrenia and a heart condition")
    assert slots.disease and "schizophren" in slots.disease.lower()
    assert any("heart" in c.lower() for c in slots.comorbidities)
    assert slots.sex == "male"
    assert slots.age_years == 18


def test_llm_off_uses_rules():
    slots = parse_with_optional_llm(
        "depression with obesity", llm="rules"
    )
    assert slots.parse_method == "rules"
    assert slots.disease


def test_slots_to_biomarker_kwargs():
    slots = QuerySlots(
        disease="depression",
        comorbidities=["obesity"],
        location="brain",
        age_years=24,
        sex="male",
        top=15,
        comorbidity_weight=0.55,
    )
    kw = slots.to_biomarker_kwargs()
    assert kw["disease"] == "depression"
    assert kw["comorbidities"] == ["obesity"]
    assert kw["location"] == "brain"
    assert kw["age_years"] == 24
    assert kw["top_n"] == 15


def test_fill_slots_interactively(monkeypatch):
    answers = iter(
        [
            "depression",  # disease
            "obesity",  # comorbidities
            "brain",  # location
            "24",  # age
            "male",  # sex
            "",  # stage
            "",  # genes
            "",  # smoking
            "0.65",  # comorbidity weight
            "hybrid",  # mode
            "20",  # top
        ]
    )
    monkeypatch.setattr("builtins.input", lambda _: next(answers))
    slots = fill_slots_interactively()
    assert slots.disease == "depression"
    assert slots.comorbidities == ["obesity"]
    assert slots.location == "brain"
    assert slots.age_years == 24
    assert slots.sex == "male"
    assert slots.mode == "hybrid"
    assert slots.top == 20
    assert slots.parse_method == "ask"


def test_nl_cli_show_slots():
    from typer.testing import CliRunner

    from exhalepath.cli import app

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "nl",
            "24yo obese male with depression",
            "--llm",
            "rules",
            "--show-slots",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "depression" in result.output.lower()
    assert "obesity" in result.output.lower()


def test_nl_cli_runs_biomarker():
    from typer.testing import CliRunner

    from exhalepath.cli import app

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "nl",
            "24yo obese male with depression",
            "--llm",
            "rules",
            "--yes",
            "--no-explain",
            "--no-opentargets",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "ExhalePath Biomarker" in result.output or "depression" in result.output.lower()
    assert "VOC" in result.output or "ppb" in result.output.lower()
