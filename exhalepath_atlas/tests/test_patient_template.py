"""Naturalistic patient input → PatientTemplate."""

from __future__ import annotations

from exhalepath.nl import PatientTemplate, parse_patient_template


def test_compact_schizophrenia_vignette():
    tpl = parse_patient_template(
        "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations and paranoia",
        llm="rules",
    )
    assert tpl.primary_disease and "schizophren" in tpl.primary_disease.lower()
    assert tpl.age_years == 35
    assert tpl.sex == "male"
    assert tpl.smoking == "current"
    assert tpl.bmi == 32
    assert "obesity" in [c.lower() for c in tpl.comorbidities]
    assert "olanzapine" in tpl.medications
    assert "hallucinations" in tpl.symptoms
    assert tpl.location == "brain"
    kw = tpl.to_stack_kwargs()
    assert kw["disease"]
    assert kw["description"] and "olanzapine" in kw["description"]


def test_clinic_note_sections():
    note = """
CC: worsening cough and dyspnea
HPI: 62 year old female former smoker with weight loss
PMH: COPD, hypertension
Meds: albuterol, lisinopril
Location: lung
"""
    tpl = parse_patient_template(note, llm="rules")
    assert tpl.chief_complaint
    assert tpl.age_years == 62
    assert tpl.sex == "female"
    assert tpl.smoking == "former"
    assert any("copd" in c.lower() for c in tpl.comorbidities) or (
        tpl.primary_disease and "copd" in tpl.primary_disease.lower()
    )
    assert "albuterol" in tpl.medications or any("albuterol" in m for m in tpl.medications)


def test_json_and_kv_input():
    raw = '{"disease": "depression", "age": 24, "sex": "male", "comorbidities": ["obesity"], "location": "brain"}'
    tpl = parse_patient_template(raw, llm="rules")
    assert tpl.primary_disease and "depress" in tpl.primary_disease.lower()
    assert tpl.age_years == 24
    assert "obesity" in tpl.comorbidities

    kv = """
disease: lung adenocarcinoma
age: 58
sex: female
genes: KRAS, TP53
stage: II
location: left lower lobe
smoking: former
"""
    tpl2 = parse_patient_template(kv, llm="rules")
    assert "lung" in (tpl2.primary_disease or "").lower()
    assert "KRAS" in tpl2.genes
    assert tpl2.stage in {"II", "2"} or (tpl2.stage or "").upper().startswith("II")


def test_patient_cli_show_template():
    from typer.testing import CliRunner

    from exhalepath.cli import app

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "patient",
            "24yo obese male with depression and insomnia on sertraline",
            "--llm",
            "rules",
            "--show-template",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "depression" in result.output.lower()
    assert "sertraline" in result.output.lower()


def test_stack_nl_flag():
    from typer.testing import CliRunner

    from exhalepath.cli import app

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "stack",
            "--nl",
            "18 year old male with schizophrenia and a heart condition",
            "--llm",
            "rules",
            "--show-template",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "schizophren" in result.output.lower()
