"""Interactive questionnaire to fill QuerySlots.

When natural-language parsing leaves fields blank (or the user skips NL),
prompt for disease, comorbidities, location, demographics, stage, genes, etc.
"""

from __future__ import annotations

from typing import Optional, Sequence

from exhalepath.nl.slots import QuerySlots


def _prompt(label: str, default: str = "", required: bool = False) -> str:
    suffix = f" [{default}]" if default else ""
    while True:
        try:
            raw = input(f"{label}{suffix}: ").strip()
        except EOFError:
            raw = ""
        if not raw and default:
            return default
        if raw or not required:
            return raw
        print("  (required)")


def _prompt_yes_no(label: str, default: bool = True) -> bool:
    d = "Y/n" if default else "y/N"
    raw = _prompt(f"{label} ({d})", default="y" if default else "n")
    if not raw:
        return default
    return raw.lower() in {"y", "yes", "1", "true"}


def _csv(raw: str) -> list[str]:
    if not raw:
        return []
    return [p.strip() for p in raw.replace(";", ",").split(",") if p.strip()]


def fill_slots_interactively(
    slots: Optional[QuerySlots] = None,
    *,
    disease_catalog: Optional[Sequence[str]] = None,
) -> QuerySlots:
    """Prompt for missing / confirmable fields. Returns a complete-enough QuerySlots."""
    s = (slots or QuerySlots()).model_copy(deep=True)
    print()
    print("ExhalePath Atlas — clinical query")
    print("Leave blank to skip optional fields. Comma-separate lists.")
    if disease_catalog:
        print(f"(known diseases: {len(disease_catalog)} — try e.g. lung adenocarcinoma, depression)")
    print()

    s.disease = _prompt(
        "Primary disease / condition",
        default=s.disease or "",
        required=True,
    ) or s.disease

    raw = _prompt(
        "Comorbidities (comma-separated)",
        default=", ".join(s.comorbidities) if s.comorbidities else "",
    )
    if raw:
        s.comorbidities = _csv(raw)
    elif not s.comorbidities:
        s.comorbidities = []

    s.location = _prompt(
        "Anatomic location / tissue (e.g. lung, brain, left lower lobe)",
        default=s.location or "",
    ) or s.location

    age_default = str(int(s.age_years)) if s.age_years is not None else ""
    age_raw = _prompt("Age (years)", default=age_default)
    if age_raw:
        try:
            s.age_years = float(age_raw)
        except ValueError:
            print("  (ignored invalid age)")

    sex_raw = _prompt("Sex (male/female/other)", default=s.sex or "")
    if sex_raw:
        sex = sex_raw.lower()
        if sex in {"male", "female", "other"}:
            s.sex = sex  # type: ignore[assignment]
        elif sex in {"m", "man", "boy"}:
            s.sex = "male"
        elif sex in {"f", "woman", "girl"}:
            s.sex = "female"

    stage_raw = _prompt("Stage (I–IV)", default=s.stage or "")
    if stage_raw:
        s.stage = stage_raw.upper()

    genes_raw = _prompt(
        "Driver genes (comma-separated)",
        default=", ".join(s.genes) if s.genes else "",
    )
    if genes_raw:
        s.genes = [g.upper() for g in _csv(genes_raw)]

    smoke_raw = _prompt(
        "Smoking (never/former/current)",
        default=s.smoking or "",
    )
    if smoke_raw:
        sm = smoke_raw.lower()
        if sm in {"never", "former", "current"}:
            s.smoking = sm  # type: ignore[assignment]

    if s.comorbidities:
        w_raw = _prompt(
            "Comorbidity weight 0–1",
            default=str(s.comorbidity_weight),
        )
        try:
            s.comorbidity_weight = max(0.0, min(1.5, float(w_raw)))
        except ValueError:
            pass

    mode_raw = _prompt("Mode (physiology|hybrid|legacy)", default=s.mode or "hybrid")
    if mode_raw and mode_raw.lower() in {"physiology", "hybrid", "legacy"}:
        s.mode = mode_raw.lower()  # type: ignore[assignment]

    top_raw = _prompt("Top-N VOCs", default=str(s.top))
    if top_raw:
        try:
            s.top = int(top_raw)
        except ValueError:
            pass

    s.parse_method = (s.parse_method or "ask") + "+ask" if s.parse_method else "ask"
    return s


def confirm_slots(slots: QuerySlots) -> bool:
    print()
    print("Parsed query:")
    for line in slots.summary_lines():
        print(f"  {line}")
    print()
    return _prompt_yes_no("Run biomarker prediction?", default=True)
