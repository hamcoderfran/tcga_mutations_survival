"""Type 2 diabetes clinical VOC profile across metabolic control states.

Phases model prediabetes → controlled T2D → poorly controlled → ketotic
exacerbation, with optional obesity comorbidity. Evaluation uses published
breath acetone / ketone literature (PMID:21903721, J Chromatogr B 2013).
"""

from __future__ import annotations

import copy
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from exhalepath.biomarker import ExhaleBiomarkerEngine

# Classic diabetes breath ketone / alcohol panel
KETONE_PANEL = [
    "acetone",
    "isopropanol",
    "2_butanone",
    "2_pentanone",
    "ethanol",
    "isoprene",
    "ammonia",
    "trimethylamine",
    "pentane",
    "hexanal",
]

DELTA_FLOOR = {
    "acetone": 40.0,
    "isopropanol": 2.0,
    "2_butanone": 0.3,
    "2_pentanone": 0.25,
    "ethanol": 10.0,
    "isoprene": 8.0,
    "ammonia": 20.0,
    "trimethylamine": 0.3,
    "pentane": 1.0,
    "hexanal": 0.5,
}

HINT_FOLD = 1.20
OBVIOUS_FOLD = 1.60
# Literature gates (public_breath / literature_benchmarks)
LIT_ACETONE_MIN_FOLD = 1.8
LIT_EXPECT_ELEVATED = ["acetone", "isopropanol", "2_butanone"]


@dataclass
class Phase:
    id: str
    label: str
    clinical_context: str
    clinical: str
    stage: str
    affected_fraction: float
    affected_activity: float
    prior_scale: float = 1.0
    comorbidities: list[str] = field(default_factory=list)
    comorbidity_weight: float = 0.55
    disease_id: str = "type_2_diabetes"


PHASES = [
    Phase(
        id="prediabetes",
        label="Prediabetes / insulin resistance",
        clinical_context="HbA1c ~5.7–6.4%; fasting glucose elevated; no ketosis",
        clinical=(
            "Insulin resistance with mild fatty-acid oxidation stress; "
            "breath acetone may be subtly up vs healthy."
        ),
        stage="I",
        affected_fraction=0.12,
        affected_activity=1.0,
        prior_scale=0.35,
    ),
    Phase(
        id="controlled_t2d",
        label="Controlled T2D",
        clinical_context="HbA1c ~6.5–7.5% on diet/metformin; euglycemic most days",
        clinical=(
            "Established type 2 diabetes with reasonable control; "
            "ketone-body VOC axis mildly activated."
        ),
        stage="II",
        affected_fraction=0.25,
        affected_activity=1.25,
        prior_scale=0.65,
    ),
    Phase(
        id="poorly_controlled",
        label="Poorly controlled T2D",
        clinical_context="HbA1c ≥8.5%; hyperglycemia; possible mild ketonemia",
        clinical=(
            "Chronic hyperglycemia with stronger ketone-body and FAO pathway bias; "
            "acetone / 2-butanone / isopropanol expected clearly elevated."
        ),
        stage="III",
        affected_fraction=0.45,
        affected_activity=1.7,
        prior_scale=1.0,
    ),
    Phase(
        id="t2d_obesity",
        label="T2D + obesity comorbidity",
        clinical_context="Poorly controlled T2D with BMI ≥30; adipose + hepatic burden",
        clinical=(
            "Comorbid obesity fuses adipose lipolysis / FAO priors onto the T2D "
            "ketone axis (common clinical phenotype)."
        ),
        stage="III",
        affected_fraction=0.50,
        affected_activity=1.8,
        prior_scale=1.0,
        comorbidities=["obesity"],
        comorbidity_weight=0.65,
    ),
    Phase(
        id="ketotic",
        label="Ketotic exacerbation / DKA-adjacent",
        clinical_context="Marked ketosis (illness, SGLT2i, or T1D-like DKA)",
        clinical=(
            "Strong ketone-body metabolism activation; breath acetone is the "
            "classic marker (literature folds often ≫2× healthy)."
        ),
        stage="IV",
        affected_fraction=0.70,
        affected_activity=2.2,
        prior_scale=1.25,
    ),
    Phase(
        id="type1_reference",
        label="Type 1 diabetes (reference)",
        clinical_context="Absolute insulin deficiency; high ketosis risk",
        clinical=(
            "T1D atlas entry for comparison — acetone / ketone alcohols elevated; "
            "not the primary T2D vignette."
        ),
        stage="III",
        affected_fraction=0.55,
        affected_activity=1.9,
        prior_scale=1.0,
        disease_id="type1_diabetes",
    ),
]


def _scale_disease_priors(disease: dict[str, Any], scale: float) -> dict[str, Any]:
    d = copy.deepcopy(disease)
    d["voc_log2fc_prior"] = {
        k: float(v) * scale for k, v in (d.get("voc_log2fc_prior") or {}).items()
    }
    pb = {}
    for k, v in (d.get("pathway_bias") or {}).items():
        pb[k] = 1.0 + (float(v) - 1.0) * scale
    d["pathway_bias"] = pb
    return d


def _pred_map(report) -> dict[str, Any]:
    out = {}
    for p in report.result.bundle.predictions:
        out[p.voc_id] = {
            "name": p.name,
            "healthy_ppb": p.healthy_ppb,
            "predicted_ppb": p.predicted_ppb,
            "delta_ppb": p.delta_ppb,
            "fold_change": p.fold_change,
            "log2_fold_change": p.log2_fold_change,
            "confidence": p.confidence,
        }
    return out


def _panel_hits(preds: dict[str, Any]) -> dict[str, Any]:
    hints: list[str] = []
    obvious: list[str] = []
    rows = []
    for vid in KETONE_PANEL:
        if vid not in preds:
            continue
        r = preds[vid]
        fold = float(r["fold_change"])
        delta = float(r["delta_ppb"])
        floor = DELTA_FLOOR.get(vid, 0.5)
        rows.append(
            {
                "voc_id": vid,
                "name": r["name"],
                "healthy_ppb": r["healthy_ppb"],
                "predicted_ppb": r["predicted_ppb"],
                "delta_ppb": delta,
                "fold_change": fold,
            }
        )
        if vid == "isoprene":
            mag = abs(delta)
            if fold <= 0.70 and mag >= floor:
                obvious.append(vid)
            elif fold <= (1.0 / HINT_FOLD) and mag >= floor * 0.6:
                hints.append(vid)
        else:
            if fold >= OBVIOUS_FOLD and abs(delta) >= floor:
                obvious.append(vid)
            elif fold >= HINT_FOLD and abs(delta) >= floor * 0.6:
                hints.append(vid)

    ketone_folds = [
        abs(math.log2(row["fold_change"]))
        for row in rows
        if row["voc_id"] in {"acetone", "isopropanol", "2_butanone", "2_pentanone"}
        and row["fold_change"] > 0
    ]
    mean_ketone = sum(ketone_folds) / len(ketone_folds) if ketone_folds else 0.0
    acetone_fold = next(
        (row["fold_change"] for row in rows if row["voc_id"] == "acetone"), None
    )

    if (
        (acetone_fold is not None and acetone_fold >= LIT_ACETONE_MIN_FOLD and "acetone" in obvious)
        or len(obvious) >= 3
        or mean_ketone >= 0.7
    ):
        verdict = "obvious"
    elif len(hints) + len(obvious) >= 2 or mean_ketone >= 0.25:
        verdict = "meaningful_hint"
    else:
        verdict = "subtle_or_absent"

    return {
        "panel": rows,
        "hint_vocs": hints,
        "obvious_vocs": obvious,
        "mean_abs_log2fc_ketone": round(mean_ketone, 4),
        "acetone_fold": acetone_fold,
        "verdict": verdict,
    }


def _literature_eval(preds: dict[str, Any]) -> dict[str, Any]:
    """Score against public breath / literature diabetes expectations."""
    checks = []
    for voc in LIT_EXPECT_ELEVATED:
        if voc not in preds:
            checks.append({"voc_id": voc, "status": "missing"})
            continue
        fold = float(preds[voc]["fold_change"])
        ok = fold >= (LIT_ACETONE_MIN_FOLD if voc == "acetone" else HINT_FOLD)
        checks.append(
            {
                "voc_id": voc,
                "fold_change": fold,
                "threshold": LIT_ACETONE_MIN_FOLD if voc == "acetone" else HINT_FOLD,
                "passed": ok,
                "status": "pass" if ok else "fail",
            }
        )
    acetone = preds.get("acetone")
    acetone_pass = bool(
        acetone and float(acetone["fold_change"]) >= LIT_ACETONE_MIN_FOLD
    )
    directional = [
        c for c in checks if c.get("status") in {"pass", "fail"}
    ]
    dir_acc = (
        sum(1 for c in directional if c["passed"]) / len(directional)
        if directional
        else None
    )
    return {
        "refs": [
            "PMID:21903721",
            "10.1016/j.jchromb.2012.12.008",
            "public_breath_benchmarks:lit_t2d_acetone",
            "literature_benchmarks:t2d_zero_shot_metabolic",
        ],
        "expect_elevated": LIT_EXPECT_ELEVATED,
        "acetone_min_fold": LIT_ACETONE_MIN_FOLD,
        "acetone_pass": acetone_pass,
        "checks": checks,
        "directional_accuracy": dir_acc,
        "passed": bool(acetone_pass and dir_acc is not None and dir_acc >= 0.67),
    }


def run_profile(
    *,
    age: float = 55,
    sex: str = "male",
    genes: list[str] | None = None,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    genes = genes or ["TCF7L2", "PPARG"]
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    out_dir = Path(out_dir or "runs/diabetes_clinical_profile")
    out_dir.mkdir(parents=True, exist_ok=True)

    vignette = {
        "disease": "Type 2 diabetes mellitus",
        "patient": {
            "age_years": age,
            "sex": sex,
            "genes": genes,
            "phenotype": (
                "Adult-onset hyperglycemia with progressive metabolic control "
                "states; optional obesity comorbidity; T1D reference arm."
            ),
            "location": "systemic (hepatic / adipose / muscle ketone axis)",
        },
        "caveat": (
            "Research hypothesis only — not a diagnostic or glycemic monitor. "
            "Breath acetone tracks ketosis / metabolic stress, not HbA1c directly. "
            "Absolute ppb are model estimates."
        ),
        "literature_gates": {
            "acetone_min_fold": LIT_ACETONE_MIN_FOLD,
            "expect_elevated": LIT_EXPECT_ELEVATED,
        },
    }

    # Cache base disease objects
    base_by_id = {
        "type_2_diabetes": engine.kb.resolve_disease("type_2_diabetes"),
        "type1_diabetes": engine.kb.resolve_disease("type1_diabetes"),
    }

    phases_out = []
    first_hint = None
    first_obvious = None
    for phase in PHASES:
        base = base_by_id[phase.disease_id]
        scaled = _scale_disease_priors(base, phase.prior_scale)
        engine.kb.diseases[scaled["disease_id"]] = scaled
        report = engine.predict(
            phase.disease_id,
            location="systemic",
            top_n=50,
            stage=phase.stage,
            genes=genes if phase.disease_id == "type_2_diabetes" else ["INS", "HLA-DQB1"],
            affected_fraction=phase.affected_fraction,
            affected_activity=phase.affected_activity,
            mode="hybrid",
            sex=sex,
            age_years=age,
            explain=True,
            comorbidities=phase.comorbidities or None,
            comorbidity_weight=phase.comorbidity_weight,
        )
        preds = _pred_map(report)
        hits = _panel_hits(preds)
        lit = _literature_eval(preds)
        top10 = [
            {
                "rank": i,
                "voc_id": p.voc_id,
                "name": p.name,
                "delta_ppb": p.delta_ppb,
                "fold_change": p.fold_change,
                "healthy_ppb": p.healthy_ppb,
                "predicted_ppb": p.predicted_ppb,
            }
            for i, p in enumerate(report.top_vocs[:10], 1)
        ]
        mechs = [m.get("why") for m in (report.mechanisms or [])[:6]]
        if hits["verdict"] in {"meaningful_hint", "obvious"} and first_hint is None:
            first_hint = phase.id
        if hits["verdict"] == "obvious" and first_obvious is None:
            first_obvious = phase.id
        phases_out.append(
            {
                **asdict(phase),
                "top10": top10,
                "ketone_panel": hits,
                "literature_eval": lit,
                "mechanisms": mechs,
                "notes": report.notes[:5],
                "model_version": report.model_version,
            }
        )

    # restore
    for did, d in base_by_id.items():
        engine.kb.diseases[did] = d

    # Primary evaluation phase: poorly controlled T2D
    primary = next(p for p in phases_out if p["id"] == "poorly_controlled")
    obese = next(p for p in phases_out if p["id"] == "t2d_obesity")
    pre = next(p for p in phases_out if p["id"] == "prediabetes")
    ket = next(p for p in phases_out if p["id"] == "ketotic")

    summary = {
        "first_meaningful_hint_phase": first_hint,
        "first_obvious_phase": first_obvious,
        "primary_literature_eval": primary["literature_eval"],
        "obesity_comorbidity_literature_eval": obese["literature_eval"],
        "tempo": {
            p["id"]: {
                "verdict": p["ketone_panel"]["verdict"],
                "acetone_fold": p["ketone_panel"]["acetone_fold"],
                "mean_abs_log2fc_ketone": p["ketone_panel"]["mean_abs_log2fc_ketone"],
                "literature_passed": p["literature_eval"]["passed"],
            }
            for p in phases_out
        },
        "interpretation": _interpret(
            first_hint, first_obvious, pre, primary, obese, ket
        ),
        "overall_passed": bool(
            primary["literature_eval"]["passed"]
            and (primary["ketone_panel"]["acetone_fold"] or 0) >= LIT_ACETONE_MIN_FOLD
        ),
    }

    payload = {"vignette": vignette, "phases": phases_out, "summary": summary}
    (out_dir / "diabetes_clinical_profile.json").write_text(json.dumps(payload, indent=2))
    (out_dir / "diabetes_clinical_profile.md").write_text(_markdown(payload))
    return payload


def _interpret(first_hint, first_obvious, pre, primary, obese, ket) -> str:
    lines = []
    lines.append(
        f"Prediabetes acetone fold={pre['ketone_panel']['acetone_fold']:.2f}× "
        f"({pre['ketone_panel']['verdict']})."
    )
    if first_hint:
        lines.append(f"First meaningful ketone VOC hint: **{first_hint}**.")
    if first_obvious:
        lines.append(f"Becomes obvious: **{first_obvious}**.")
    lit = primary["literature_eval"]
    lines.append(
        f"Poorly controlled T2D literature gate: "
        f"{'PASS' if lit['passed'] else 'FAIL'} "
        f"(acetone {primary['ketone_panel']['acetone_fold']:.2f}× vs ≥{LIT_ACETONE_MIN_FOLD}×; "
        f"dir_acc={lit['directional_accuracy']})."
    )
    lines.append(
        f"With obesity comorbidity acetone={obese['ketone_panel']['acetone_fold']:.2f}×; "
        f"ketotic exacerbation acetone={ket['ketone_panel']['acetone_fold']:.2f}×."
    )
    lines.append(
        "Breath acetone reflects ketosis / FAO stress — complementary to, not a "
        "replacement for, glucose/HbA1c."
    )
    return " ".join(lines)


def _markdown(payload: dict[str, Any]) -> str:
    v = payload["vignette"]
    s = payload["summary"]
    lines = [
        "# Type 2 diabetes — clinical VOC profile",
        "",
        f"**Patient:** {v['patient']['age_years']}yo {v['patient']['sex']}, "
        f"genes={','.join(v['patient']['genes'])}, location=systemic.",
        "",
        f"> {v['caveat']}",
        "",
        "## Bottom line",
        "",
        s["interpretation"],
        "",
        f"- Overall literature gate: "
        f"**{'PASSED' if s['overall_passed'] else 'FAILED'}**",
        f"- First meaningful hint: `{s['first_meaningful_hint_phase']}`",
        f"- First obvious: `{s['first_obvious_phase']}`",
        "",
        "## Phase tempo",
        "",
        "| Phase | Verdict | Acetone fold | Lit pass |",
        "|---|---|---:|:---:|",
    ]
    for pid, t in s["tempo"].items():
        af = t["acetone_fold"]
        af_s = f"{af:.2f}×" if af is not None else "—"
        lines.append(
            f"| {pid} | {t['verdict']} | {af_s} | "
            f"{'yes' if t['literature_passed'] else 'no'} |"
        )
    lines += ["", "## Phases", ""]
    for p in payload["phases"]:
        kp = p["ketone_panel"]
        lit = p["literature_eval"]
        lines += [
            f"### {p['label']}",
            "",
            f"*{p['clinical_context']}*",
            "",
            p["clinical"],
            "",
            f"- Burden: stage {p['stage']}, frac={p['affected_fraction']}, "
            f"activity={p['affected_activity']}, prior_scale={p['prior_scale']}",
            f"- Comorbidities: {', '.join(p['comorbidities']) or '—'}",
            f"- Ketone panel: **{kp['verdict']}** "
            f"(acetone {kp['acetone_fold']:.2f}×, "
            f"mean |log2FC| ketone={kp['mean_abs_log2fc_ketone']})",
            f"- Literature eval: "
            f"{'PASS' if lit['passed'] else 'FAIL'} "
            f"(dir_acc={lit['directional_accuracy']})",
            "",
            "| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |",
            "|---:|---|---:|---:|---:|---:|",
        ]
        for t in p["top10"]:
            lines.append(
                f"| {t['rank']} | {t['name']} | {t['healthy_ppb']:.2f} | "
                f"{t['predicted_ppb']:.2f} | {t['delta_ppb']:+.2f} | "
                f"{t['fold_change']:.2f}x |"
            )
        lines.append("")
        if p.get("mechanisms"):
            lines.append("Mechanisms:")
            for m in p["mechanisms"][:3]:
                lines.append(f"- {m}")
            lines.append("")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    profile = run_profile()
    print(profile["summary"]["interpretation"])
    print("overall_passed:", profile["summary"]["overall_passed"])
    print("tempo:", json.dumps(profile["summary"]["tempo"], indent=2))
    print("Wrote runs/diabetes_clinical_profile/")
