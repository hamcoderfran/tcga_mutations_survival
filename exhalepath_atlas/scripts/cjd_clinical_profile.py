"""Creutzfeldt–Jakob disease clinical VOC profile across disease tempo.

Models a typical sporadic CJD course (incubating → weeks → months) with
increasing affected-cell density/activity. There is no published exhaled-VOC
gold standard for CJD; signals are mechanistic hypotheses (lipid peroxidation,
neuroinflammation, necrosis), analogized from AD/PD/TBI oxidative panels.
"""

from __future__ import annotations

import copy
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from exhalepath.biomarker import ExhaleBiomarkerEngine

OXIDATIVE_PANEL = [
    "hexanal",
    "heptanal",
    "nonanal",
    "octanal",
    "pentane",
    "ethane",
    "acetone",
    "ammonia",
    "hydrogen_sulfide",
    "acetaldehyde",
    "isoprene",
]

# Absolute Δppb floors roughly above typical PTR-MS/GC noise for these species
DELTA_FLOOR = {
    "hexanal": 0.8,
    "heptanal": 0.5,
    "nonanal": 0.8,
    "octanal": 0.4,
    "pentane": 1.5,
    "ethane": 0.4,
    "acetone": 40.0,
    "ammonia": 25.0,
    "hydrogen_sulfide": 0.8,
    "acetaldehyde": 2.0,
    "isoprene": 8.0,  # decrease magnitude
}

HINT_FOLD = 1.20
OBVIOUS_FOLD = 1.60
OBVIOUS_ISOPRENE_FOLD = 0.65


@dataclass
class Phase:
    id: str
    label: str
    weeks_from_onset: str
    clinical: str
    stage: str
    affected_fraction: float
    affected_activity: float
    prior_scale: float = 1.0  # <1 for incubating / very early prior tempering


PHASES = [
    Phase(
        id="incubating",
        label="Incubating / preclinical",
        weeks_from_onset="months–years before onset (silent)",
        clinical=(
            "Prion propagation without frank dementia; MRI/EEG normal; "
            "CSF RT-QuIC may convert late in incubation. Breath VOCs expected near baseline."
        ),
        stage="I",
        affected_fraction=0.02,
        affected_activity=0.7,
        prior_scale=0.15,
    ),
    Phase(
        id="prodromal",
        label="Prodromal / very early clinical",
        weeks_from_onset="0–4 weeks from clinical onset",
        clinical=(
            "Subtle cognitive slowing, anxiety/depression, sleep disturbance; "
            "MRI/EEG often still nondiagnostic; RT-QuIC may already be positive."
        ),
        stage="I",
        affected_fraction=0.10,
        affected_activity=1.0,
        prior_scale=0.55,
    ),
    Phase(
        id="early_clinical",
        label="Early clinical",
        weeks_from_onset="4–8 weeks",
        clinical=(
            "Rapid cognitive decline, ataxia, emerging myoclonus; "
            "diffusion MRI cortical/basal-ganglia ribboning often appears here."
        ),
        stage="II",
        affected_fraction=0.25,
        affected_activity=1.4,
        prior_scale=0.85,
    ),
    Phase(
        id="fulminant",
        label="Fulminant",
        weeks_from_onset="8–16 weeks",
        clinical=(
            "Dementia, startle myoclonus, pyramidal/extrapyramidal signs; "
            "periodic sharp waves on EEG; clinical diagnosis usually clear."
        ),
        stage="III",
        affected_fraction=0.50,
        affected_activity=1.95,
        prior_scale=1.0,
    ),
    Phase(
        id="terminal",
        label="Terminal / akinetic mute",
        weeks_from_onset="16–26 weeks (typical sCJD)",
        clinical=(
            "Akinetic mutism with severe neuronal loss and gliosis; "
            "median sCJD survival ~4–6 months from clinical onset."
        ),
        stage="IV",
        affected_fraction=0.75,
        affected_activity=2.4,
        prior_scale=1.1,
    ),
]


def _scale_disease_priors(disease: dict[str, Any], scale: float) -> dict[str, Any]:
    d = copy.deepcopy(disease)
    d["voc_log2fc_prior"] = {
        k: float(v) * scale for k, v in (d.get("voc_log2fc_prior") or {}).items()
    }
    # Soften pathway bias toward 1.0 when scale < 1
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
            "confidence": p.confidence,
        }
    return out


def _panel_hits(preds: dict[str, Any]) -> dict[str, Any]:
    hints: list[str] = []
    obvious: list[str] = []
    rows = []
    for vid in OXIDATIVE_PANEL:
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
            if fold <= OBVIOUS_ISOPRENE_FOLD and mag >= floor:
                obvious.append(vid)
            elif fold <= (1.0 / HINT_FOLD) and mag >= floor * 0.6:
                hints.append(vid)
        else:
            if (fold >= OBVIOUS_FOLD or fold <= 1.0 / OBVIOUS_FOLD) and abs(delta) >= floor:
                obvious.append(vid)
            elif (fold >= HINT_FOLD or fold <= 1.0 / HINT_FOLD) and abs(delta) >= floor * 0.6:
                hints.append(vid)

    log2s = [
        abs(math.log2(row["fold_change"]))
        for row in rows
        if row["voc_id"] != "isoprene" and row["fold_change"] > 0
    ]
    mean_abs_log2fc = sum(log2s) / len(log2s) if log2s else 0.0

    if len(obvious) >= 3 or mean_abs_log2fc >= 0.7:
        verdict = "obvious"
    elif (len(hints) + len(obvious) >= 2 and mean_abs_log2fc >= 0.18) or mean_abs_log2fc >= 0.28:
        verdict = "meaningful_hint"
    else:
        verdict = "subtle_or_absent"
    return {
        "panel": rows,
        "hint_vocs": hints,
        "obvious_vocs": obvious,
        "mean_abs_log2fc_oxidative": round(mean_abs_log2fc, 4),
        "verdict": verdict,
    }


def run_profile(
    *,
    age: float = 62,
    sex: str = "female",
    genes: list[str] | None = None,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    genes = genes or ["PRNP"]
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    base = engine.kb.resolve_disease("creutzfeldt_jakob")
    out_dir = Path(out_dir or "runs/cjd_clinical_profile")
    out_dir.mkdir(parents=True, exist_ok=True)

    vignette = {
        "disease": "Creutzfeldt-Jakob disease (sporadic)",
        "patient": {
            "age_years": age,
            "sex": sex,
            "genes": genes,
            "phenotype": (
                "Rapidly progressive dementia with myoclonus and cerebellar ataxia; "
                "no cancer comorbidity modeled."
            ),
            "location": "brain (cortex / basal ganglia burden proxy)",
        },
        "caveat": (
            "Research hypothesis only — not a diagnostic test. No peer-reviewed "
            "exhaled-VOC signature is established for CJD; gold-standard workup "
            "remains MRI DWI, EEG, CSF RT-QuIC / 14-3-3 / t-tau. "
            "Model absolute folds are upper-bound estimates; interpret direction and tempo."
        ),
        "thresholds": {
            "hint_fold": HINT_FOLD,
            "obvious_fold": OBVIOUS_FOLD,
            "delta_floor_ppb": DELTA_FLOOR,
        },
    }

    phases_out = []
    first_hint = None
    first_obvious = None
    for phase in PHASES:
        # Temporarily swap scaled disease entry into the knowledge base
        scaled = _scale_disease_priors(base, phase.prior_scale)
        engine.kb.diseases[scaled["disease_id"]] = scaled
        report = engine.predict(
            "creutzfeldt_jakob",
            location="brain",
            top_n=50,
            stage=phase.stage,
            genes=genes,
            affected_fraction=phase.affected_fraction,
            affected_activity=phase.affected_activity,
            mode="hybrid",
            sex=sex,
            age_years=age,
            explain=True,
        )
        preds = _pred_map(report)
        hits = _panel_hits(preds)
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
                "oxidative_panel": hits,
                "mechanisms": mechs,
                "notes": report.notes[:4],
                "model_version": report.model_version,
            }
        )

    # restore full priors
    engine.kb.diseases[base["disease_id"]] = base

    # AD comparison at early_clinical burden
    ad = engine.predict(
        "alzheimer_disease",
        location="brain",
        top_n=20,
        stage="II",
        genes=["APP", "PSEN1"],
        affected_fraction=0.25,
        affected_activity=1.4,
        mode="hybrid",
        sex=sex,
        age_years=age,
        explain=False,
    )
    ad_hits = _panel_hits(_pred_map(ad))
    cjd_early = next(p for p in phases_out if p["id"] == "early_clinical")
    cjd_prod = next(p for p in phases_out if p["id"] == "prodromal")
    cjd_inc = next(p for p in phases_out if p["id"] == "incubating")

    summary = {
        "first_meaningful_hint_phase": first_hint,
        "first_obvious_phase": first_obvious,
        "interpretation": _interpret(first_hint, first_obvious, phases_out),
        "tempo": {
            "incubating_verdict": cjd_inc["oxidative_panel"]["verdict"],
            "prodromal_verdict": cjd_prod["oxidative_panel"]["verdict"],
            "early_clinical_verdict": cjd_early["oxidative_panel"]["verdict"],
            "mean_abs_log2fc": {
                p["id"]: p["oxidative_panel"]["mean_abs_log2fc_oxidative"]
                for p in phases_out
            },
        },
        "vs_alzheimer_same_burden": {
            "cjd_early_clinical_mean_abs_log2fc": cjd_early["oxidative_panel"][
                "mean_abs_log2fc_oxidative"
            ],
            "ad_matched_mean_abs_log2fc": ad_hits["mean_abs_log2fc_oxidative"],
            "cjd_verdict": cjd_early["oxidative_panel"]["verdict"],
            "ad_verdict": ad_hits["verdict"],
        },
    }

    payload = {"vignette": vignette, "phases": phases_out, "summary": summary}
    (out_dir / "cjd_clinical_profile.json").write_text(json.dumps(payload, indent=2))
    (out_dir / "cjd_clinical_profile.md").write_text(_markdown(payload))
    return payload


def _interpret(first_hint, first_obvious, phases) -> str:
    by_id = {p["id"]: p for p in phases}
    lines = []
    inc = by_id.get("incubating")
    if inc:
        lines.append(
            f"Incubating phase stays `{inc['oxidative_panel']['verdict']}` "
            f"(mean |log2FC|={inc['oxidative_panel']['mean_abs_log2fc_oxidative']})."
        )
    if first_hint is None:
        lines.append("No meaningful oxidative VOC hint across modeled clinical phases.")
    else:
        p = by_id[first_hint]
        lines.append(
            f"First meaningful VOC hint: **{p['label']}** ({p['weeks_from_onset']}), "
            f"verdict `{p['oxidative_panel']['verdict']}`."
        )
    if first_obvious:
        p = by_id[first_obvious]
        lines.append(
            f"Becomes obvious on the oxidative breath panel: **{p['label']}** "
            f"({p['weeks_from_onset']})."
        )
    else:
        lines.append("Did not cross the 'obvious' threshold in this run.")
    lines.append(
        "Even when model VOCs move early after clinical onset, they are a "
        "nonspecific oxidative/necrosis signature — not prion-specific. "
        "MRI DWI + CSF RT-QuIC remain the appropriate early diagnostic path."
    )
    return " ".join(lines)


def _markdown(payload: dict[str, Any]) -> str:
    v = payload["vignette"]
    s = payload["summary"]
    lines = [
        "# Creutzfeldt–Jakob disease — clinical VOC profile",
        "",
        f"**Patient:** {v['patient']['age_years']}yo {v['patient']['sex']}, "
        f"genes={','.join(v['patient']['genes'])}, location=brain.",
        "",
        f"> {v['caveat']}",
        "",
        "## Bottom line",
        "",
        s["interpretation"],
        "",
        f"- First meaningful hint phase: `{s['first_meaningful_hint_phase']}`",
        f"- First obvious phase: `{s['first_obvious_phase']}`",
        f"- CJD vs AD (matched early burden) mean |log2FC|: "
        f"{s['vs_alzheimer_same_burden']['cjd_early_clinical_mean_abs_log2fc']} vs "
        f"{s['vs_alzheimer_same_burden']['ad_matched_mean_abs_log2fc']}",
        "",
        "## Phases",
        "",
    ]
    for p in payload["phases"]:
        ox = p["oxidative_panel"]
        lines += [
            f"### {p['label']} ({p['weeks_from_onset']})",
            "",
            p["clinical"],
            "",
            f"- Model burden: stage {p['stage']}, affected_fraction={p['affected_fraction']}, "
            f"activity={p['affected_activity']}, prior_scale={p['prior_scale']}",
            f"- Oxidative panel verdict: **{ox['verdict']}** "
            f"(mean |log2FC|={ox['mean_abs_log2fc_oxidative']})",
            f"- Hint VOCs: {', '.join(ox['hint_vocs']) or '—'}",
            f"- Obvious VOCs: {', '.join(ox['obvious_vocs']) or '—'}",
            "",
            "| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |",
            "|---:|---|---:|---:|---:|---:|",
        ]
        for t in p["top10"]:
            lines.append(
                f"| {t['rank']} | {t['name']} | {t['healthy_ppb']:.2f} | "
                f"{t['predicted_ppb']:.2f} | {t['delta_ppb']:+.2f} | {t['fold_change']:.2f}x |"
            )
        lines.append("")
        if p.get("mechanisms"):
            lines.append("Mechanisms:")
            for m in p["mechanisms"][:4]:
                lines.append(f"- {m}")
            lines.append("")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    profile = run_profile()
    print(profile["summary"]["interpretation"])
    print("tempo:", json.dumps(profile["summary"]["tempo"], indent=2))
    print("Wrote runs/cjd_clinical_profile/")
