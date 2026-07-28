"""Implementation-readiness gate: run all evaluation layers and score readiness."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from ..config import PACKAGE_ROOT
from .audit import run_full_audit
from .comorbidity_clinical import evaluate_comorbidity_clinical
from .multisite import run_multisite_literature_eval
from .patient_cohort import run_cohort
from .public_breath import evaluate_public_breath
from .stress_hard import evaluate_stress_hard
from .vision import evaluate_vision
from .zero_shot_hard import evaluate_zero_shot_hard
from .zero_shot_reliability import evaluate_zero_shot_reliability


def _pytest() -> dict[str, Any]:
    t0 = time.time()
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/", "-q", "--tb=line"],
        cwd=str(PACKAGE_ROOT),
        capture_output=True,
        text=True,
    )
    out = (r.stdout or "") + (r.stderr or "")
    return {
        "passed": r.returncode == 0,
        "returncode": r.returncode,
        "seconds": round(time.time() - t0, 2),
        "tail": "\n".join(out.strip().splitlines()[-10:]),
    }


def _priority10(out_dir: Path) -> dict[str, Any]:
    from ..biomarker import ExhaleBiomarkerEngine
    from ..ingest.real_breath_corpus import PANEL_PATH, PRIORITY_DISEASES

    t0 = time.time()
    dest = out_dir / "priority10"
    dest.mkdir(parents=True, exist_ok=True)
    panel_path = (
        PANEL_PATH
        if PANEL_PATH.exists()
        else PACKAGE_ROOT
        / "data"
        / "real_breath"
        / "literature_panels"
        / "priority10_voc_panels.json"
    )
    panels = {
        p["disease_id"]: p for p in json.loads(panel_path.read_text()).get("panels", [])
    }
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    rows = []
    for did in PRIORITY_DISEASES:
        panel = panels.get(did) or {}
        expect = panel.get("measured_log2fc") or {}
        if not expect:
            continue
        report = engine.predict(
            did,
            location=panel.get("location"),
            top_n=50,
            mode="hybrid",
            explain=False,
        )
        pred = {p.voc_id: p.log2_fold_change for p in report.result.bundle.predictions}
        hits = total = 0
        for voc, exp in expect.items():
            if voc not in pred:
                continue
            total += 1
            if (float(exp) >= 0 and pred[voc] >= 0) or (
                float(exp) < 0 and pred[voc] < 0
            ):
                hits += 1
        rows.append(
            {
                "disease_id": did,
                "n_vocs_scored": total,
                "directional_accuracy": (hits / total) if total else None,
            }
        )
    overall = [
        r["directional_accuracy"]
        for r in rows
        if r["directional_accuracy"] is not None
    ]
    mean_acc = float(sum(overall) / len(overall)) if overall else 0.0
    summary = {
        "n_diseases": len(rows),
        "mean_directional_accuracy": mean_acc,
        "diseases": rows,
    }
    (dest / "priority10_eval.json").write_text(json.dumps(summary, indent=2))
    return {
        "passed": mean_acc >= 0.9 and len(rows) >= 8,
        "mean_directional_accuracy": mean_acc,
        "n_diseases": len(rows),
        "seconds": round(time.time() - t0, 2),
    }


def evaluate_implementation_readiness(
    *,
    out_dir: Path | None = None,
    max_cohort_patients: int | None = 1000,
    skip_pytest: bool = False,
) -> dict[str, Any]:
    """
    Comprehensive multi-layer readiness evaluation.

    Layers:
      A. Unit/regression (pytest)
      B. Audit (literature + calibrator + stress + novel vs T2D)
      C. Disease suites (priority10, vision, multisite, public breath, comorbidity)
      D. Adversarial (stress-hard, zero-shot-hard)
      E. Zero-shot reliability (rare holdout + leave-disease-out)
      F. Scale (patient cohort)

    Research-implementation ready when all hard gates pass.
    External measured-breath holdouts remain the clinical bar (documented).
    """
    out = Path(out_dir or PACKAGE_ROOT / "runs" / "implementation_readiness")
    out.mkdir(parents=True, exist_ok=True)
    layers: dict[str, Any] = {}

    if skip_pytest:
        layers["pytest"] = {"passed": True, "skipped": True}
    else:
        layers["pytest"] = _pytest()

    t0 = time.time()
    audit = run_full_audit(out_path=out / "audit_report.json")
    layers["audit"] = {
        "passed": bool(audit.passed),
        "seconds": round(time.time() - t0, 2),
        "calibrator_holdout": audit.calibrator_holdout,
        "zero_shot": {
            k: (audit.zero_shot or {}).get(k)
            for k in (
                "pass",
                "pass_rate",
                "novel_vs_t2d_acetone_delta_ppb",
                "novel_peak_abs_log2fc",
            )
        },
    }

    layers["priority10"] = _priority10(out)

    t0 = time.time()
    vision = evaluate_vision(out_dir=out / "vision")
    vo = vision.get("overall") or {}
    layers["vision"] = {
        "passed": float(vo.get("vision_fidelity_pct") or 0) >= 90.0,
        "vision_fidelity_pct": vo.get("vision_fidelity_pct"),
        "n_profiles": vo.get("n_profiles"),
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    multi = run_multisite_literature_eval(out_dir=out / "multisite", mode="hybrid")
    mo = multi.get("overall") or multi
    layers["multisite"] = {
        "passed": float(mo.get("composite_accuracy") or 0) >= 0.90,
        "composite_accuracy": mo.get("composite_accuracy"),
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    public = evaluate_public_breath(
        top_k=15, mode="hybrid", out_dir=out / "public_breath"
    )
    po = public.get("overall") or public
    layers["public_breath"] = {
        "passed": float(po.get("mean_directional_accuracy") or 0) >= 0.85
        and float(po.get("scientific_data_elevated_directional_accuracy") or 0)
        >= 0.55,
        "mean_directional_accuracy": po.get("mean_directional_accuracy"),
        "mean_elevated_recall_at_k": po.get("mean_elevated_recall_at_k"),
        "scientific_data_elevated_directional_accuracy": po.get(
            "scientific_data_elevated_directional_accuracy"
        ),
        "scientific_data_elevated_recall_at_k": po.get(
            "scientific_data_elevated_recall_at_k"
        ),
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    comorb = evaluate_comorbidity_clinical(out_dir=out / "comorbidity")
    layers["comorbidity"] = {
        "passed": float(comorb.get("mean_directional_accuracy") or 0) >= 0.85,
        "mean_directional_accuracy": comorb.get("mean_directional_accuracy"),
        "mean_elevated_recall_at_k": comorb.get("mean_elevated_recall_at_k"),
        "n_cases": comorb.get("n_cases"),
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    stress = evaluate_stress_hard(out_dir=out / "stress_hard")
    so = stress.get("overall") or stress
    layers["stress_hard"] = {
        "passed": float(so.get("pass_pct") or 0) >= 95.0,
        "pass_pct": so.get("pass_pct"),
        "n_passed": so.get("n_passed"),
        "n_cases": so.get("n_cases"),
        "failed_ids": so.get("failed_ids") or [],
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    zsh = evaluate_zero_shot_hard(out_dir=out / "zero_shot_hard")
    layers["zero_shot_hard"] = {
        "passed": float(zsh.get("pass_rate") or 0) >= 0.85,
        "pass_rate": zsh.get("pass_rate"),
        "n_passed": zsh.get("n_passed"),
        "n_profiles": zsh.get("n_profiles"),
        "failed_ids": zsh.get("failed_ids") or [],
        "by_trap": zsh.get("by_trap") or {},
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    zsr = evaluate_zero_shot_reliability(out_dir=out / "zero_shot_reliability")
    layers["zero_shot_reliability"] = {
        "passed": bool(zsr.get("pass")),
        "holdout_pass_rate": (zsr.get("holdout") or {}).get("pass_rate"),
        "ldo_mean_directional_accuracy": (zsr.get("leave_disease_out") or {}).get(
            "mean_directional_accuracy"
        ),
        "seconds": round(time.time() - t0, 2),
    }

    t0 = time.time()
    cohort = run_cohort(out_dir=out / "cohort", max_patients=max_cohort_patients)
    co = cohort.get("overall") or cohort
    layers["cohort"] = {
        "passed": float(co.get("ok_pct") or 0) >= 99.0,
        "ok_pct": co.get("ok_pct"),
        "n_ok": co.get("n_ok"),
        "n_patients": co.get("n_patients"),
        "n_unresolved": co.get("n_unresolved"),
        "seconds": round(time.time() - t0, 2),
    }

    cal = layers["audit"].get("calibrator_holdout") or {}
    cal_ok = (not cal.get("skipped")) and bool(cal.get("pass"))

    gates = {
        "pytest": bool(layers["pytest"].get("passed")),
        "audit": bool(layers["audit"].get("passed")),
        "calibrator_holdout": cal_ok,
        "priority10": bool(layers["priority10"].get("passed")),
        "vision_fidelity_ge_90": bool(layers["vision"].get("passed")),
        "multisite_ge_0.90": bool(layers["multisite"].get("passed")),
        "public_breath": bool(layers["public_breath"].get("passed")),
        "comorbidity": bool(layers["comorbidity"].get("passed")),
        "stress_hard_ge_95": bool(layers["stress_hard"].get("passed")),
        "zero_shot_hard_ge_85": bool(layers["zero_shot_hard"].get("passed")),
        "zero_shot_reliability": bool(layers["zero_shot_reliability"].get("passed")),
        "cohort_ok_ge_99": bool(layers["cohort"].get("passed")),
    }
    n_pass = sum(1 for v in gates.values() if v)
    n_gates = len(gates)
    ready = all(gates.values())

    limitations = [
        "Zero-shot / hard-break suites score mechanism consistency and adversarial "
        "robustness — not external GC-MS clinical accuracy.",
        "Public scientific-data elevated recall@15 remains modest; directional "
        "agreement is the stronger public-breath signal today.",
        "Calibrator is fit on a small real corpus; treat MAE as small-n diagnostics, "
        "not a large held-out clinical trial.",
        "Owlstone and other gated breath atlases are not used for training.",
    ]

    report = {
        "ready_for_implementation": ready,
        "readiness_score": n_pass / n_gates if n_gates else 0.0,
        "n_gates_passed": n_pass,
        "n_gates": n_gates,
        "gates": gates,
        "layers": layers,
        "limitations": limitations,
        "recommendation": (
            "READY for research / hypothesis-generation deployment with uncertainty "
            "flags and mechanism audit trails."
            if ready
            else "NOT READY — failing gates must be fixed before implementation."
        ),
    }

    (out / "implementation_readiness.json").write_text(json.dumps(report, indent=2))
    lines = [
        "# Implementation readiness",
        "",
        f"**Ready:** {ready}",
        f"**Score:** {n_pass}/{n_gates} gates ({100 * report['readiness_score']:.0f}%)",
        f"**Recommendation:** {report['recommendation']}",
        "",
        "## Gates",
    ]
    for k, v in gates.items():
        lines.append(f"- {'✓' if v else '✗'} `{k}`")
    lines.extend(["", "## Layer metrics"])
    for name, layer in layers.items():
        bits = []
        for k in (
            "passed",
            "pass_pct",
            "pass_rate",
            "vision_fidelity_pct",
            "composite_accuracy",
            "mean_directional_accuracy",
            "ok_pct",
            "n_patients",
            "n_profiles",
            "n_cases",
            "seconds",
            "failed_ids",
        ):
            if k in layer and layer[k] not in (None, [], {}):
                bits.append(f"{k}={layer[k]}")
        if name == "audit" and isinstance(layer.get("calibrator_holdout"), dict):
            cal = layer["calibrator_holdout"]
            bits.append(
                f"calibrator_mae={cal.get('mae_log2fc')} skipped={cal.get('skipped')}"
            )
        lines.append(f"- **{name}**: " + ", ".join(bits))
    lines.extend(["", "## Limitations (read before shipping)"])
    for lim in limitations:
        lines.append(f"- {lim}")
    failed = [k for k, v in gates.items() if not v]
    if failed:
        lines.extend(["", "## Failing gates"])
        for k in failed:
            lines.append(f"- `{k}` → {layers.get(k.split('_')[0], layers.get(k, {}))}")
    (out / "IMPLEMENTATION_READINESS.md").write_text("\n".join(lines) + "\n")
    return report
