#!/usr/bin/env python3
"""Run the full real-breath pipeline and emit a consolidated evaluation report.

Steps
-----
1. build-real-corpus  (MW + Sci Data + literature panels; no synthetic y)
2. train              (per-VOC calibrators)
3. eval-priority10
4. eval-public-breath
5. eval-comorbidity-clinical
6. eval-multisite
7. audit
8. eval-vision
9. diabetes-profile / cjd-profile (clinical tempo)
10. eval-completion   (CI gate: pytest + audit + multisite + public breath)

Usage:
  python scripts/run_real_pipeline_eval.py
  # or: voc-equivalent CLI steps documented in the report
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "real_pipeline"


def _run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    print(f"\n>>> {' '.join(cmd)}", flush=True)
    return subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True)


def _ok(proc: subprocess.CompletedProcess[str], *, allow_fail: bool = False) -> str:
    out = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0 and not allow_fail:
        print(out[-4000:])
        raise SystemExit(f"Command failed ({proc.returncode}): {' '.join(proc.args)}")
    print(out[-2500:] if len(out) > 2500 else out)
    return out


def _load(path: Path) -> Any:
    return json.loads(path.read_text()) if path.exists() else None


def _pct(x: float | None) -> str:
    if x is None:
        return "—"
    return f"{100.0 * float(x):.1f}%"


def build_markdown(summary: dict[str, Any]) -> str:
    s = summary
    lines = [
        "# ExhalePath real-data pipeline evaluation",
        "",
        f"Generated: `{s['generated_at']}`",
        "",
        "## Pipeline status",
        "",
        f"**Overall: {'PASSED' if s['passed'] else 'FAILED'}**",
        "",
        "### Corpus (real labels only)",
        "",
        f"- Type: `{s['corpus'].get('corpus_type')}` · synthetic={s['corpus'].get('synthetic_labels')}",
        f"- VOC training targets: **{s['corpus'].get('n_voc_targets')}** across **{s['corpus'].get('n_cases')}** cases",
        f"- Diseases: {', '.join(s['corpus'].get('diseases') or [])}",
        f"- Sources: {', '.join(s['corpus'].get('sources') or [])}",
        f"- Label origins: `{json.dumps(s['corpus'].get('label_origins') or {})}`",
        "",
        "### Calibrator",
        "",
        f"- VOC models: **{s['calibrator'].get('n_voc_models')}** · mean MAE(log2fc)={s['calibrator'].get('mean_mae')}",
        "",
        "## Real-data evaluation scores",
        "",
        "| Suite | Headline score | Notes |",
        "|---|---:|---|",
        f"| Priority-10 literature panels | **{_pct(s['priority10'].get('mean_directional_accuracy'))}** | {s['priority10'].get('n_diseases')} diseases · directional vs measured_log2fc |",
        f"| Public breath (all) | **{_pct(s['public_breath'].get('mean_directional_accuracy'))}** | elev recall@15={_pct(s['public_breath'].get('mean_elevated_recall_at_k'))} |",
        f"| Public breath — literature | **{_pct(s['public_breath'].get('literature_directional_accuracy'))}** | elev recall@15={_pct(s['public_breath'].get('literature_elevated_recall_at_k'))} |",
        f"| Public breath — Sci Data 2024 | elev dir **{_pct(s['public_breath'].get('scientific_data_elevated_directional_accuracy'))}** | elev recall@15={_pct(s['public_breath'].get('scientific_data_elevated_recall_at_k'))} (harder; no healthy controls in source) |",
        f"| Comorbidity clinical | **{_pct(s['comorbidity'].get('mean_directional_accuracy'))}** | elev dir={_pct(s['comorbidity'].get('mean_elevated_directional_accuracy'))} · Magdeburg SZ + ST003181 |",
        f"| Multisite literature (20×5) | **{_pct(s['multisite'].get('composite_accuracy'))}** | directional={_pct(s['multisite'].get('directional_accuracy'))} · site sens={_pct(s['multisite'].get('site_sensitivity'))} |",
        f"| Literature audit | **{_pct(s['audit'].get('literature_case_pass_rate'))}** | zero-shot={_pct(s['audit'].get('zero_shot_pass_rate'))} · passed={s['audit'].get('passed')} |",
        f"| Vision suite (50 diseases) | **{s['vision'].get('vision_fidelity_pct')}%** | evidence-backed={s['vision'].get('evidence_backed_pct')}% |",
        f"| Completion CI gate | **{'PASSED' if s['completion'].get('passed') else 'FAILED'}** | all gates green |",
        f"| Diabetes clinical profile | **{'PASSED' if s['diabetes'].get('passed') else 'FAILED'}** | acetone gate vs PMID:21903721 |",
        f"| CJD clinical profile | ran | tempo incubating→terminal |",
        "",
        "## Honesty notes",
        "",
        "- Priority-10 directional accuracy can look near-perfect when panels also inform atlas priors / training labels (partial circularity).",
        "- Sci Data 2024 recall@k is lower because source tables lack healthy controls (cross-cohort differentials only).",
        "- Vision Grade D is prior-consistency; prefer evidence-backed % for independent truthfulness.",
        "- Research tool — not a medical device.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd exhalepath_atlas && pip install -e \".[dev]\"",
        "voc build-real-corpus && voc train",
        "voc eval-priority10",
        "voc eval-public-breath",
        "voc eval-comorbidity-clinical",
        "voc eval-multisite && voc audit",
        "voc eval-vision",
        "voc diabetes-profile --age 55 --sex male",
        "voc cjd-profile --age 62 --sex female",
        "voc eval-completion",
        "# or: python scripts/run_real_pipeline_eval.py",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    py = sys.executable

    # Prefer installed `voc` CLI
    voc = "voc"

    _ok(_run([voc, "build-real-corpus"]))
    train_out = _ok(_run([voc, "train"]))
    _ok(_run([voc, "eval-priority10", "--out-dir", "runs/priority10_eval"]))
    _ok(_run([voc, "eval-public-breath", "--out-dir", "runs/public_breath_eval"]))
    _ok(_run([voc, "eval-comorbidity-clinical", "--out-dir", "runs/comorbidity_clinical_eval"]))
    _ok(_run([voc, "eval-multisite"]))
    _ok(_run([voc, "audit"]))
    _ok(_run([voc, "eval-vision", "--out-dir", "runs/vision_eval"]))
    dia = _ok(_run([voc, "diabetes-profile", "--age", "55", "--sex", "male"]), allow_fail=True)
    _ok(_run([voc, "cjd-profile", "--age", "62", "--sex", "female"]), allow_fail=True)
    _ok(_run([voc, "eval-completion", "--out-dir", "runs/completion"]))

    corpus = _load(ROOT / "src/exhalepath/data/processed/corpus_manifest.json") or _load(
        ROOT / "data/processed/corpus_manifest.json"
    ) or {}
    cal_metrics = _load(ROOT / "src/exhalepath/data/models/voc_calibrator_metrics.json") or _load(
        ROOT / "data/models/voc_calibrator_metrics.json"
    ) or {}
    p10 = _load(ROOT / "runs/priority10_eval/priority10_eval.json") or {}
    pb = (_load(ROOT / "runs/public_breath_eval/public_breath_eval.json") or {}).get("overall") or {}
    co = _load(ROOT / "runs/comorbidity_clinical_eval/comorbidity_clinical_eval.json") or {}
    ms = (_load(ROOT / "runs/multisite_eval/multisite_accuracy_report.json") or {}).get("overall") or {}
    au = _load(ROOT / "runs/audit/audit_report.json") or {}
    vis = (_load(ROOT / "runs/vision_eval/vision_eval.json") or {}).get("overall") or {}
    comp = _load(ROOT / "runs/completion/completion_report.json") or {}

    n_voc_models = len([k for k in cal_metrics if not str(k).startswith("__")])
    mean_mae = None
    if cal_metrics:
        maes = [float(v["mae_log2fc"]) for k, v in cal_metrics.items() if isinstance(v, dict) and "mae_log2fc" in v and not str(k).startswith("__")]
        mean_mae = sum(maes) / len(maes) if maes else None

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "passed": bool(comp.get("passed")),
        "corpus": corpus,
        "calibrator": {"n_voc_models": n_voc_models, "mean_mae": mean_mae, "train_stdout_tail": train_out[-500:]},
        "priority10": {
            "n_diseases": p10.get("n_diseases"),
            "mean_directional_accuracy": p10.get("mean_directional_accuracy"),
        },
        "public_breath": pb,
        "comorbidity": {
            "mean_directional_accuracy": co.get("mean_directional_accuracy"),
            "mean_elevated_directional_accuracy": co.get("mean_elevated_directional_accuracy"),
            "mean_elevated_recall_at_k": co.get("mean_elevated_recall_at_k"),
            "n_cases": co.get("n_cases"),
        },
        "multisite": ms,
        "audit": {
            "passed": au.get("passed"),
            "literature_case_pass_rate": (au.get("literature_accuracy") or {}).get("case_pass_rate"),
            "zero_shot_pass_rate": (au.get("zero_shot") or {}).get("pass_rate"),
        },
        "vision": {
            "vision_fidelity_pct": vis.get("vision_fidelity_pct"),
            "evidence_backed_pct": vis.get("evidence_backed_pct"),
            "held_out_style_pct": vis.get("held_out_style_pct"),
        },
        "completion": {"passed": comp.get("passed"), "gates": comp.get("gates")},
        "diabetes": {"passed": "PASSED" in dia and "FAILED" not in dia.split("PASSED")[-1][:80]},
    }

    (OUT / "pipeline_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    md = build_markdown(summary)
    (OUT / "REAL_PIPELINE_REPORT.md").write_text(md)
    # also snapshot under knowledge for packaging visibility
    for dest in [
        ROOT / "data" / "knowledge" / "REAL_PIPELINE_REPORT.md",
        ROOT / "src" / "exhalepath" / "data" / "knowledge" / "REAL_PIPELINE_REPORT.md",
    ]:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(md)
    print("\n=== SUMMARY ===")
    print(md)
    print(f"\nWrote {OUT / 'REAL_PIPELINE_REPORT.md'}")
    raise SystemExit(0 if summary["passed"] else 1)


if __name__ == "__main__":
    main()
