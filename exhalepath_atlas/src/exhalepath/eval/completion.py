"""End-to-end dataset completion gate for ExhalePath Atlas."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from ..config import PACKAGE_ROOT
from .audit import run_full_audit
from .multisite import run_multisite_literature_eval
from .public_breath import evaluate_public_breath


def _pytest_ok() -> dict[str, Any]:
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
        "tail": "\n".join(out.strip().splitlines()[-8:]),
    }


def run_completion_suite(*, out_dir: Path | None = None) -> dict[str, Any]:
    """
    Run pytest + audit + multisite + public breath and apply pass gates.

    Gates (research-ready CI bar):
      - pytest all green
      - audit.passed
      - multisite composite_accuracy >= 0.90
      - literature public-breath directional >= 0.90 and elevated_recall@15 >= 0.85
      - sci_data public-breath elevated_directional >= 0.60
      - atlas diseases >= 100, datasources fused >= 12
    """
    out_dir = Path(out_dir or PACKAGE_ROOT / "runs" / "completion")
    out_dir.mkdir(parents=True, exist_ok=True)

    pytest_res = _pytest_ok()
    audit = run_full_audit(out_path=out_dir / "audit_report.json")
    multi = run_multisite_literature_eval(out_dir=out_dir / "multisite_eval", mode="hybrid")
    public = evaluate_public_breath(top_k=15, mode="hybrid", out_dir=out_dir / "public_breath_eval")

    lit_cases = [c for c in public["cases"] if str(c.get("source", "")).startswith("literature")]
    sci_cases = [c for c in public["cases"] if str(c.get("source", "")).startswith("scientific_data")]

    def _mean(vals):
        vals = [v for v in vals if v is not None]
        return float(sum(vals) / len(vals)) if vals else None

    lit_dir = _mean([c["directional_accuracy"] for c in lit_cases])
    lit_rec = _mean([c["elevated_recall_at_k"] for c in lit_cases])
    sci_elev_dir = _mean([c["elevated_directional_accuracy"] for c in sci_cases])
    sci_rec = _mean([c["elevated_recall_at_k"] for c in sci_cases])

    # capability
    cap_path = PACKAGE_ROOT / "data" / "knowledge" / "atlas_capability.json"
    capability = json.loads(cap_path.read_text()) if cap_path.exists() else {}

    gates = {
        "pytest": bool(pytest_res["passed"]),
        "audit": bool(audit.passed),
        "multisite_composite_ge_0.90": float(multi["overall"]["composite_accuracy"]) >= 0.90,
        "lit_public_directional_ge_0.90": (lit_dir or 0) >= 0.90,
        "lit_public_recall_ge_0.85": (lit_rec or 0) >= 0.85,
        "sci_data_elev_directional_ge_0.60": (sci_elev_dir or 0) >= 0.60,
        "diseases_ge_100": int(capability.get("n_diseases") or 0) >= 100,
        "datasources_ge_12": int(capability.get("n_datasources_integrated") or 0) >= 12,
    }
    passed = all(gates.values())

    report = {
        "passed": passed,
        "gates": gates,
        "pytest": pytest_res,
        "audit": {
            "passed": audit.passed,
            "literature_case_pass_rate": audit.literature_accuracy.get("case_pass_rate"),
            "directional_accuracy": audit.literature_accuracy.get("directional_accuracy"),
            "zero_shot_pass_rate": audit.zero_shot.get("pass_rate"),
        },
        "multisite": multi["overall"],
        "public_breath": {
            "overall": public["overall"],
            "literature": {
                "n": len(lit_cases),
                "mean_directional_accuracy": lit_dir,
                "mean_elevated_recall_at_k": lit_rec,
            },
            "scientific_data": {
                "n": len(sci_cases),
                "mean_elevated_directional_accuracy": sci_elev_dir,
                "mean_elevated_recall_at_k": sci_rec,
            },
            "cases": public["cases"],
        },
        "capability": capability,
    }
    (out_dir / "completion_report.json").write_text(json.dumps(report, indent=2))
    return report
