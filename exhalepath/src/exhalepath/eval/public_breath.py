"""Evaluate ExhalePath against public breath VOC benchmarks (top-k / directional)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR


def _load_cases(path: Path | None = None) -> list[dict[str, Any]]:
    p = Path(path or KNOWLEDGE_DIR / "public_breath_benchmarks.json")
    if not p.exists():
        raise FileNotFoundError(
            f"Missing {p}. Run: python -m exhalepath harvest-public-breath"
        )
    return list(json.loads(p.read_text())["cases"])


def evaluate_public_breath(
    *,
    top_k: int = 15,
    mode: str = "hybrid",
    out_dir: Path | None = None,
) -> dict[str, Any]:
    """
    Score each public case:
      - elevated_recall@k: fraction of expected-elevated VOCs in predicted top-k by |Δppb|
      - elevated_precision-like: among expected set, fraction correctly predicted fold>1
      - directional_accuracy on union of elevated+suppressed
    """
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    cases_out = []
    for case in _load_cases():
        report = engine.predict(
            case.get("disease") or case.get("disease_id"),
            location=case.get("location"),
            top_n=50,
            mode=mode,
        )
        ranked_ids = [p.voc_id for p in report.top_vocs]
        top_ids = set(ranked_ids[:top_k])
        by_id = {p.voc_id: p for p in report.top_vocs}
        elev = list(case.get("expect_elevated") or [])
        supp = list(case.get("expect_suppressed") or [])

        elev_in_top = [v for v in elev if v in top_ids]
        elev_dir_hits = 0
        elev_dir_tot = 0
        for v in elev:
            if v not in by_id:
                continue
            elev_dir_tot += 1
            if by_id[v].fold_change > 1.0:
                elev_dir_hits += 1
        supp_dir_hits = 0
        supp_dir_tot = 0
        for v in supp:
            if v not in by_id:
                continue
            supp_dir_tot += 1
            if by_id[v].fold_change < 1.0:
                supp_dir_hits += 1

        dir_tot = elev_dir_tot + supp_dir_tot
        dir_hits = elev_dir_hits + supp_dir_hits
        row = {
            "case_id": case["case_id"],
            "disease": case.get("disease"),
            "disease_id": report.disease_id,
            "source": case.get("source"),
            "n_expect_elevated": len(elev),
            "elevated_recall_at_k": (len(elev_in_top) / len(elev)) if elev else None,
            "elevated_in_top_k": elev_in_top,
            "elevated_directional_accuracy": (elev_dir_hits / elev_dir_tot)
            if elev_dir_tot
            else None,
            "suppressed_directional_accuracy": (supp_dir_hits / supp_dir_tot)
            if supp_dir_tot
            else None,
            "directional_accuracy": (dir_hits / dir_tot) if dir_tot else None,
            "top_k": top_k,
            "predicted_top": ranked_ids[:top_k],
        }
        cases_out.append(row)

    recalls = [c["elevated_recall_at_k"] for c in cases_out if c["elevated_recall_at_k"] is not None]
    dirs = [c["directional_accuracy"] for c in cases_out if c["directional_accuracy"] is not None]
    overall = {
        "n_cases": len(cases_out),
        "mean_elevated_recall_at_k": float(sum(recalls) / len(recalls)) if recalls else None,
        "mean_directional_accuracy": float(sum(dirs) / len(dirs)) if dirs else None,
        "top_k": top_k,
        "mode": mode,
    }
    report = {"overall": overall, "cases": cases_out}
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "public_breath_eval.json").write_text(json.dumps(report, indent=2))
    return report
