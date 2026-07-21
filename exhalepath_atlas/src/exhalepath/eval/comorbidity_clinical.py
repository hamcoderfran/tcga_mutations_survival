"""Evaluate comorbidity-aware predictions against clinical breath/metabolome benchmarks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR


def _load_cases(path: Path | None = None) -> dict[str, Any]:
    p = Path(path or KNOWLEDGE_DIR / "comorbidity_clinical_benchmarks.json")
    if not p.exists():
        raise FileNotFoundError(
            f"Missing {p}. Run: python -m exhalepath harvest-clinical-comorbidity"
        )
    return json.loads(p.read_text())


def evaluate_comorbidity_clinical(
    *,
    top_k: int = 15,
    mode: str = "hybrid",
    comorbidity_weight: float = 0.65,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    """
    Score each clinical comorbidity case for elevated/suppressed VOC concordance.
    """
    doc = _load_cases()
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    cases_out = []

    for case in doc.get("cases") or []:
        report = engine.predict(
            case.get("disease") or case.get("disease_id"),
            location=case.get("location") or "brain",
            top_n=50,
            mode=mode,
            age_years=case.get("age_years"),
            sex=case.get("sex"),
            comorbidities=list(case.get("comorbidities") or []),
            comorbidity_weight=comorbidity_weight,
            explain=False,
        )
        ranked_ids = [p.voc_id for p in report.top_vocs]
        top_ids = set(ranked_ids[:top_k])
        by_id = {p.voc_id: p for p in report.result.bundle.predictions}

        elev = list(case.get("expect_elevated") or [])
        supp = list(case.get("expect_suppressed") or [])

        elev_in_top = [v for v in elev if v in top_ids]
        elev_dir_hits = elev_dir_tot = 0
        for v in elev:
            if v not in by_id:
                continue
            elev_dir_tot += 1
            if by_id[v].fold_change > 1.0:
                elev_dir_hits += 1
        supp_dir_hits = supp_dir_tot = 0
        for v in supp:
            if v not in by_id:
                continue
            supp_dir_tot += 1
            if by_id[v].fold_change < 1.0:
                supp_dir_hits += 1

        dir_tot = elev_dir_tot + supp_dir_tot
        dir_hits = elev_dir_hits + supp_dir_hits
        detail = []
        for v in elev:
            p = by_id.get(v)
            detail.append(
                {
                    "voc_id": v,
                    "role": "elevated",
                    "fold": None if not p else float(p.fold_change),
                    "delta_ppb": None if not p else float(p.delta_ppb),
                    "in_top_k": v in top_ids,
                    "direction_ok": bool(p and p.fold_change > 1.0),
                }
            )
        for v in supp:
            p = by_id.get(v)
            detail.append(
                {
                    "voc_id": v,
                    "role": "suppressed",
                    "fold": None if not p else float(p.fold_change),
                    "delta_ppb": None if not p else float(p.delta_ppb),
                    "in_top_k": v in top_ids,
                    "direction_ok": bool(p and p.fold_change < 1.0),
                }
            )

        cases_out.append(
            {
                "case_id": case["case_id"],
                "disease": case.get("disease"),
                "comorbidities": case.get("comorbidities") or [],
                "source": case.get("source"),
                "n_expect_elevated": len(elev),
                "n_expect_suppressed": len(supp),
                "elevated_recall_at_k": (len(elev_in_top) / len(elev)) if elev else None,
                "elevated_in_top_k": elev_in_top,
                "elevated_directional_accuracy": (elev_dir_hits / elev_dir_tot)
                if elev_dir_tot
                else None,
                "suppressed_directional_accuracy": (supp_dir_hits / supp_dir_tot)
                if supp_dir_tot
                else None,
                "directional_accuracy": (dir_hits / dir_tot) if dir_tot else None,
                "fused_comorbidities": (report.result.bundle.metadata or {}).get(
                    "comorbidities"
                ),
                "detail": detail,
            }
        )

    def _mean(vals: list[float | None]) -> float | None:
        xs = [float(v) for v in vals if v is not None]
        return sum(xs) / len(xs) if xs else None

    summary = {
        "version": "1.0.0",
        "mode": mode,
        "top_k": top_k,
        "comorbidity_weight": comorbidity_weight,
        "n_cases": len(cases_out),
        "mean_elevated_recall_at_k": _mean(
            [c["elevated_recall_at_k"] for c in cases_out]
        ),
        "mean_elevated_directional_accuracy": _mean(
            [c["elevated_directional_accuracy"] for c in cases_out]
        ),
        "mean_suppressed_directional_accuracy": _mean(
            [c["suppressed_directional_accuracy"] for c in cases_out]
        ),
        "mean_directional_accuracy": _mean(
            [c["directional_accuracy"] for c in cases_out]
        ),
        "datasets": doc.get("datasets"),
        "cases": cases_out,
    }

    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "comorbidity_clinical_eval.json").write_text(
            json.dumps(summary, indent=2)
        )
    return summary
