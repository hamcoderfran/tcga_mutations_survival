"""50-disease vision evaluation suite for ExhalePath Atlas.

Scores each patient profile on:
  - directional VOC accuracy (must/should elevate & suppress)
  - min-fold gates
  - top-k recall of hallmark VOCs
  - pathway explainability
  - cell-state alignment
  - anatomic site alignment

Headline metrics
----------------
vision_fidelity_pct     — all 50 profiles (includes Grade D prior-aligned)
evidence_backed_pct     — Grade A+B only
held_out_style_pct      — Grade A+B+C (excludes prior-only Grade D)
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR

_TOKEN = re.compile(r"[a-z0-9]+")


def _norm(s: str) -> str:
    return " ".join(_TOKEN.findall(str(s).lower()))


def _fuzzy_hit(expected: str, candidates: list[str]) -> bool:
    e = _norm(expected)
    if not e:
        return False
    for c in candidates:
        n = _norm(c)
        if not n:
            continue
        if e in n or n in e:
            return True
        # token overlap ≥ 50% of expected tokens
        et = set(e.split())
        ct = set(n.split())
        if et and len(et & ct) / len(et) >= 0.5:
            return True
    return False


def _load_profiles(path: Path | None = None) -> dict[str, Any]:
    candidates = [
        path,
        KNOWLEDGE_DIR / "vision_50_profiles.json",
        Path("data/knowledge/vision_50_profiles.json"),
    ]
    for p in candidates:
        if p and Path(p).exists():
            return json.loads(Path(p).read_text())
    raise FileNotFoundError(
        "Missing vision_50_profiles.json. Run: python scripts/build_vision_50_profiles.py"
    )


def _direction_score(fold: float | None, *, want_up: bool, soft: bool) -> float:
    if fold is None:
        return 0.0
    if want_up:
        if fold >= 1.15:
            return 1.0
        if fold >= 1.05:
            return 0.55 if not soft else 0.7
        if fold > 1.0:
            return 0.25 if not soft else 0.4
        return 0.0
    # suppress
    if fold <= 0.87:
        return 1.0
    if fold <= 0.95:
        return 0.55 if not soft else 0.7
    if fold < 1.0:
        return 0.25 if not soft else 0.4
    return 0.0


def _score_direction(gt: dict[str, Any], by_id: dict[str, Any]) -> dict[str, Any]:
    parts: list[tuple[float, float]] = []  # (score, weight)
    detail: list[dict[str, Any]] = []

    def _check(vocs: list[str], *, want_up: bool, soft: bool, weight: float) -> None:
        for v in vocs:
            pred = by_id.get(v)
            fold = float(pred.fold_change) if pred is not None else None
            s = _direction_score(fold, want_up=want_up, soft=soft)
            parts.append((s, weight))
            detail.append(
                {
                    "voc_id": v,
                    "want": "elevate" if want_up else "suppress",
                    "soft": soft,
                    "fold_change": fold,
                    "score": s,
                    "weight": weight,
                }
            )

    _check(list(gt.get("must_elevate") or []), want_up=True, soft=False, weight=1.0)
    _check(list(gt.get("should_elevate") or []), want_up=True, soft=True, weight=0.5)
    _check(list(gt.get("must_suppress") or []), want_up=False, soft=False, weight=1.0)
    _check(list(gt.get("should_suppress") or []), want_up=False, soft=True, weight=0.5)

    if not parts:
        return {"score": None, "detail": detail}
    num = sum(s * w for s, w in parts)
    den = sum(w for _, w in parts)
    return {"score": float(num / den), "detail": detail}


def _score_fold(gt: dict[str, Any], by_id: dict[str, Any]) -> dict[str, Any]:
    min_fold = dict(gt.get("min_fold") or {})
    if not min_fold:
        return {"score": None, "detail": []}
    detail = []
    scores = []
    for voc, thr in min_fold.items():
        thr_f = float(thr)
        pred = by_id.get(voc)
        fold = float(pred.fold_change) if pred is not None else None
        if fold is None:
            s = 0.0
        elif fold >= thr_f:
            s = 1.0
        else:
            # partial credit toward threshold (from 1.0 baseline)
            s = max(0.0, min(1.0, (fold - 1.0) / max(1e-6, thr_f - 1.0)))
        scores.append(s)
        detail.append({"voc_id": voc, "threshold": thr_f, "fold_change": fold, "score": s})
    return {"score": float(sum(scores) / len(scores)), "detail": detail}


def _score_topk(gt: dict[str, Any], ranked_ids: list[str], top_k: int) -> dict[str, Any]:
    need = list(gt.get("top_k_should_include") or [])
    if not need:
        return {"score": None, "in_top": [], "missing": []}
    top = set(ranked_ids[:top_k])
    hits = [v for v in need if v in top]
    return {
        "score": float(len(hits) / len(need)),
        "in_top": hits,
        "missing": [v for v in need if v not in top],
    }


def _score_pathways(gt: dict[str, Any], report) -> dict[str, Any]:
    need = list(gt.get("pathways_should_include") or [])
    if not need:
        return {"score": None, "hits": [], "missing": []}

    present: set[str] = set()
    # Bundle pathway scores (disease-biased)
    bundle = getattr(getattr(report, "result", None), "bundle", None)
    for ps in list(getattr(bundle, "pathway_scores", None) or []):
        pid = getattr(ps, "pathway_id", None) or (ps.get("pathway_id") if isinstance(ps, dict) else None)
        score = getattr(ps, "score", None)
        bias = getattr(ps, "disease_bias", None)
        if score is None and isinstance(ps, dict):
            score = ps.get("score")
            bias = ps.get("disease_bias")
        if pid and ((score or 0) > 0 or (bias or 1.0) > 1.05):
            present.add(str(pid))

    # Mechanism pathways with contribution
    for m in list(getattr(report, "mechanisms", None) or []):
        if not isinstance(m, dict):
            continue
        for pw in m.get("pathways") or []:
            if not isinstance(pw, dict):
                continue
            pid = pw.get("pathway_id")
            if pid and float(pw.get("contribution") or pw.get("score") or 0) > 0:
                present.add(str(pid))

    hits = [p for p in need if p in present]
    return {
        "score": float(len(hits) / len(need)),
        "hits": hits,
        "missing": [p for p in need if p not in present],
        "present": sorted(present),
    }


def _score_cells(gt: dict[str, Any], report) -> dict[str, Any]:
    need = list(gt.get("cell_states_should_include") or [])
    if not need:
        return {"score": None, "hits": [], "missing": []}

    names: list[str] = []
    bundle = getattr(getattr(report, "result", None), "bundle", None)
    for cs in list(getattr(bundle, "cell_states", None) or [])[:12]:
        names.append(getattr(cs, "name", None) or (cs.get("name") if isinstance(cs, dict) else "") or "")
        names.append(getattr(cs, "state_id", None) or (cs.get("state_id") if isinstance(cs, dict) else "") or "")
    for m in list(getattr(report, "mechanisms", None) or [])[:8]:
        if not isinstance(m, dict):
            continue
        for cs in m.get("cell_states") or []:
            if isinstance(cs, dict):
                names.append(str(cs.get("name") or ""))
                names.append(str(cs.get("state_id") or ""))

    hits = [e for e in need if _fuzzy_hit(e, names)]
    return {
        "score": float(len(hits) / len(need)),
        "hits": hits,
        "missing": [e for e in need if e not in hits],
        "observed_top": [n for n in names if n][:10],
    }


def _score_sites(gt: dict[str, Any], report) -> dict[str, Any]:
    need = list(gt.get("sites_should_include") or [])
    if not need:
        return {"score": None, "hits": [], "missing": []}

    observed: list[str] = []
    loc = getattr(report, "resolved_location", None) or getattr(report, "location", None)
    if isinstance(loc, dict):
        observed.extend([str(loc.get(k) or "") for k in ("site", "tissue", "name", "id")])
    elif loc:
        observed.append(str(loc))
    tissue = getattr(report, "tissue", None)
    if tissue:
        observed.append(str(tissue))

    bundle = getattr(getattr(report, "result", None), "bundle", None)
    for cs in list(getattr(bundle, "cell_states", None) or [])[:12]:
        t = getattr(cs, "tissue", None) or (cs.get("tissue") if isinstance(cs, dict) else None)
        if t:
            observed.append(str(t))
    for m in list(getattr(report, "mechanisms", None) or [])[:8]:
        if not isinstance(m, dict):
            continue
        for cs in m.get("cell_states") or []:
            if isinstance(cs, dict) and cs.get("tissue"):
                observed.append(str(cs["tissue"]))

    hits = [e for e in need if _fuzzy_hit(e, observed)]
    return {
        "score": float(len(hits) / len(need)),
        "hits": hits,
        "missing": [e for e in need if e not in hits],
        "observed": sorted({_norm(x) for x in observed if x}),
    }


def _composite(parts: dict[str, float | None], weights: dict[str, float]) -> float | None:
    num = den = 0.0
    for k, w in weights.items():
        v = parts.get(k)
        if v is None:
            continue
        num += float(v) * w
        den += w
    return float(num / den) if den else None


DEFAULT_WEIGHTS = {
    "direction": 1.0,
    "fold": 1.0,
    "topk": 0.9,
    "pathway": 0.75,
    "cell": 0.5,
    "site": 0.45,
}


def evaluate_vision(
    *,
    top_k: int = 15,
    mode: str = "hybrid",
    out_dir: Path | None = None,
    profiles_path: Path | None = None,
    limit: int | None = None,
) -> dict[str, Any]:
    """Run the full 50-disease vision suite and optionally write JSON + Markdown."""
    doc = _load_profiles(profiles_path)
    profiles = list(doc["profiles"])
    if limit is not None:
        profiles = profiles[: int(limit)]

    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    cases: list[dict[str, Any]] = []

    for prof in profiles:
        patient = dict(prof.get("patient") or {})
        gt = dict(prof.get("ground_truth") or {})
        report = engine.predict(
            prof.get("query_name") or prof.get("name") or prof["id"],
            location=prof.get("location"),
            top_n=50,
            mode=mode,
            sex=patient.get("sex"),
            age_years=patient.get("age"),
            smoking_status=patient.get("smoking_status"),
            genes=list(patient.get("genes") or []),
            comorbidities=list(patient.get("comorbidities") or []),
            explain=True,
        )
        ranked = [p.voc_id for p in report.top_vocs]
        by_id = {p.voc_id: p for p in report.top_vocs}

        direction = _score_direction(gt, by_id)
        fold = _score_fold(gt, by_id)
        topk = _score_topk(gt, ranked, top_k)
        pathway = _score_pathways(gt, report)
        cell = _score_cells(gt, report)
        site = _score_sites(gt, report)

        parts = {
            "direction": direction["score"],
            "fold": fold["score"],
            "topk": topk["score"],
            "pathway": pathway["score"],
            "cell": cell["score"],
            "site": site["score"],
        }
        comp = _composite(parts, DEFAULT_WEIGHTS)
        cases.append(
            {
                "id": prof["id"],
                "name": prof.get("name"),
                "grade": prof.get("grade"),
                "location": prof.get("location"),
                "resolved_disease_id": report.disease_id,
                "evidence_basis": prof.get("evidence_basis"),
                "patient": patient,
                "scores": {**parts, "composite": comp},
                "composite_pct": None if comp is None else round(100.0 * comp, 2),
                "direction_detail": direction["detail"],
                "fold_detail": fold["detail"],
                "topk": topk,
                "pathway": pathway,
                "cell": cell,
                "site": site,
                "predicted_top": ranked[:top_k],
            }
        )

    def _mean_comp(rows: list[dict[str, Any]]) -> float | None:
        vals = [c["scores"]["composite"] for c in rows if c["scores"].get("composite") is not None]
        return float(sum(vals) / len(vals)) if vals else None

    def _mean_metric(rows: list[dict[str, Any]], key: str) -> float | None:
        vals = [c["scores"].get(key) for c in rows if c["scores"].get(key) is not None]
        return float(sum(vals) / len(vals)) if vals else None

    all_rows = cases
    ab = [c for c in cases if c.get("grade") in ("A", "B")]
    abc = [c for c in cases if c.get("grade") in ("A", "B", "C")]
    d_only = [c for c in cases if c.get("grade") == "D"]
    by_grade: dict[str, Any] = {}
    for g in ("A", "B", "C", "D"):
        rows = [c for c in cases if c.get("grade") == g]
        by_grade[g] = {
            "n": len(rows),
            "mean_composite": _mean_comp(rows),
            "mean_composite_pct": None if _mean_comp(rows) is None else round(100.0 * _mean_comp(rows), 2),
            "mean_direction": _mean_metric(rows, "direction"),
            "mean_fold": _mean_metric(rows, "fold"),
            "mean_topk": _mean_metric(rows, "topk"),
            "mean_pathway": _mean_metric(rows, "pathway"),
        }

    vision = _mean_comp(all_rows)
    evidence = _mean_comp(ab)
    held = _mean_comp(abc)
    prior_only = _mean_comp(d_only)

    overall = {
        "n_profiles": len(cases),
        "top_k": top_k,
        "mode": mode,
        "vision_fidelity": vision,
        "vision_fidelity_pct": None if vision is None else round(100.0 * vision, 2),
        "evidence_backed": evidence,
        "evidence_backed_pct": None if evidence is None else round(100.0 * evidence, 2),
        "held_out_style": held,
        "held_out_style_pct": None if held is None else round(100.0 * held, 2),
        "grade_D_prior_consistency": prior_only,
        "grade_D_prior_consistency_pct": None
        if prior_only is None
        else round(100.0 * prior_only, 2),
        "metric_means": {
            "direction": _mean_metric(all_rows, "direction"),
            "fold": _mean_metric(all_rows, "fold"),
            "topk": _mean_metric(all_rows, "topk"),
            "pathway": _mean_metric(all_rows, "pathway"),
            "cell": _mean_metric(all_rows, "cell"),
            "site": _mean_metric(all_rows, "site"),
        },
        "by_grade": by_grade,
        "weights": DEFAULT_WEIGHTS,
        "scoring_notes": doc.get("scoring_notes"),
        "grade_counts": doc.get("grade_counts"),
    }

    # Pass/fail helper vs vision targets
    overall["verdict"] = {
        "vision_fidelity_pct": overall["vision_fidelity_pct"],
        "primary_score_for_original_vision": overall["vision_fidelity_pct"],
        "interpretation": (
            "vision_fidelity_pct is the headline % closeness to the original vision "
            "(disease+location+patient → VOC ranking + pathway/cell/site explainability) "
            "across all 50 profiles. Prefer evidence_backed_pct for independent truthfulness; "
            "Grade D inflates fidelity via prior circularity."
        ),
    }

    result = {"overall": overall, "cases": cases}
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "vision_eval.json").write_text(json.dumps(result, indent=2) + "\n")
        (out_dir / "VISION_EVAL.md").write_text(_markdown_report(result))
    return result


def _markdown_report(result: dict[str, Any]) -> str:
    o = result["overall"]
    lines = [
        "# ExhalePath vision evaluation (50 diseases)",
        "",
        f"**Vision fidelity (all 50): {o['vision_fidelity_pct']}%**",
        f"**Evidence-backed (Grade A+B): {o['evidence_backed_pct']}%**",
        f"**Held-out style (Grade A+B+C): {o['held_out_style_pct']}%**",
        f"**Grade D prior consistency: {o['grade_D_prior_consistency_pct']}%**",
        "",
        "## Metric means (all profiles)",
        "",
    ]
    for k, v in (o.get("metric_means") or {}).items():
        if v is None:
            lines.append(f"- {k}: n/a")
        else:
            lines.append(f"- {k}: {100.0 * v:.1f}%")
    lines += ["", "## By evidence grade", ""]
    for g, row in (o.get("by_grade") or {}).items():
        lines.append(
            f"- **Grade {g}** (n={row['n']}): composite {row['mean_composite_pct']}% "
            f"(dir={_pct(row.get('mean_direction'))}, fold={_pct(row.get('mean_fold'))}, "
            f"topk={_pct(row.get('mean_topk'))}, pathway={_pct(row.get('mean_pathway'))})"
        )
    lines += ["", "## Per-disease composite", ""]
    lines.append("| ID | Grade | Composite % | Direction | Fold | Top-k | Pathway | Cell | Site |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for c in result["cases"]:
        s = c["scores"]
        lines.append(
            f"| {c['id']} | {c.get('grade')} | {c.get('composite_pct')} | "
            f"{_pct(s.get('direction'))} | {_pct(s.get('fold'))} | {_pct(s.get('topk'))} | "
            f"{_pct(s.get('pathway'))} | {_pct(s.get('cell'))} | {_pct(s.get('site'))} |"
        )
    lines += [
        "",
        "## Notes",
        "",
        "- Grade **A**: quantified measured fold changes.",
        "- Grade **B**: strong directional literature / priority panels.",
        "- Grade **C**: held-out / zero-shot literature expectations.",
        "- Grade **D**: atlas-prior–aligned approximate GT (consistency, not independent accuracy).",
        "",
        str((o.get("verdict") or {}).get("interpretation") or ""),
        "",
    ]
    return "\n".join(lines)


def _pct(v: float | None) -> str:
    if v is None:
        return "—"
    return f"{100.0 * float(v):.0f}%"
