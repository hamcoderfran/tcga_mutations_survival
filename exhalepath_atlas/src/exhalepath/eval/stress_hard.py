"""Adversarial 50-case stress suite — try to break ExhalePath."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..body.tissues import default_body_map, resolve_location
from ..config import KNOWLEDGE_DIR
from ..knowledge.loader import clear_knowledge_cache


def _load(path: Path | None = None) -> dict[str, Any]:
    for p in [
        path,
        KNOWLEDGE_DIR / "stress_hard_50.json",
        Path("data/knowledge/stress_hard_50.json"),
    ]:
        if p and Path(p).exists():
            return json.loads(Path(p).read_text())
    raise FileNotFoundError("Missing stress_hard_50.json — run scripts/build_stress_hard_50.py")


def _finite_report(report) -> bool:
    for p in list(report.top_vocs or []):
        for attr in ("predicted_ppb", "healthy_ppb", "fold_change", "log2_fold_change", "delta_ppb"):
            v = getattr(p, attr, None)
            if v is None:
                continue
            if not math.isfinite(float(v)):
                return False
            if attr in ("predicted_ppb", "healthy_ppb") and float(v) <= 0:
                return False
    return True


def _by_voc(report) -> dict[str, Any]:
    out: dict[str, Any] = {}
    preds = list(getattr(getattr(report.result, "bundle", None), "predictions", None) or [])
    for p in preds:
        out[str(p.voc_id)] = p
    for p in report.top_vocs:
        out.setdefault(p.voc_id, p)
    return out


def _max_abs_log2(report) -> float:
    vals = [abs(float(p.log2_fold_change)) for p in report.top_vocs]
    return float(max(vals)) if vals else 0.0


def _check_case(case: dict[str, Any], engine: ExhaleBiomarkerEngine) -> dict[str, Any]:
    exp = dict(case.get("expect") or {})
    q = dict(case.get("query") or {})
    checks: list[dict[str, Any]] = []
    ok_all = True
    report = None
    err = None

    kw = {
        "disease": q.get("disease") if q.get("disease") is not None else "",
        "location": q.get("location"),
        "top_n": 50,
        "mode": "hybrid",
        "explain": False,
        "genes": q.get("genes"),
        "comorbidities": q.get("comorbidities"),
        "age_years": q.get("age_years"),
        "sex": q.get("sex"),
        "smoking_status": q.get("smoking_status"),
        "stage": q.get("stage"),
        "metastatic": bool(q.get("metastatic") or False),
    }
    if "comorbidity_weight" in q:
        kw["comorbidity_weight"] = q["comorbidity_weight"]

    try:
        report = engine.predict(**kw)
    except Exception as e:  # noqa: BLE001 — stress must catch everything
        err = f"{type(e).__name__}: {e}"

    def add(name: str, passed: bool, detail: Any = None) -> None:
        nonlocal ok_all
        checks.append({"check": name, "pass": bool(passed), "detail": detail})
        if not passed:
            ok_all = False

    if exp.get("no_exception", True):
        add("no_exception", err is None, err)
    if err is not None:
        return {
            "id": case["id"],
            "category": case.get("category"),
            "difficulty": case.get("difficulty"),
            "passed": False,
            "error": err,
            "checks": checks,
        }

    assert report is not None
    add("finite", _finite_report(report))

    did = report.disease_id
    unresolved = did.startswith("custom::") or bool(getattr(report, "notes", None) and "not in the curated" in str(report.notes))
    # also check bundle notes
    notes = " ".join(str(n) for n in (getattr(getattr(report.result, "bundle", None), "notes", None) or []))
    if "not in the curated" in notes or did.startswith("custom::"):
        unresolved = True

    if exp.get("unresolved") is True:
        add("unresolved", unresolved, did)
    if exp.get("resolved_disease_id"):
        add("resolved_disease_id", did == exp["resolved_disease_id"], {"got": did, "want": exp["resolved_disease_id"]})
    if exp.get("resolved_in"):
        add("resolved_in", did in set(exp["resolved_in"]), {"got": did, "want": exp["resolved_in"]})
    if exp.get("unresolved_or_not_breast"):
        add(
            "unresolved_or_not_breast",
            unresolved or did != "breast_invasive_carcinoma",
            did,
        )

    # location expectations (BiomarkerReport.location is the tissue dict)
    loc = getattr(report, "location", None)
    if isinstance(loc, dict) and loc:
        tid = str(loc.get("tissue_id") or "")
        matched = bool(loc.get("matched"))
    else:
        default_body_map.cache_clear()
        loc_r = resolve_location(q.get("location"))
        tid = str(loc_r.get("tissue_id") or "")
        matched = bool(loc_r.get("matched"))
        loc = loc_r

    if exp.get("location_unmatched") is True:
        add("location_unmatched", matched is False, loc)
    if exp.get("location_matched") is True:
        add("location_matched", matched is True, loc)
    if exp.get("location_tissue_id"):
        add("location_tissue_id", tid == exp["location_tissue_id"], {"got": tid, "want": exp["location_tissue_id"]})
    if exp.get("location_contains"):
        add("location_contains", exp["location_contains"] in tid, {"got": tid})

    by = _by_voc(report)
    if exp.get("near_healthy"):
        thr = float(exp.get("max_abs_log2fc") or 0.4)
        m = _max_abs_log2(report)
        add("near_healthy", m <= thr, {"max_abs_log2fc": m, "threshold": thr})
    elif exp.get("max_abs_log2fc") is not None:
        thr = float(exp["max_abs_log2fc"])
        m = _max_abs_log2(report)
        add("max_abs_log2fc", m <= thr, {"max_abs_log2fc": m, "threshold": thr})

    for voc in exp.get("must_elevate") or []:
        p = by.get(voc)
        fold = float(p.fold_change) if p is not None else None
        add(f"elevate:{voc}", fold is not None and fold > 1.0, fold)
    for voc in exp.get("must_suppress") or []:
        p = by.get(voc)
        fold = float(p.fold_change) if p is not None else None
        add(f"suppress:{voc}", fold is not None and fold < 1.0, fold)
    for voc, thr in (exp.get("min_fold") or {}).items():
        p = by.get(voc)
        fold = float(p.fold_change) if p is not None else None
        add(f"min_fold:{voc}", fold is not None and fold >= float(thr), {"fold": fold, "thr": thr})

    if exp.get("finite"):
        add("finite_explicit", _finite_report(report))

    return {
        "id": case["id"],
        "category": case.get("category"),
        "difficulty": case.get("difficulty"),
        "passed": ok_all,
        "resolved_disease_id": did,
        "location": loc if isinstance(loc, dict) else {"tissue_id": tid, "matched": matched},
        "max_abs_log2fc": _max_abs_log2(report),
        "top_vocs": [p.voc_id for p in report.top_vocs[:8]],
        "checks": checks,
        "notes": case.get("notes"),
    }


def evaluate_stress_hard(
    *,
    out_dir: Path | None = None,
    profiles_path: Path | None = None,
) -> dict[str, Any]:
    clear_knowledge_cache()
    default_body_map.cache_clear()
    doc = _load(profiles_path)
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)

    results = [_check_case(c, engine) for c in doc["cases"]]
    n = len(results)
    n_pass = sum(1 for r in results if r["passed"])
    by_cat: dict[str, dict[str, Any]] = {}
    for r in results:
        cat = r.get("category") or "other"
        by_cat.setdefault(cat, {"n": 0, "passed": 0})
        by_cat[cat]["n"] += 1
        by_cat[cat]["passed"] += int(r["passed"])
    for cat, row in by_cat.items():
        row["pass_rate"] = float(row["passed"] / row["n"]) if row["n"] else None

    failed = [r for r in results if not r["passed"]]
    # failure mode tallies
    fail_checks: dict[str, int] = {}
    for r in failed:
        for ch in r.get("checks") or []:
            if not ch.get("pass"):
                key = str(ch.get("check") or "?")
                # normalize elevate/suppress
                if key.startswith("elevate:"):
                    key = "wrong_or_missing_elevate"
                elif key.startswith("suppress:"):
                    key = "wrong_or_missing_suppress"
                elif key.startswith("min_fold:"):
                    key = "min_fold_miss"
                fail_checks[key] = fail_checks.get(key, 0) + 1

    overall = {
        "n_cases": n,
        "n_passed": n_pass,
        "n_failed": n - n_pass,
        "pass_rate": float(n_pass / n) if n else None,
        "pass_pct": None if not n else round(100.0 * n_pass / n, 2),
        "by_category": by_cat,
        "failure_modes": fail_checks,
        "failed_ids": [r["id"] for r in failed],
    }
    report = {"overall": overall, "cases": results, "spec_version": doc.get("version")}
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "stress_hard_eval.json").write_text(json.dumps(report, indent=2) + "\n")
        (out_dir / "STRESS_HARD.md").write_text(_md(report))
    return report


def _md(report: dict[str, Any]) -> str:
    o = report["overall"]
    lines = [
        "# ExhalePath adversarial stress suite (50 hard cases)",
        "",
        f"**Pass rate: {o['pass_pct']}%** ({o['n_passed']}/{o['n_cases']})",
        "",
        "## By category",
        "",
    ]
    for cat, row in (o.get("by_category") or {}).items():
        lines.append(f"- **{cat}**: {row['passed']}/{row['n']} ({100*row['pass_rate']:.0f}%)")
    lines += ["", "## Failure modes", ""]
    for k, v in sorted((o.get("failure_modes") or {}).items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{k}`: {v}")
    lines += ["", "## Failed cases", ""]
    for r in report["cases"]:
        if r["passed"]:
            continue
        bad = [c["check"] for c in r.get("checks") or [] if not c.get("pass")]
        lines.append(f"- **{r['id']}** ({r.get('category')}): {', '.join(bad)}")
        if r.get("error"):
            lines.append(f"  - error: `{r['error']}`")
    lines += ["", "## All cases", ""]
    lines.append("| ID | Cat | Pass | Disease | max\\|log2\\| |")
    lines.append("|---|---|---|---|---:|")
    for r in report["cases"]:
        lines.append(
            f"| {r['id']} | {r.get('category')} | {'✓' if r['passed'] else '✗'} | "
            f"{r.get('resolved_disease_id','')} | {r.get('max_abs_log2fc', 0):.3f} |"
        )
    lines.append("")
    return "\n".join(lines)
