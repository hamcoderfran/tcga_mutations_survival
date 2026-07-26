"""Zero-shot reliability: rare held-out diseases + leave-disease-out prior audit."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR, PROCESSED_DIR
from ..knowledge.loader import KnowledgeBase
from ..model.predict import ExhalePathPredictor
from ..schemas import DiseaseQuery, TumorContext


def _load_profiles(path: Path | None = None) -> list[dict[str, Any]]:
    candidates = [
        path,
        KNOWLEDGE_DIR / "zero_shot_holdout_profiles.json",
        Path("data/knowledge/zero_shot_holdout_profiles.json"),
    ]
    for p in candidates:
        if p and Path(p).exists():
            doc = json.loads(Path(p).read_text())
            return list(doc.get("profiles") or [])
    raise FileNotFoundError("zero_shot_holdout_profiles.json not found")


def _pred_map(report) -> dict[str, float]:
    return {p.voc_id: float(p.log2_fold_change) for p in report.result.bundle.predictions}


def _pathway_map(report) -> dict[str, float]:
    return {p.pathway_id: float(p.score) for p in report.result.bundle.pathway_scores}


def evaluate_holdout_profiles(
    *,
    profiles_path: Path | None = None,
    mode: str = "hybrid",
) -> dict[str, Any]:
    engine = ExhaleBiomarkerEngine(use_opentargets=False)
    profiles = _load_profiles(profiles_path)
    rows: list[dict[str, Any]] = []
    n_pass = 0

    for prof in profiles:
        genes = [str(g).upper() for g in (prof.get("genes") or [])]
        report = engine.predict(
            prof["disease"],
            location=prof.get("location"),
            genes=genes or None,
            description=prof.get("description"),
            pathway_overrides=prof.get("pathway_overrides"),
            mode=mode,
            explain=False,
        )
        meta = report.result.bundle.metadata or {}
        pred = _pred_map(report)
        pw = _pathway_map(report)
        peak = max((abs(v) for v in pred.values()), default=0.0)
        checks: list[dict[str, Any]] = []

        def add(name: str, ok: bool, detail: Any = None) -> None:
            checks.append({"name": name, "pass": bool(ok), "detail": detail})

        unresolved = bool(meta.get("unresolved") or str(report.disease_id).startswith("custom::"))
        if prof.get("expect_unresolved") is not None:
            add("unresolved", unresolved == bool(prof["expect_unresolved"]), report.disease_id)

        mech_c = float(meta.get("mechanism_confidence") or 0.0)
        zs_mode = meta.get("zero_shot_mode")
        if prof.get("expect_mechanism") is True:
            add(
                "mechanism_mode",
                zs_mode == "mechanism" or mech_c >= float(prof.get("min_mechanism_confidence") or 0.5),
                {"mode": zs_mode, "mechanism_confidence": mech_c},
            )
        if prof.get("expect_mechanism") is False:
            add(
                "near_healthy_fallback",
                zs_mode == "near_healthy_fallback" or peak < float(prof.get("max_abs_log2fc") or 0.45),
                {"mode": zs_mode, "peak": peak},
            )

        if prof.get("min_mechanism_confidence") is not None and prof.get("expect_mechanism"):
            add(
                "min_mechanism_confidence",
                mech_c >= float(prof["min_mechanism_confidence"]),
                mech_c,
            )

        if prof.get("near_healthy"):
            thr = float(prof.get("max_abs_log2fc") or 0.45)
            add("near_healthy", peak <= thr, {"peak": peak, "threshold": thr})

        if prof.get("min_peak_abs_log2fc") is not None:
            add(
                "min_peak_abs_log2fc",
                peak >= float(prof["min_peak_abs_log2fc"]),
                peak,
            )

        for voc in prof.get("must_elevate") or []:
            add(f"elevate:{voc}", float(pred.get(voc, 0.0)) > 0.05, pred.get(voc))

        for voc in prof.get("must_suppress") or []:
            add(f"suppress:{voc}", float(pred.get(voc, 0.0)) < -0.05, pred.get(voc))

        if prof.get("pathway_any"):
            hit = any(float(pw.get(p, 0.0)) > 0.05 for p in prof["pathway_any"])
            add("pathway_any", hit, {p: pw.get(p) for p in prof["pathway_any"]})

        # Empty VOC prior for true zero-shot unresolved
        if unresolved:
            d = engine.kb.resolve_disease(prof["disease"], location_hint=prof.get("location"))
            add("no_voc_prior_leak", not bool(d.get("voc_log2fc_prior")), d.get("voc_log2fc_prior"))

        passed = all(c["pass"] for c in checks) if checks else False
        if passed:
            n_pass += 1
        top = sorted(pred.items(), key=lambda x: abs(x[1]), reverse=True)[:6]
        rows.append(
            {
                "id": prof.get("id"),
                "disease": prof["disease"],
                "disease_id": report.disease_id,
                "location": report.resolved_location,
                "passed": passed,
                "peak_abs_log2fc": peak,
                "mechanism_confidence": mech_c,
                "zero_shot_mode": zs_mode,
                "top_vocs": [{"voc_id": v, "log2fc": fc} for v, fc in top],
                "checks": checks,
                "notes": report.notes[:6],
            }
        )

    n = len(rows)
    return {
        "n_profiles": n,
        "n_passed": n_pass,
        "pass_rate": n_pass / n if n else 0.0,
        "pass": (n_pass / n if n else 0.0) >= 0.75,
        "profiles": rows,
    }


def evaluate_leave_disease_out_priors(
    *,
    voc_targets_path: Path | None = None,
    max_diseases: int | None = None,
) -> dict[str, Any]:
    """
    Leave-disease-out style audit: for each atlas disease with VOC priors,
    predict with voc_log2fc_prior masked to {} and score directional agreement
    vs the curated prior signs (mechanism/pathway generalization).
    """
    kb = KnowledgeBase()
    predictor = ExhalePathPredictor(knowledge=kb, use_opentargets=False)
    rows = []
    diseases = list(kb.diseases.values())
    if max_diseases is not None:
        diseases = diseases[: max(1, int(max_diseases))]

    for d in diseases:
        prior = dict(d.get("voc_log2fc_prior") or {})
        if len(prior) < 3:
            continue
        # Mask VOC prior on a shallow copy used only via temporary monkeypatch of resolve
        masked = dict(d)
        masked["voc_log2fc_prior"] = {}
        masked["_ldo_masked"] = True

        original = kb.diseases[d["disease_id"]]
        kb.diseases[d["disease_id"]] = masked
        try:
            res = predictor.predict(
                DiseaseQuery(
                    disease=d["name"],
                    tumor=TumorContext(primary_site=d.get("default_site")),
                    mode="hybrid",
                )
            )
        finally:
            kb.diseases[d["disease_id"]] = original

        pred = {p.voc_id: float(p.log2_fold_change) for p in res.bundle.predictions}
        signs = []
        for voc, truth in prior.items():
            if abs(float(truth)) < 0.15:
                continue
            if voc not in pred:
                continue
            signs.append(float(np.sign(pred[voc])) == float(np.sign(truth)))
        if len(signs) < 3:
            continue
        acc = float(np.mean(signs))
        rows.append(
            {
                "disease_id": d["disease_id"],
                "category": d.get("category"),
                "n_signed_priors": len(signs),
                "directional_accuracy": acc,
                "pass": acc >= 0.55,
            }
        )

    if not rows:
        return {"skipped": True, "reason": "no diseases with usable VOC priors"}

    mean_acc = float(np.mean([r["directional_accuracy"] for r in rows]))
    n_pass = sum(1 for r in rows if r["pass"])
    return {
        "skipped": False,
        "n_diseases": len(rows),
        "n_passed": n_pass,
        "pass_rate": n_pass / len(rows),
        "mean_directional_accuracy": mean_acc,
        "pass": mean_acc >= 0.55 and (n_pass / len(rows)) >= 0.6,
        "diseases": rows,
        "voc_targets_path": str(voc_targets_path or PROCESSED_DIR / "voc_training_targets.csv"),
    }


def evaluate_zero_shot_reliability(
    *,
    out_dir: Path | None = None,
    profiles_path: Path | None = None,
    mode: str = "hybrid",
    max_ldo_diseases: int | None = None,
) -> dict[str, Any]:
    holdout = evaluate_holdout_profiles(profiles_path=profiles_path, mode=mode)
    ldo = evaluate_leave_disease_out_priors(max_diseases=max_ldo_diseases)
    report = {
        "holdout": holdout,
        "leave_disease_out": ldo,
        "pass": bool(holdout.get("pass")) and bool(ldo.get("pass") or ldo.get("skipped")),
    }

    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "zero_shot_reliability.json").write_text(json.dumps(report, indent=2))
        lines = [
            "# Zero-shot reliability",
            "",
            f"**Passed:** {report['pass']}",
            "",
            "## Held-out rare diseases (no VOC labels)",
            f"- profiles: {holdout['n_passed']}/{holdout['n_profiles']} "
            f"({100 * holdout['pass_rate']:.0f}%)",
            "",
        ]
        for r in holdout["profiles"]:
            mark = "✓" if r["passed"] else "✗"
            lines.append(
                f"- {mark} **{r['id']}** `{r['disease']}` peak|log2fc|={r['peak_abs_log2fc']:.3f} "
                f"mech={r['mechanism_confidence']:.2f} mode={r['zero_shot_mode']}"
            )
        lines.extend(["", "## Leave-disease-out prior direction"])
        if ldo.get("skipped"):
            lines.append(f"- skipped: {ldo.get('reason')}")
        else:
            lines.append(
                f"- diseases: {ldo['n_passed']}/{ldo['n_diseases']} "
                f"(mean dir={ldo['mean_directional_accuracy']:.3f})"
            )
        (out / "ZERO_SHOT_RELIABILITY.md").write_text("\n".join(lines) + "\n")
    return report
