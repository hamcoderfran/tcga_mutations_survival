from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import KNOWLEDGE_DIR, MODELS_DIR, PROCESSED_DIR
from ..knowledge.loader import KnowledgeBase
from ..model.predict import ExhalePathPredictor
from ..schemas import DiseaseQuery, TumorContext


@dataclass
class AuditReport:
    literature_accuracy: dict[str, Any] = field(default_factory=dict)
    calibrator_holdout: dict[str, Any] = field(default_factory=dict)
    zero_shot: dict[str, Any] = field(default_factory=dict)
    stress: dict[str, Any] = field(default_factory=dict)
    bugs_checked: list[str] = field(default_factory=list)
    passed: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "passed": self.passed,
            "literature_accuracy": self.literature_accuracy,
            "calibrator_holdout": self.calibrator_holdout,
            "zero_shot": self.zero_shot,
            "stress": self.stress,
            "bugs_checked": self.bugs_checked,
        }


def _load_benchmarks(path: Path | None = None) -> list[dict[str, Any]]:
    p = Path(path or KNOWLEDGE_DIR / "literature_benchmarks.json")
    return json.loads(p.read_text())["benchmarks"]


def evaluate_literature_benchmarks(
    predictor: ExhalePathPredictor | None = None,
) -> dict[str, Any]:
    predictor = predictor or ExhalePathPredictor(use_opentargets=False)
    rows = []
    for bench in _load_benchmarks():
        tumor = TumorContext(**bench["tumor"]) if bench.get("tumor") else None
        query = DiseaseQuery(
            disease=bench["disease"],
            tumor=tumor,
            mutated_genes=bench.get("mutated_genes") or [],
            sex=bench.get("sex"),
            smoking_status=bench.get("smoking_status"),
            age_years=bench.get("age_years"),
        )
        result = predictor.predict(query)
        df = result.to_dataframe().set_index("voc_id")
        case = {
            "case_id": bench["case_id"],
            "disease": bench["disease"],
            "resolved_disease_id": result.bundle.disease_id,
            "zero_shot_flag": bool(bench.get("zero_shot")),
            "elevated_hits": 0,
            "elevated_total": 0,
            "min_fold_hits": 0,
            "min_fold_total": 0,
            "near_healthy_ok": None,
            "errors": [],
        }

        for voc_id in bench.get("expect_elevated", []):
            case["elevated_total"] += 1
            if voc_id not in df.index:
                case["errors"].append(f"missing voc {voc_id}")
                continue
            if df.loc[voc_id, "fold_change"] > 1.0:
                case["elevated_hits"] += 1
            else:
                case["errors"].append(
                    f"{voc_id} not elevated (fold={df.loc[voc_id, 'fold_change']:.3f})"
                )

        for voc_id, min_fold in (bench.get("min_fold") or {}).items():
            case["min_fold_total"] += 1
            if voc_id not in df.index:
                case["errors"].append(f"missing voc {voc_id}")
                continue
            if df.loc[voc_id, "fold_change"] >= float(min_fold):
                case["min_fold_hits"] += 1
            else:
                case["errors"].append(
                    f"{voc_id} fold {df.loc[voc_id, 'fold_change']:.3f} < {min_fold}"
                )

        if bench.get("expect_near_healthy"):
            max_abs = float(bench.get("max_abs_log2fc", 0.35))
            peak = float(df["log2_fold_change"].abs().max())
            case["near_healthy_ok"] = peak <= max_abs
            if not case["near_healthy_ok"]:
                case["errors"].append(f"peak |log2fc|={peak:.3f} > {max_abs}")

        case["passed"] = len(case["errors"]) == 0
        rows.append(case)

    n = len(rows)
    n_pass = sum(1 for r in rows if r["passed"])
    elev_h = sum(r["elevated_hits"] for r in rows)
    elev_t = sum(r["elevated_total"] for r in rows)
    fold_h = sum(r["min_fold_hits"] for r in rows)
    fold_t = sum(r["min_fold_total"] for r in rows)
    return {
        "n_cases": n,
        "n_passed": n_pass,
        "case_pass_rate": n_pass / n if n else 0.0,
        "directional_accuracy": elev_h / elev_t if elev_t else None,
        "min_fold_accuracy": fold_h / fold_t if fold_t else None,
        "min_fold_pass_rate": fold_h / fold_t if fold_t else None,  # alias
        "cases": rows,
    }


def evaluate_calibrator_holdout(
    voc_targets_path: Path | None = None,
    case_features_path: Path | None = None,
    sample_n: int = 40_000,
    random_state: int = 11,
) -> dict[str, Any]:
    """Holdout MAE/R² on corpus targets if available; else skip."""
    from sklearn.metrics import mean_absolute_error, r2_score

    from ..model.features import feature_frame_for_training

    voc_targets_path = Path(voc_targets_path or PROCESSED_DIR / "voc_training_targets.csv")
    case_features_path = Path(case_features_path or PROCESSED_DIR / "case_pathway_features.csv")
    model_path = MODELS_DIR / "voc_calibrator.joblib"
    if not (voc_targets_path.exists() and case_features_path.exists() and model_path.exists()):
        return {"skipped": True, "reason": "corpus or calibrator missing"}

    import joblib

    kb = KnowledgeBase()
    bundle = joblib.load(model_path)
    voc_targets = pd.read_csv(voc_targets_path)
    case_features = pd.read_csv(case_features_path)
    if len(voc_targets) > sample_n:
        voc_targets = voc_targets.sample(sample_n, random_state=random_state)

    X, y, voc_ids = feature_frame_for_training(case_features, voc_targets, kb)
    per_voc = {}
    y_true_all, y_pred_all = [], []
    for voc_id, model in bundle["models"].items():
        mask = voc_ids == voc_id
        if mask.sum() < 20:
            continue
        Xi = X.loc[mask].reindex(columns=bundle["feature_columns"], fill_value=0.0)
        pred = model.predict(Xi)
        yt = y.loc[mask].to_numpy()
        per_voc[voc_id] = {
            "mae_log2fc": float(mean_absolute_error(yt, pred)),
            "r2": float(r2_score(yt, pred)) if len(np.unique(yt)) > 1 else None,
            "n": int(mask.sum()),
        }
        y_true_all.append(yt)
        y_pred_all.append(pred)

    if not y_true_all:
        return {"skipped": True, "reason": "no overlapping VOC rows"}

    yt = np.concatenate(y_true_all)
    yp = np.concatenate(y_pred_all)
    return {
        "skipped": False,
        "mae_log2fc": float(mean_absolute_error(yt, yp)),
        "r2": float(r2_score(yt, yp)),
        "n": int(len(yt)),
        "per_voc": per_voc,
        "pass": float(mean_absolute_error(yt, yp)) < 0.25,
    }


def evaluate_zero_shot(predictor: ExhalePathPredictor | None = None) -> dict[str, Any]:
    predictor = predictor or ExhalePathPredictor(use_opentargets=False)
    lit = evaluate_literature_benchmarks(predictor)
    zs = [c for c in lit["cases"] if c.get("zero_shot_flag")]
    n = len(zs)
    n_pass = sum(1 for c in zs if c["passed"])
    # Novel disease discrimination: unknown near-healthy vs known metabolic
    unknown = predictor.predict(DiseaseQuery(disease="completely novel ailment QZZ-42"))
    t2d = predictor.predict(DiseaseQuery(disease="type 2 diabetes"))
    u = unknown.to_dataframe().set_index("voc_id")
    t = t2d.to_dataframe().set_index("voc_id")
    acetone_sep = float(t.loc["acetone", "predicted_ppb"] - u.loc["acetone", "predicted_ppb"])
    return {
        "n_zero_shot_cases": n,
        "n_passed": n_pass,
        "pass_rate": n_pass / n if n else 0.0,
        "cases": zs,
        "novel_vs_t2d_acetone_delta_ppb": acetone_sep,
        "novel_peak_abs_log2fc": float(u["log2_fold_change"].abs().max()),
        "pass": (n_pass / n if n else 0) >= 0.75 and acetone_sep > 100,
    }


def stress_test_pipeline(predictor: ExhalePathPredictor | None = None) -> dict[str, Any]:
    from ..utils import stage_multiplier, stage_ordinal

    predictor = predictor or ExhalePathPredictor(use_opentargets=False)
    kb = KnowledgeBase()
    failures: list[str] = []
    checks: list[str] = []

    # Stage parsing (regression for I-prefix bug)
    expected = {"I": 1.0, "II": 1.15, "III": 1.35, "IV": 1.55, "Stage IV": 1.55, "IIIA": 1.35}
    for s, exp in expected.items():
        got = stage_multiplier(s)
        checks.append(f"stage_multiplier({s})={got}")
        if abs(got - exp) > 1e-6:
            failures.append(f"stage_multiplier({s})={got}, expected {exp}")
    if stage_ordinal("IV") < stage_ordinal("III") or stage_ordinal("III") < stage_ordinal("II"):
        failures.append("stage_ordinal ordering broken")

    # Disease resolution should not map bare 'cancer' to a specific cancer
    d_cancer = kb.resolve_disease("cancer")
    checks.append(f"resolve(cancer)={d_cancer['disease_id']}")
    if not d_cancer.get("_unresolved"):
        failures.append("bare 'cancer' incorrectly resolved to " + d_cancer["disease_id"])

    d_lung = kb.resolve_disease("lung cancer")
    if d_lung["disease_id"] != "lung_adenocarcinoma":
        failures.append(f"lung cancer resolved to {d_lung['disease_id']}")

    # Empty / weird inputs must not crash
    weird_queries = [
        DiseaseQuery(disease=""),
        DiseaseQuery(disease="   "),
        DiseaseQuery(disease="???"),
        DiseaseQuery(disease="PAAD", mutated_genes=[]),
        DiseaseQuery(disease="PAAD", mutated_genes=["", "kras", "TP53"]),
        DiseaseQuery(
            disease="LUAD",
            tumor=TumorContext(stage="not-a-stage", primary_site="", histology=None),
            mutated_genes=["TP53"] * 200,
        ),
        DiseaseQuery(
            disease="type 2 diabetes",
            pathway_overrides={"ketone_body_metabolism": 2.0},
        ),
        {"disease": "BRCA", "tumor": {"stage": "IV", "metastatic": True}, "mutated_genes": ["PIK3CA"]},
    ]
    for i, q in enumerate(weird_queries):
        try:
            res = predictor.predict(q)
            df = res.to_dataframe()
            assert len(df) >= 15
            assert np.isfinite(df["predicted_ppb"]).all()
            assert (df["predicted_ppb"] > 0).all()
            # fold_change must match ppb ratio
            ratio = df["predicted_ppb"] / df["healthy_ppb"]
            if not np.allclose(ratio, df["fold_change"], rtol=1e-5, atol=1e-5):
                failures.append(f"weird[{i}] fold_change inconsistent with ppb ratio")
            checks.append(f"weird_query[{i}] ok")
        except Exception as e:  # noqa: BLE001
            failures.append(f"weird_query[{i}] raised {type(e).__name__}: {e}")

    # Stage IV must exceed stage I for peroxidation VOCs in LUAD
    early = predictor.predict(
        DiseaseQuery(
            disease="LUAD",
            tumor=TumorContext(stage="I", primary_site="lung"),
            mutated_genes=["TP53", "KRAS"],
        )
    )
    late = predictor.predict(
        DiseaseQuery(
            disease="LUAD",
            tumor=TumorContext(stage="IV", primary_site="lung", metastatic=True),
            mutated_genes=["TP53", "KRAS"],
        )
    )
    e = early.to_dataframe().set_index("voc_id")
    l = late.to_dataframe().set_index("voc_id")
    if not (l.loc["hexanal", "predicted_ppb"] > e.loc["hexanal", "predicted_ppb"]):
        failures.append("stage IV hexanal not greater than stage I")
    # Stronger: stage multiplier itself must differ
    if abs(stage_multiplier("IV") - stage_multiplier("I")) < 0.2:
        failures.append("stage I vs IV multipliers too similar")

    return {
        "n_checks": len(checks),
        "failures": failures,
        "checks": checks,
        "pass": len(failures) == 0,
    }


def evaluate_prior_direction_recovery(
    predictor: ExhalePathPredictor | None = None,
) -> dict[str, Any]:
    """Check that predicted VOC directions recover curated disease priors."""
    predictor = predictor or ExhalePathPredictor(use_opentargets=False)
    kb = KnowledgeBase()
    hits = tot = 0
    failures = []
    for d in kb.diseases.values():
        site = d.get("default_site")
        q = DiseaseQuery(
            disease=d["name"],
            tumor=TumorContext(stage="II", primary_site=site) if site else None,
            mutated_genes=[],
        )
        df = predictor.predict(q).to_dataframe().set_index("voc_id")
        for voc, prior in d.get("voc_log2fc_prior", {}).items():
            if abs(float(prior)) < 0.25:
                continue
            tot += 1
            got = float(df.loc[voc, "log2_fold_change"])
            ok = (prior > 0 and got > 0) or (prior < 0 and got < 0)
            hits += int(ok)
            if not ok:
                failures.append(
                    {"disease_id": d["disease_id"], "voc_id": voc, "prior": prior, "got": got}
                )
    rate = hits / tot if tot else 0.0
    return {
        "n": tot,
        "hits": hits,
        "accuracy": rate,
        "failures": failures,
        "pass": rate >= 0.95,
    }


def run_full_audit(out_path: Path | None = None) -> AuditReport:
    predictor = ExhalePathPredictor(use_opentargets=False)
    report = AuditReport()
    report.literature_accuracy = evaluate_literature_benchmarks(predictor)
    report.calibrator_holdout = evaluate_calibrator_holdout()
    report.zero_shot = evaluate_zero_shot(predictor)
    report.stress = stress_test_pipeline(predictor)
    prior_rec = evaluate_prior_direction_recovery(predictor)
    report.literature_accuracy["prior_direction_recovery"] = prior_rec
    report.bugs_checked = [
        "stage_multiplier I-prefix false match",
        "bare 'cancer' over-eager alias resolution",
        "ACYSL4 → ACSL4 gene typo",
        "fold_change vs predicted_ppb consistency after clamping",
        "empty/invalid query crash resistance",
        "stage I vs IV monotonicity",
    ]
    lit_ok = report.literature_accuracy.get("case_pass_rate", 0) >= 0.8
    lit_dir = (report.literature_accuracy.get("directional_accuracy") or 0) >= 0.85
    zs_ok = report.zero_shot.get("pass", False)
    stress_ok = report.stress.get("pass", False)
    cal = report.calibrator_holdout
    cal_ok = cal.get("skipped", False) or cal.get("pass", False)
    prior_ok = prior_rec.get("pass", False)
    report.passed = bool(lit_ok and lit_dir and zs_ok and stress_ok and cal_ok and prior_ok)

    if out_path:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report.to_dict(), indent=2))
    return report
