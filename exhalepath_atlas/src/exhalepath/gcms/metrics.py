"""Patient-level diagnostic metrics for GC-MS / breath VOC research."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)

from .schema import DiagnosticMetrics


def _as_arrays(
    y_true: list[int] | np.ndarray, scores: list[float] | np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(y_true, dtype=int)
    s = np.asarray(scores, dtype=float)
    return y, s


def bootstrap_auroc_ci(
    y_true: np.ndarray,
    scores: np.ndarray,
    *,
    n_boot: int = 500,
    seed: int = 42,
) -> Optional[tuple[float, float]]:
    rng = np.random.default_rng(seed)
    if len(np.unique(y_true)) < 2:
        return None
    vals = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yb, sb = y_true[idx], scores[idx]
        if len(np.unique(yb)) < 2:
            continue
        vals.append(float(roc_auc_score(yb, sb)))
    if len(vals) < 20:
        return None
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def youden_threshold(y_true: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return float(np.median(scores))
    fpr, tpr, thr = roc_curve(y_true, scores)
    j = tpr - fpr
    return float(thr[int(np.argmax(j))])


def binary_operating_point(
    y_true: np.ndarray, scores: np.ndarray, threshold: float
) -> dict[str, Any]:
    pred = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) else None
    spec = tn / (tn + fp) if (tn + fp) else None
    ppv = tp / (tp + fp) if (tp + fp) else None
    npv = tn / (tn + fn) if (tn + fn) else None
    return {
        "sensitivity": sens,
        "specificity": spec,
        "ppv": ppv,
        "npv": npv,
        "confusion": {
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
        },
    }


def compute_diagnostic_metrics(
    y_true: list[int] | np.ndarray,
    scores: list[float] | np.ndarray,
    *,
    n_boot: int = 500,
    seed: int = 42,
    threshold: float | None = None,
) -> DiagnosticMetrics:
    y, s = _as_arrays(y_true, scores)
    notes: list[str] = []
    if len(y) < 4:
        notes.append("n<4 — metrics unstable")
    if len(np.unique(y)) < 2:
        return DiagnosticMetrics(
            n_subjects=int(len(y)),
            n_positive=int((y == 1).sum()),
            n_negative=int((y == 0).sum()),
            notes=notes + ["single class — AUROC undefined"],
        )

    auroc = float(roc_auc_score(y, s))
    try:
        auprc = float(average_precision_score(y, s))
    except Exception:  # noqa: BLE001
        auprc = None
    # map scores to [0,1] for Brier via rank-based sigmoid-ish
    s_min, s_max = float(s.min()), float(s.max())
    if s_max > s_min:
        prob = (s - s_min) / (s_max - s_min)
    else:
        prob = np.full_like(s, 0.5, dtype=float)
    try:
        brier = float(brier_score_loss(y, prob))
    except Exception:  # noqa: BLE001
        brier = None

    thr = youden_threshold(y, s) if threshold is None else float(threshold)
    op = binary_operating_point(y, s, thr)
    ci = bootstrap_auroc_ci(y, s, n_boot=n_boot, seed=seed)

    return DiagnosticMetrics(
        n_subjects=int(len(y)),
        n_positive=int((y == 1).sum()),
        n_negative=int((y == 0).sum()),
        auroc=auroc,
        auroc_ci95=ci,
        auprc=auprc,
        sensitivity=op["sensitivity"],
        specificity=op["specificity"],
        ppv=op["ppv"],
        npv=op["npv"],
        youden_threshold=thr,
        confusion=op["confusion"],
        brier=brier,
        notes=notes,
    )


def confounder_proxy_auroc(
    feature_matrix: np.ndarray,
    confounder: np.ndarray,
) -> Optional[float]:
    """AUROC of a quick logistic probe: can VOCs predict the confounder?

    High values mean the feature space carries strong confounder signal.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler

    y = np.asarray(confounder)
    if len(np.unique(y)) < 2 or len(y) < 8:
        return None
    X = StandardScaler().fit_transform(np.asarray(feature_matrix, dtype=float))
    clf = LogisticRegression(max_iter=500, class_weight="balanced")
    try:
        scores = cross_val_score(clf, X, y, cv=min(5, int(y.sum()), int((1 - y).sum()) or 2), scoring="roc_auc")
        return float(np.nanmean(scores))
    except Exception:  # noqa: BLE001
        return None


def stratified_auroc(
    y_true: list[int] | np.ndarray,
    scores: list[float] | np.ndarray,
    strata: dict[str, list[Any] | np.ndarray],
    *,
    min_n: int = 8,
) -> dict[str, Any]:
    """Compute AUROC within strata (sex, smoking, age tertile, …) when possible.

    strata maps stratum_name → per-subject values aligned with y_true/scores.
    Missing strata or single-class subsets are reported as unavailable.
    """
    y, s = _as_arrays(y_true, scores)
    out: dict[str, Any] = {"overall": None, "strata": {}}
    if len(np.unique(y)) >= 2:
        out["overall"] = float(roc_auc_score(y, s))

    for name, values in strata.items():
        vals = list(values)
        if len(vals) != len(y):
            out["strata"][name] = {"error": "length_mismatch"}
            continue
        # derive age tertiles if numeric age
        series = pd_series_safe(vals)
        if name.lower() == "age" and series is not None:
            try:
                cats = pd_qcut_tertiles(series)
            except Exception:  # noqa: BLE001
                cats = [str(v) if v is not None else "NA" for v in vals]
        else:
            cats = [("NA" if v is None or (isinstance(v, float) and np.isnan(v)) else str(v)) for v in vals]

        bucket: dict[str, dict[str, Any]] = {}
        for level in sorted(set(cats)):
            idx = [i for i, c in enumerate(cats) if c == level]
            if len(idx) < min_n:
                bucket[level] = {"n": len(idx), "auroc": None, "note": f"n<{min_n}"}
                continue
            yb = y[idx]
            sb = s[idx]
            if len(np.unique(yb)) < 2:
                bucket[level] = {"n": len(idx), "auroc": None, "note": "single_class"}
                continue
            bucket[level] = {
                "n": len(idx),
                "n_positive": int((yb == 1).sum()),
                "n_negative": int((yb == 0).sum()),
                "auroc": float(roc_auc_score(yb, sb)),
            }
        out["strata"][name] = bucket
    return out


def pd_series_safe(vals: list[Any]):
    try:
        import pandas as pd

        s = pd.to_numeric(pd.Series(vals), errors="coerce")
        if s.notna().sum() >= 6:
            return s
    except Exception:  # noqa: BLE001
        return None
    return None


def pd_qcut_tertiles(series) -> list[str]:
    import pandas as pd

    cats = pd.qcut(series, q=3, labels=["age_T1", "age_T2", "age_T3"], duplicates="drop")
    return [str(c) if pd.notna(c) else "NA" for c in cats]
