"""Nested sparse / feature-selected logistic models for patient GC-MS diagnostics.

Outer folds evaluate AUROC; inner folds choose features (SelectKBest or L1).
Reports transferable signature metrics separately from fit-on-cohort sparse ML.
"""

from __future__ import annotations

from typing import Any, Literal, Optional

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .locked_split import SplitManifest
from .patient_matrix import PatientVOCMatrix


def nested_sparse_logistic(
    matrix: PatientVOCMatrix,
    manifest: SplitManifest,
    *,
    max_features: int = 8,
    selector: Literal["kbest", "l1"] = "kbest",
    seed: int = 42,
    inner_splits: int = 3,
) -> dict[str, Any]:
    """Outer locked folds; inner FS + logistic. Returns nested AUROC + OOF scores."""
    X_all = matrix.matrix.to_numpy(dtype=float)
    y_all = matrix.labels.to_numpy(dtype=int)
    id_to_i = {sid: i for i, sid in enumerate(matrix.subject_ids)}
    feat_names = list(matrix.matrix.columns.astype(str))

    fold_rows = []
    oof = pd.Series(index=matrix.subject_ids, dtype=float)
    selected_per_fold: list[list[str]] = []

    for fold in manifest.folds:
        tr = [id_to_i[s] for s in fold.train_subject_ids]
        te = [id_to_i[s] for s in fold.test_subject_ids]
        Xtr, ytr = X_all[tr], y_all[tr]
        Xte, yte = X_all[te], y_all[te]
        k = int(min(max_features, Xtr.shape[1], max(2, ytr.sum()), max(2, (1 - ytr).sum())))
        if selector == "l1":
            pipe = Pipeline(
                [
                    ("scaler", StandardScaler()),
                    (
                        "clf",
                        LogisticRegression(
                            solver="saga",
                            l1_ratio=1.0,
                            C=0.5,
                            max_iter=2000,
                            class_weight="balanced",
                            random_state=seed,
                        ),
                    ),
                ]
            )
            # Inner C selection
            best_c, best_score = 0.5, -1.0
            if len(np.unique(ytr)) >= 2 and len(ytr) >= inner_splits * 2:
                inner = StratifiedKFold(
                    n_splits=min(inner_splits, int(ytr.sum()), int((1 - ytr).sum()) or 2),
                    shuffle=True,
                    random_state=seed,
                )
                for C in (0.1, 0.5, 1.0, 2.0):
                    scores = []
                    for itr, iva in inner.split(Xtr, ytr):
                        p = Pipeline(
                            [
                                ("scaler", StandardScaler()),
                                (
                                    "clf",
                                    LogisticRegression(
                                        solver="saga",
                                        l1_ratio=1.0,
                                        C=C,
                                        max_iter=2000,
                                        class_weight="balanced",
                                        random_state=seed,
                                    ),
                                ),
                            ]
                        )
                        p.fit(Xtr[itr], ytr[itr])
                        pr = p.predict_proba(Xtr[iva])[:, 1]
                        if len(np.unique(ytr[iva])) >= 2:
                            scores.append(float(roc_auc_score(ytr[iva], pr)))
                    if scores and float(np.mean(scores)) > best_score:
                        best_score = float(np.mean(scores))
                        best_c = C
            pipe.set_params(clf__C=best_c)
            pipe.fit(Xtr, ytr)
            coef = np.abs(pipe.named_steps["clf"].coef_.ravel())
            keep_idx = np.argsort(coef)[::-1][:k]
            selected = [feat_names[i] for i in keep_idx if coef[i] > 0]
            if len(selected) < 2:
                selected = [feat_names[i] for i in keep_idx[:k]]
            proba = pipe.predict_proba(Xte)[:, 1]
        else:
            # SelectKBest on train only, then logistic
            k = min(k, Xtr.shape[1])
            sel = SelectKBest(f_classif, k=k)
            clf = LogisticRegression(
                max_iter=800, class_weight="balanced", solver="lbfgs", random_state=seed
            )
            pipe = Pipeline([("scaler", StandardScaler()), ("sel", sel), ("clf", clf)])
            pipe.fit(Xtr, ytr)
            mask = pipe.named_steps["sel"].get_support()
            selected = [feat_names[i] for i, m in enumerate(mask) if m]
            proba = pipe.predict_proba(Xte)[:, 1]

        selected_per_fold.append(selected)
        for sid, p in zip(fold.test_subject_ids, proba):
            oof.loc[sid] = float(p)
        auc = (
            float(roc_auc_score(yte, proba))
            if len(np.unique(yte)) >= 2 and len(yte) >= 2
            else None
        )
        fold_rows.append(
            {
                "fold_id": fold.fold_id,
                "n_test": len(te),
                "auroc": auc,
                "n_features": len(selected),
                "features": selected,
                "model": f"nested_sparse_{selector}",
            }
        )

    aucs = [r["auroc"] for r in fold_rows if r["auroc"] is not None]
    # Non-nested optimistic: FS + fit on all
    k_all = int(min(max_features, X_all.shape[1], max(2, int(y_all.sum()))))
    if selector == "l1":
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "clf",
                    LogisticRegression(
                        solver="saga",
                        l1_ratio=1.0,
                        C=0.5,
                        max_iter=2000,
                        class_weight="balanced",
                        random_state=seed,
                    ),
                ),
            ]
        )
    else:
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("sel", SelectKBest(f_classif, k=k_all)),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=800, class_weight="balanced", solver="lbfgs", random_state=seed
                    ),
                ),
            ]
        )
    pipe.fit(X_all, y_all)
    proba_all = pipe.predict_proba(X_all)[:, 1]
    non_nested = (
        float(roc_auc_score(y_all, proba_all)) if len(np.unique(y_all)) >= 2 else None
    )
    nested = float(np.mean(aucs)) if aucs else None
    # consensus features
    from collections import Counter

    cnt = Counter(f for fold in selected_per_fold for f in fold)
    consensus = [f for f, _ in cnt.most_common(max_features)]

    return {
        "selector": selector,
        "max_features": max_features,
        "mean_test_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": (non_nested - nested) if (non_nested is not None and nested is not None) else None,
        "folds": fold_rows,
        "consensus_features": consensus,
        "oof_scores": oof,
    }


def learning_curve_nested_auroc(
    matrix: PatientVOCMatrix,
    *,
    sizes: list[int] | None = None,
    n_repeats: int = 20,
    seed: int = 42,
    max_features: int = 8,
) -> dict[str, Any]:
    """Subsample subjects and estimate nested-ish AUROC (JBR n-curve style)."""
    rng = np.random.default_rng(seed)
    n = len(matrix.subject_ids)
    sizes = sizes or [s for s in (12, 16, 20, 24, 28, 32, n) if s <= n]
    y = matrix.labels.to_numpy(dtype=int)
    X = matrix.matrix.to_numpy(dtype=float)
    rows = []
    for size in sizes:
        vals = []
        for _ in range(n_repeats):
            # stratified subsample
            pos = np.where(y == 1)[0]
            neg = np.where(y == 0)[0]
            n_pos = max(3, int(round(size * (y.mean()))))
            n_neg = size - n_pos
            if n_pos > len(pos) or n_neg > len(neg) or n_pos < 3 or n_neg < 3:
                continue
            idx = np.concatenate(
                [
                    rng.choice(pos, n_pos, replace=False),
                    rng.choice(neg, n_neg, replace=False),
                ]
            )
            rng.shuffle(idx)
            Xs, ys = X[idx], y[idx]
            # 3-fold stratified
            try:
                skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=int(rng.integers(0, 1e9)))
                fold_aucs = []
                for tr, te in skf.split(Xs, ys):
                    k = min(max_features, Xs.shape[1], max(2, int(ys[tr].sum())))
                    pipe = Pipeline(
                        [
                            ("scaler", StandardScaler()),
                            ("sel", SelectKBest(f_classif, k=k)),
                            (
                                "clf",
                                LogisticRegression(
                                    max_iter=600, class_weight="balanced", solver="lbfgs"
                                ),
                            ),
                        ]
                    )
                    pipe.fit(Xs[tr], ys[tr])
                    pr = pipe.predict_proba(Xs[te])[:, 1]
                    if len(np.unique(ys[te])) >= 2:
                        fold_aucs.append(float(roc_auc_score(ys[te], pr)))
                if fold_aucs:
                    vals.append(float(np.mean(fold_aucs)))
            except Exception:  # noqa: BLE001
                continue
        rows.append(
            {
                "n": size,
                "mean_auroc": float(np.mean(vals)) if vals else None,
                "std_auroc": float(np.std(vals)) if vals else None,
                "n_repeats_ok": len(vals),
            }
        )
    return {
        "description": "Learning curve: nested-ish sparse logistic vs subsample size (JBR-style)",
        "rows": rows,
        "note": "Gains typically flatten near n≈50 in breath ML; n=35 remains wide-CI territory.",
    }


__all__ = ["learning_curve_nested_auroc", "nested_sparse_logistic"]
