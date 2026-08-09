"""Confounder residualization — disease signal after removing age/sex/smoking.

Field gap: most breath-ML papers report disease AUROC without showing what
remains after demographic regression. This module residualizes each VOC
column on covariates (train-fit only when a split is provided).
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from .patient_matrix import PatientVOCMatrix


def _design_matrix(covariates: pd.DataFrame, columns: list[str]) -> tuple[np.ndarray, list[str]]:
    """Build numeric design matrix from selected covariate columns."""
    pieces = []
    names: list[str] = []
    for col in columns:
        if col not in covariates.columns:
            continue
        s = covariates[col]
        if col.lower() in {"sex", "gender"}:
            male = s.astype(str).str.lower().isin(["male", "m"]).astype(float)
            pieces.append(male.to_numpy().reshape(-1, 1))
            names.append("sex_male")
        elif col.lower() in {"smoking", "smoking_status"}:
            # current smoker = 1; else 0 (ex/never/missing → 0)
            cur = s.astype(str).str.lower().str.contains("current").astype(float)
            pieces.append(cur.to_numpy().reshape(-1, 1))
            names.append("smoking_current")
        else:
            v = pd.to_numeric(s, errors="coerce")
            med = float(v.median()) if v.notna().any() else 0.0
            v = v.fillna(med)
            pieces.append(v.to_numpy().reshape(-1, 1))
            names.append(col)
    if not pieces:
        raise ValueError("No usable covariates for residualization")
    X = np.hstack(pieces)
    # drop all-NaN / constant columns
    keep = []
    keep_names = []
    for j, name in enumerate(names):
        col = X[:, j]
        if np.nanstd(col) > 1e-12:
            keep.append(j)
            keep_names.append(name)
    if not keep:
        raise ValueError("Covariates are constant; cannot residualize")
    return X[:, keep], keep_names


def residualize_matrix(
    feature_matrix: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    covariate_cols: list[str] | None = None,
    train_index: list[str] | None = None,
) -> dict[str, Any]:
    """Residualize each feature on covariates; fit on ``train_index`` when given.

    Returns residual DataFrame aligned to ``feature_matrix`` index plus metadata.
    """
    cov_cols = covariate_cols or [
        c for c in ("age", "sex", "smoking", "smoking_status", "bmi") if c in covariates.columns
    ]
    common = [i for i in feature_matrix.index.astype(str) if i in set(covariates.index.astype(str))]
    if len(common) < 8:
        raise ValueError(f"Need ≥8 overlapping samples; got {len(common)}")
    X = feature_matrix.loc[common].astype(float)
    cov = covariates.loc[common]
    Z, z_names = _design_matrix(cov, cov_cols)
    fit_ids = [i for i in (train_index or common) if i in set(common)]
    fit_mask = np.array([i in set(fit_ids) for i in common])
    if fit_mask.sum() < 5:
        raise ValueError("Too few train rows for residualization fit")

    resid = pd.DataFrame(index=common, columns=X.columns, dtype=float)
    r2: dict[str, float] = {}
    for col in X.columns:
        y = X[col].to_numpy(dtype=float)
        model = LinearRegression()
        model.fit(Z[fit_mask], y[fit_mask])
        pred = model.predict(Z)
        resid[col] = y - pred
        ss_tot = float(np.sum((y[fit_mask] - y[fit_mask].mean()) ** 2)) + 1e-12
        ss_res = float(np.sum((y[fit_mask] - model.predict(Z[fit_mask])) ** 2))
        r2[str(col)] = max(0.0, 1.0 - ss_res / ss_tot)

    mean_r2 = float(np.mean(list(r2.values()))) if r2 else None
    return {
        "residual_matrix": resid,
        "covariate_names": z_names,
        "per_voc_r2": r2,
        "mean_covariate_r2": mean_r2,
        "n_samples": len(common),
        "n_train_fit": int(fit_mask.sum()),
        "note": (
            "Residuals = VOC − E[VOC|covariates]. High mean R² means demographics "
            "explain much of the VOC variance; disease models should use residuals."
        ),
    }


def residualize_patient_matrix(
    matrix: PatientVOCMatrix,
    covariates: pd.DataFrame,
    *,
    covariate_cols: list[str] | None = None,
) -> PatientVOCMatrix:
    """Return a copy of ``matrix`` with residualized intensities."""
    cov = covariates.copy()
    if cov.index.name != "sample_id" and "sample_id" in cov.columns:
        cov = cov.set_index("sample_id")
    cov.index = cov.index.astype(str)
    out = residualize_matrix(
        matrix.matrix, cov, covariate_cols=covariate_cols
    )
    resid = out["residual_matrix"]
    return PatientVOCMatrix(
        study_id=matrix.study_id + "_resid",
        disease_id=matrix.disease_id,
        disease_name=matrix.disease_name,
        modality=matrix.modality,
        unit=f"residual({matrix.unit})",
        matrix=resid.astype(float),
        labels=matrix.labels.loc[resid.index].astype(int),
        samples=[s for s in matrix.samples if s.subject_id in set(resid.index)],
        metadata={
            **matrix.metadata,
            "residualized": True,
            "covariate_names": out["covariate_names"],
            "mean_covariate_r2": out["mean_covariate_r2"],
        },
    )


__all__ = ["residualize_matrix", "residualize_patient_matrix"]
