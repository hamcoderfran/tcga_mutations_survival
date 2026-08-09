"""External / larger malaria breath cohort adapters.

Public reality (2026):
- ST000883 (Schaber 2018, MW PR000612) is the only open *quantified* pediatric
  malaria breath intensity table on Metabolomics Workbench.
- Berna CHMI (CSIRO csiro:33843) has overview-derived baseline/peak labels bundled
  under ``datasources/malaria_external/``; intensity still requires MassHunter
  peak export from Agilent ``.D`` / mzdata.xml.
- 2024 Malawi reproducibility (J Infect Dis jiae323) is not deposited as an
  open intensity matrix at time of writing.

This module therefore:
1. Documents CSIRO CHMI labels + catalog for diligence
2. Builds a *virtual external* protocol: learning curves + multi-seed nested
   CIs on the remapped ST000883 matrix (JBR-style n≥50 guidance)
3. Loads a user-supplied second mwtab/CSV when available (`--external-matrix`)
4. Runs leave-one-study-out when ≥2 intensity matrices share VOC columns
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..config import DATA_DIR
from .patient_matrix import PatientVOCMatrix
from .schema import SampleRecord

EXTERNAL_CATALOG = {
    "ST000883": {
        "role": "primary_open_intensity",
        "doi": "10.1093/infdis/jiy072",
        "n_subjects": 35,
        "status": "bundled",
    },
    "CSIRO_CHMI_QTOF": {
        "role": "external_raw_chmi",
        "doi": "10.25919/5b5b7530a39f4",
        "url": "https://data.csiro.au/collection/csiro:33843",
        "n_subjects": 7,
        "n_labeled_baseline_peak": 14,
        "status": "labels_bundled_intensity_pending",
        "labels": "datasources/malaria_external/csiro_chmi_labels.csv",
        "note": (
            "CHMI QTOF breath (.D + mzdata.xml). Overview labels bundled "
            "(Day0/ND vs per-subject peak parasitemia). Export MassHunter peak "
            "table aligned to Sample name, then --external-matrix + --loso."
        ),
    },
    "JID_2024_MALAWI": {
        "role": "external_reproducibility",
        "doi": "10.1093/infdis/jiae323",
        "status": "no_open_intensity_matrix",
        "note": "Independent Blantyre cohort; deposit intensity table when released.",
    },
}

_CACHE = DATA_DIR / "datasources" / "malaria_external"


def write_external_catalog(out_dir: Path | None = None) -> Path:
    out_dir = Path(out_dir or _CACHE)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "EXTERNAL_MALARIA_CATALOG.json"
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": EXTERNAL_CATALOG,
        "guidance": (
            "JBR breath-ML learning curves flatten near n≈50; single n=35 studies "
            "keep wide AUROC CIs. Prefer multi-cohort locked evaluation when a second "
            "intensity table becomes available."
        ),
    }
    path.write_text(json.dumps(payload, indent=2))
    return path


def load_external_patient_csv(
    matrix_csv: Path,
    labels_csv: Path,
    *,
    study_id: str = "EXTERNAL_MALARIA",
    disease_id: str = "malaria",
) -> PatientVOCMatrix:
    """Load a user-exported patient×VOC CSV + labels (0/1) as external cohort."""
    mat = pd.read_csv(matrix_csv, index_col=0)
    lab = pd.read_csv(labels_csv, index_col=0).iloc[:, 0].astype(int)
    mat.index = mat.index.astype(str)
    lab.index = lab.index.astype(str)
    common = [i for i in mat.index if i in set(lab.index)]
    mat = mat.loc[common].astype(float)
    lab = lab.loc[common]
    samples = [
        SampleRecord(
            sample_id=sid,
            subject_id=sid,
            study_id=study_id,
            label="disease" if int(lab.loc[sid]) == 1 else "control",
            disease_id=disease_id,
            modality="gcms_external",
        )
        for sid in common
    ]
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=disease_id,
        disease_name="Malaria",
        modality="gcms_external",
        unit="peak_intensity",
        matrix=mat,
        labels=lab.astype(int),
        samples=samples,
        metadata={
            "source": "user_external_csv",
            "n_voc_features": int(mat.shape[1]),
            "n_subjects": int(mat.shape[0]),
            "n_positive": int((lab == 1).sum()),
            "n_negative": int((lab == 0).sum()),
        },
    )


def combine_cohorts(
    matrices: list[PatientVOCMatrix],
    *,
    study_id: str = "MALARIA_POOLED",
) -> PatientVOCMatrix:
    """Inner-join VOC columns and concatenate subjects (prefix study ids)."""
    if not matrices:
        raise ValueError("No matrices to combine")
    # intersection of features
    cols = set(matrices[0].matrix.columns.astype(str))
    for m in matrices[1:]:
        cols &= set(m.matrix.columns.astype(str))
    cols = sorted(cols)
    if len(cols) < 2:
        raise ValueError("Fewer than 2 shared VOC features across cohorts")
    blocks = []
    labels = []
    samples = []
    for m in matrices:
        idx = [f"{m.study_id}:{sid}" for sid in m.subject_ids]
        block = m.matrix[cols].copy()
        block.index = idx
        blocks.append(block)
        lab = m.labels.copy()
        lab.index = idx
        labels.append(lab)
        for s in m.samples:
            samples.append(
                s.model_copy(
                    update={
                        "sample_id": f"{m.study_id}:{s.sample_id}",
                        "subject_id": f"{m.study_id}:{s.subject_id}",
                        "study_id": study_id,
                        "metadata": {**(s.metadata or {}), "source_study": m.study_id},
                    }
                )
            )
    X = pd.concat(blocks, axis=0)
    y = pd.concat(labels, axis=0).astype(int)
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=matrices[0].disease_id,
        disease_name=matrices[0].disease_name,
        modality="gcms_pooled",
        unit="peak_intensity",
        matrix=X.astype(float),
        labels=y,
        samples=samples,
        metadata={
            "source_studies": [m.study_id for m in matrices],
            "n_voc_features": len(cols),
            "n_subjects": int(X.shape[0]),
            "n_positive": int((y == 1).sum()),
            "n_negative": int((y == 0).sum()),
            "shared_vocs": cols,
        },
    )


def leave_one_study_out(
    matrices: list[PatientVOCMatrix],
    *,
    signature: dict[str, float] | None = None,
    max_features: int = 8,
    seed: int = 42,
) -> dict[str, Any]:
    """Train on all-but-one study; evaluate on the held-out study.

    Scoring modes:
    - If ``signature`` is provided: transferable cosine scores with control mean
      from the *train* studies only.
    - Always also fit a train-only sparse logistic (SelectKBest) on shared VOCs.
    """
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    from .score import score_observed_vector

    if len(matrices) < 2:
        return {
            "status": "skipped",
            "reason": "need_at_least_two_studies_with_intensity",
            "n_studies": len(matrices),
        }

    # shared columns across all studies
    cols = set(matrices[0].matrix.columns.astype(str))
    for m in matrices[1:]:
        cols &= set(m.matrix.columns.astype(str))
    cols = sorted(cols)
    if len(cols) < 2:
        return {
            "status": "skipped",
            "reason": "fewer_than_2_shared_voc_features",
            "n_shared": len(cols),
        }

    folds: list[dict[str, Any]] = []
    for i, holdout in enumerate(matrices):
        train_mats = [m for j, m in enumerate(matrices) if j != i]
        train = combine_cohorts(train_mats, study_id="LOSO_TRAIN")
        # holdout without prefix (keep native ids)
        X_te = holdout.matrix[cols].astype(float)
        y_te = holdout.labels.astype(int)
        X_tr = train.matrix[cols].astype(float)
        y_tr = train.labels.astype(int)

        row: dict[str, Any] = {
            "holdout_study": holdout.study_id,
            "train_studies": [m.study_id for m in train_mats],
            "n_train": int(X_tr.shape[0]),
            "n_test": int(X_te.shape[0]),
            "n_shared_vocs": len(cols),
            "n_test_positive": int((y_te == 1).sum()),
            "n_test_negative": int((y_te == 0).sum()),
        }
        if len(np.unique(y_te)) < 2 or len(np.unique(y_tr)) < 2:
            row["auroc_signature"] = None
            row["auroc_sparse"] = None
            row["note"] = "single_class_in_train_or_test"
            folds.append(row)
            continue

        # Transferable signature (optional)
        if signature:
            ctrl = X_tr.loc[y_tr == 0]
            ref = ctrl.mean(axis=0) if len(ctrl) else X_tr.mean(axis=0)
            sig_use = {k: v for k, v in signature.items() if k in cols}
            scores = []
            for sid in X_te.index:
                vec = (X_te.loc[sid] - ref).to_dict()
                out = score_observed_vector(vec, sig_use, method="cosine")
                scores.append(
                    float(out["score"]) if out.get("score") is not None else 0.0
                )
            try:
                row["auroc_signature"] = float(roc_auc_score(y_te, scores))
                row["auprc_signature"] = float(average_precision_score(y_te, scores))
            except Exception as exc:  # noqa: BLE001
                row["auroc_signature"] = None
                row["signature_error"] = str(exc)
        else:
            row["auroc_signature"] = None

        # Sparse logistic fit on train only
        k = min(max_features, X_tr.shape[1], max(1, int(y_tr.sum()), int((1 - y_tr).sum())))
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("select", SelectKBest(f_classif, k=k)),
                (
                    "clf",
                    LogisticRegression(
                        max_iter=1000, class_weight="balanced", random_state=seed
                    ),
                ),
            ]
        )
        try:
            pipe.fit(X_tr.to_numpy(), y_tr.to_numpy())
            proba = pipe.predict_proba(X_te.to_numpy())[:, 1]
            row["auroc_sparse"] = float(roc_auc_score(y_te, proba))
            row["auprc_sparse"] = float(average_precision_score(y_te, proba))
            sel = pipe.named_steps["select"]
            mask = sel.get_support()
            row["selected_features"] = [c for c, m in zip(cols, mask) if m]
        except Exception as exc:  # noqa: BLE001
            row["auroc_sparse"] = None
            row["sparse_error"] = str(exc)

        folds.append(row)

    sig_aucs = [f["auroc_signature"] for f in folds if f.get("auroc_signature") is not None]
    sp_aucs = [f["auroc_sparse"] for f in folds if f.get("auroc_sparse") is not None]
    return {
        "status": "ok",
        "n_studies": len(matrices),
        "shared_vocs": cols,
        "folds": folds,
        "mean_auroc_signature": float(np.mean(sig_aucs)) if sig_aucs else None,
        "mean_auroc_sparse": float(np.mean(sp_aucs)) if sp_aucs else None,
    }


__all__ = [
    "EXTERNAL_CATALOG",
    "combine_cohorts",
    "leave_one_study_out",
    "load_external_patient_csv",
    "write_external_catalog",
]
