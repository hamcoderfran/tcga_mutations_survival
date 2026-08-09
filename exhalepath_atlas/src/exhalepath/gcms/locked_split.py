"""Preregistration-style locked patient-level splits with content hashing."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.model_selection import StratifiedKFold, LeaveOneOut

from .patient_matrix import PatientVOCMatrix
from .schema import SplitFold, SplitManifest


def _content_hash(payload: dict[str, Any]) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def lock_split(
    matrix: PatientVOCMatrix,
    *,
    strategy: str = "stratified_kfold",
    n_splits: int = 5,
    seed: int = 42,
    out_path: Path | None = None,
) -> SplitManifest:
    """Create an immutable patient-level split manifest.

    Never splits at the peak/feature level — only subject_ids.
    """
    subjects = np.asarray(matrix.subject_ids)
    y = matrix.labels.loc[subjects].to_numpy(dtype=int)
    folds: list[SplitFold] = []

    if strategy == "loocv":
        loo = LeaveOneOut()
        for i, (tr, te) in enumerate(loo.split(subjects)):
            folds.append(
                SplitFold(
                    fold_id=f"loo_{i}",
                    train_subject_ids=subjects[tr].tolist(),
                    test_subject_ids=subjects[te].tolist(),
                )
            )
        n_splits = len(folds)
    elif strategy in {"stratified_kfold", "nested_kfold"}:
        # nested_kfold here locks the *outer* folds; inner HP loops are caller's job
        n_splits = min(n_splits, int(y.sum()), int((1 - y).sum()), len(subjects))
        n_splits = max(2, n_splits)
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        for i, (tr, te) in enumerate(skf.split(subjects, y)):
            folds.append(
                SplitFold(
                    fold_id=f"fold_{i}",
                    train_subject_ids=subjects[tr].tolist(),
                    test_subject_ids=subjects[te].tolist(),
                )
            )
    else:
        raise ValueError(f"Unknown strategy: {strategy}")

    pre_hash = {
        "dataset_id": matrix.study_id,
        "disease_id": matrix.disease_id,
        "seed": seed,
        "strategy": strategy,
        "n_splits": n_splits,
        "subject_ids": subjects.tolist(),
        "labels": y.tolist(),
        "folds": [f.model_dump() for f in folds],
    }
    digest = _content_hash(pre_hash)
    manifest = SplitManifest(
        dataset_id=matrix.study_id,
        seed=seed,
        strategy=strategy,
        n_splits=n_splits,
        subject_ids=subjects.tolist(),
        folds=folds,
        content_sha256=digest,
        notes=[
            "Patient-level only — do not fit scalers/feature selectors on test subjects.",
            f"disease_id={matrix.disease_id}",
            f"n_positive={int((y == 1).sum())} n_negative={int((y == 0).sum())}",
        ],
    )
    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(manifest.model_dump_json(indent=2))
    return manifest


def load_split_manifest(path: Path) -> SplitManifest:
    return SplitManifest.model_validate_json(Path(path).read_text())


def verify_split_manifest(manifest: SplitManifest, matrix: PatientVOCMatrix) -> list[str]:
    """Return list of integrity problems (empty = OK)."""
    problems: list[str] = []
    subjects = set(matrix.subject_ids)
    if set(manifest.subject_ids) != subjects:
        problems.append("subject_ids mismatch vs matrix")
    seen_test: set[str] = set()
    for fold in manifest.folds:
        tr, te = set(fold.train_subject_ids), set(fold.test_subject_ids)
        if tr & te:
            problems.append(f"{fold.fold_id}: train∩test non-empty")
        if not te <= subjects or not tr <= subjects:
            problems.append(f"{fold.fold_id}: unknown subject ids")
        seen_test |= te
    # recompute hash without created_utc
    pre = {
        "dataset_id": manifest.dataset_id,
        "disease_id": matrix.disease_id,
        "seed": manifest.seed,
        "strategy": manifest.strategy,
        "n_splits": manifest.n_splits,
        "subject_ids": manifest.subject_ids,
        "labels": matrix.labels.loc[manifest.subject_ids].astype(int).tolist(),
        "folds": [f.model_dump() for f in manifest.folds],
    }
    digest = _content_hash(pre)
    if digest != manifest.content_sha256:
        problems.append(
            f"content_sha256 mismatch (file={manifest.content_sha256[:12]}… "
            f"recomputed={digest[:12]}…)"
        )
    return problems
