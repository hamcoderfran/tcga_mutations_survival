"""Score observed patient VOC vectors against ExhalePath / stack disease signatures.

This turns the mechanism stack into a research diagnostic *template matcher*:
predicted disease log2fc directions become a signature; patient intensities
(relative to control cohort mean, or z-scored) are correlated / dotted with it.
"""

from __future__ import annotations

from typing import Any, Literal, Optional

import numpy as np
import pandas as pd

from .patient_matrix import PatientVOCMatrix


SignatureSource = Literal["hybrid", "stack", "literature", "data_loo"]


def disease_signature(
    disease_id: str,
    *,
    source: SignatureSource = "hybrid",
    top_n: int = 40,
    voc_ids: list[str] | None = None,
) -> dict[str, float]:
    """Return voc_id → expected log2fc (disease vs control direction)."""
    if source == "hybrid":
        from ..biomarker import ExhaleBiomarkerEngine

        eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=False)
        report = eng.predict(disease_id, top_n=top_n, explain=False)
        sig = {p.voc_id: float(p.log2_fold_change) for p in report.result.bundle.predictions}
    elif source == "stack":
        from ..great_stack import GreatDiseaseStack
        from ..great_stack.types import StackQuery

        stack = GreatDiseaseStack()
        result = stack.predict(StackQuery(disease=disease_id, top_n=top_n))
        sig = {v.voc_id: float(v.fused_log2fc) for v in result.fused_vocs}
        # pad with hybrid votes for coverage
        for mo in result.model_outputs:
            if mo.model_id == "exhalepath_hybrid":
                for s in mo.voc_signals:
                    sig.setdefault(s.voc_id, float(s.log2fc))
    elif source == "literature":
        from ..config import DATA_DIR
        import json

        path = DATA_DIR / "real_breath" / "literature_panels" / "priority10_voc_panels.json"
        payload = json.loads(path.read_text())
        sig = {}
        for p in payload.get("panels") or []:
            if p.get("disease_id") == disease_id:
                for voc, fc in (p.get("measured_log2fc") or {}).items():
                    sig[str(voc)] = float(fc)
                break
        if not sig:
            raise ValueError(f"No literature panel for {disease_id}")
    else:
        raise ValueError(f"Unknown signature source: {source}")

    if voc_ids is not None:
        sig = {k: v for k, v in sig.items() if k in set(voc_ids)}
    return sig


def _patient_relative_profile(
    matrix: PatientVOCMatrix,
    *,
    reference: Literal["control_mean", "cohort_mean"] = "control_mean",
) -> pd.DataFrame:
    """Patient log2(intensity / reference) per VOC — comparable to log2fc signature."""
    X = matrix.matrix.clip(lower=1e-12)
    if reference == "control_mean":
        ctrl = matrix.labels == 0
        if ctrl.sum() < 2:
            ref = X.mean(axis=0)
        else:
            ref = X.loc[ctrl].mean(axis=0).replace(0, np.nan)
            ref = ref.fillna(X.mean(axis=0))
    else:
        ref = X.mean(axis=0).replace(0, np.nan).fillna(1.0)
    # avoid divide-by-zero
    ref = ref.clip(lower=1e-12)
    return np.log2(X.div(ref, axis=1))


def score_patients(
    matrix: PatientVOCMatrix,
    signature: dict[str, float],
    *,
    method: Literal["cosine", "dot", "spearman"] = "cosine",
    reference: Literal["control_mean", "cohort_mean"] = "control_mean",
) -> pd.Series:
    """Score each patient; higher ⇒ more disease-like vs signature."""
    profiles = _patient_relative_profile(matrix, reference=reference)
    common = [c for c in profiles.columns if c in signature and abs(signature[c]) > 1e-9]
    if len(common) < 2:
        raise RuntimeError(
            f"Fewer than 2 overlapping VOCs between matrix and signature "
            f"(overlap={len(common)}). Mapped atlas VOCs may be sparse for this study."
        )
    P = profiles[common].to_numpy(dtype=float)
    s = np.asarray([signature[c] for c in common], dtype=float)

    scores = []
    for i in range(P.shape[0]):
        x = P[i]
        if method == "dot":
            scores.append(float(np.dot(x, s)))
        elif method == "cosine":
            nx, ns = np.linalg.norm(x), np.linalg.norm(s)
            scores.append(float(np.dot(x, s) / (nx * ns + 1e-12)))
        elif method == "spearman":
            # rank correlation
            rx = pd.Series(x).rank().to_numpy()
            rs = pd.Series(s).rank().to_numpy()
            rx = (rx - rx.mean()) / (rx.std() + 1e-12)
            rs = (rs - rs.mean()) / (rs.std() + 1e-12)
            scores.append(float(np.mean(rx * rs)))
        else:
            raise ValueError(method)
    return pd.Series(scores, index=profiles.index, name=f"score_{method}")


def score_observed_vector(
    observed: dict[str, float],
    signature: dict[str, float],
    *,
    method: Literal["cosine", "dot", "spearman"] = "cosine",
) -> dict[str, Any]:
    """Score a single observed VOC intensity/log2fc vector against a signature."""
    common = [k for k in observed if k in signature]
    if len(common) < 2:
        return {"score": None, "n_overlap": len(common), "method": method}
    x = np.asarray([float(observed[k]) for k in common])
    s = np.asarray([float(signature[k]) for k in common])
    if method == "dot":
        score = float(np.dot(x, s))
    elif method == "cosine":
        score = float(np.dot(x, s) / (np.linalg.norm(x) * np.linalg.norm(s) + 1e-12))
    else:
        rx = pd.Series(x).rank().to_numpy()
        rs = pd.Series(s).rank().to_numpy()
        rx = (rx - rx.mean()) / (rx.std() + 1e-12)
        rs = (rs - rs.mean()) / (rs.std() + 1e-12)
        score = float(np.mean(rx * rs))
    return {
        "score": score,
        "n_overlap": len(common),
        "method": method,
        "vocs_used": common,
    }
