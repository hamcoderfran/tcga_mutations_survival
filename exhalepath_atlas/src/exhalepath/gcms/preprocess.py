"""Blank-ratio / on-breath filters for real GC-MS peak tables."""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd


_BLANK_PAT = re.compile(r"(blank|room.?air|background|ambient|system.?blank)", re.I)


def detect_blank_sample_ids(sample_ids: list[str]) -> list[str]:
    return [s for s in sample_ids if _BLANK_PAT.search(str(s))]


def blank_ratio_filter(
    intensity: pd.DataFrame,
    *,
    blank_ids: list[str] | None = None,
    min_ratio: float = 2.0,
    min_detect_frac: float = 0.5,
) -> dict[str, Any]:
    """Filter features by median(breath)/median(blank) and detection rate.

    Parameters
    ----------
    intensity
        index=sample_id, columns=feature_id (or voc_id)
    blank_ids
        sample ids treated as blanks/background. Auto-detected if None.
    min_ratio
        keep features with median(breath) / median(blank) >= min_ratio
        (or keep all if no blanks)
    min_detect_frac
        among breath samples, fraction with intensity > 0 required to keep feature
    """
    blank_ids = blank_ids or detect_blank_sample_ids(list(intensity.index.astype(str)))
    blank_ids = [b for b in blank_ids if b in set(intensity.index.astype(str))]
    breath_ids = [i for i in intensity.index.astype(str) if i not in set(blank_ids)]
    if not breath_ids:
        return {
            "kept_features": list(intensity.columns.astype(str)),
            "dropped_features": [],
            "n_blanks": 0,
            "ratios": {},
            "note": "no breath samples",
        }

    X = intensity.copy()
    X.index = X.index.astype(str)
    breath = X.loc[breath_ids]
    ratios: dict[str, float] = {}
    kept: list[str] = []
    dropped: list[str] = []

    if blank_ids:
        blanks = X.loc[blank_ids]
        blank_med = blanks.median(axis=0).replace(0, np.nan)
    else:
        blank_med = None

    for feat in X.columns:
        b = breath[feat].to_numpy(dtype=float)
        detect = float(np.mean(b > 0)) if len(b) else 0.0
        if detect < min_detect_frac:
            dropped.append(str(feat))
            ratios[str(feat)] = float("nan")
            continue
        if blank_med is None:
            kept.append(str(feat))
            ratios[str(feat)] = float("inf")
            continue
        bm = float(blank_med.get(feat, np.nan))
        breath_med = float(np.nanmedian(b))
        if not np.isfinite(bm) or bm <= 0:
            ratio = float("inf") if breath_med > 0 else 0.0
        else:
            ratio = breath_med / bm
        ratios[str(feat)] = ratio
        if ratio >= min_ratio:
            kept.append(str(feat))
        else:
            dropped.append(str(feat))

    return {
        "kept_features": kept,
        "dropped_features": dropped,
        "n_blanks": len(blank_ids),
        "blank_ids": blank_ids,
        "min_ratio": min_ratio,
        "min_detect_frac": min_detect_frac,
        "ratios": ratios,
        "note": (
            "no blank samples detected — detection-fraction filter only"
            if not blank_ids
            else f"blank-ratio filter min_ratio={min_ratio}"
        ),
    }


def apply_feature_filter(
    intensity: pd.DataFrame, kept_features: list[str]
) -> pd.DataFrame:
    cols = [c for c in intensity.columns if str(c) in set(kept_features)]
    return intensity.loc[:, cols]
