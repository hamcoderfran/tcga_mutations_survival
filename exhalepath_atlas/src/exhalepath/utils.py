from __future__ import annotations

import re


_STAGE_ORDINAL = [
    ("IV", 4.0),
    ("IIIB", 3.5),
    ("IIIA", 3.3),
    ("III", 3.0),
    ("IIB", 2.5),
    ("IIA", 2.3),
    ("II", 2.0),
    ("IB", 1.5),
    ("IA", 1.3),
    ("I", 1.0),
    ("0", 0.0),
]

_STAGE_MULTIPLIER = {
    0.0: 0.85,
    1.0: 1.0,
    1.3: 1.0,
    1.5: 1.05,
    2.0: 1.15,
    2.3: 1.15,
    2.5: 1.2,
    3.0: 1.35,
    3.3: 1.35,
    3.5: 1.4,
    4.0: 1.55,
}


def normalize_stage_token(stage: str | None) -> str | None:
    if not stage or not isinstance(stage, str):
        return None
    s = stage.upper().strip()
    s = re.sub(r"^STAGE\s*", "", s)
    s = s.replace("STAGE", "").strip()
    # Keep leading roman / numeric token (e.g. "IVB", "IIA", "4")
    m = re.match(r"^(IV|III|II|I|0|[0-4])([A-C]?)", s)
    if not m:
        # Arabic numerals
        m2 = re.search(r"\b([0-4])\b", s)
        if m2:
            return m2.group(1)
        return s
    return m.group(1) + (m.group(2) or "")


def stage_ordinal(stage: str | None, default: float = 1.5) -> float:
    """Parse AJCC-like stage strings to an ordinal. Longest roman match wins."""
    token = normalize_stage_token(stage)
    if token is None:
        return default
    # Arabic 1-4
    if token.isdigit():
        return float(token)
    for key, num in _STAGE_ORDINAL:
        if token == key or token.startswith(key):
            return num
    return default


def stage_multiplier(stage: str | None) -> float:
    ord_ = stage_ordinal(stage, default=-1.0)
    if ord_ < 0:
        return 1.0
    # Snap to nearest known multiplier key
    if ord_ in _STAGE_MULTIPLIER:
        return _STAGE_MULTIPLIER[ord_]
    # Interpolate coarsely
    keys = sorted(_STAGE_MULTIPLIER)
    for k in reversed(keys):
        if ord_ >= k:
            return _STAGE_MULTIPLIER[k]
    return 1.0
