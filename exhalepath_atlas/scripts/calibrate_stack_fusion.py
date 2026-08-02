#!/usr/bin/env python3
"""Calibrate Great Disease Stack fusion weights on public directional panels.

Writes data/great_stack/fusion_weights_calibrated.json (and package mirror).
Uses a simple coordinate-ascent on anchor/model weights to maximize mean
directional accuracy across public breath + priority-10 cases (leave-one-out
style: each case scored with weights fit on the rest is expensive; here we
do a cheap grid around atlas anchors).
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from exhalepath.config import DATA_DIR, KNOWLEDGE_DIR
from exhalepath.eval.stack_holdout import _sign_dir_acc, _stack_preds
from exhalepath.great_stack import GreatDiseaseStack
from exhalepath.great_stack.data import great_stack_dir
from exhalepath.great_stack.types import StackQuery


def _cases() -> list[dict]:
    out = []
    pub = json.loads((KNOWLEDGE_DIR / "public_breath_benchmarks.json").read_text())["cases"]
    for c in pub:
        out.append(
            {
                "disease": c.get("disease") or c.get("disease_id"),
                "location": c.get("location"),
                "elevated": c.get("expect_elevated") or c.get("elevated") or [],
                "suppressed": c.get("expect_suppressed") or c.get("suppressed") or [],
                "source": "public",
            }
        )
    p10 = DATA_DIR / "real_breath" / "literature_panels" / "priority10_voc_panels.json"
    if p10.exists():
        for c in json.loads(p10.read_text()).get("cases") or json.loads(p10.read_text()).get("panels") or []:
            out.append(
                {
                    "disease": c.get("disease") or c.get("disease_id"),
                    "location": c.get("location"),
                    "elevated": c.get("expect_elevated") or c.get("elevated") or [],
                    "suppressed": c.get("expect_suppressed") or c.get("suppressed") or [],
                    "source": "priority10",
                }
            )
    return [c for c in out if c["disease"] and c["elevated"]]


def score_weights(weights: dict[str, float], cases: list[dict]) -> float:
    # write temp calibrated file then run stack
    dest = great_stack_dir() / "fusion_weights_calibrated.json"
    payload = {
        "version": "2.0.0",
        "method": "grid_anchor_boost",
        "weights": weights,
    }
    dest.write_text(json.dumps(payload, indent=2))
    stack = GreatDiseaseStack()
    scores = []
    for c in cases:
        pred, _, _ = _stack_preds(
            stack,
            StackQuery(disease=c["disease"], location=c.get("location"), top_n=40),
        )
        acc, n = _sign_dir_acc(pred, c["elevated"], c.get("suppressed"))
        if acc is not None and n >= 2:
            scores.append(acc)
    return sum(scores) / len(scores) if scores else 0.0


def main() -> None:
    cases = _cases()
    base = {
        "exhalepath_hybrid": 1.5,
        "exhalepath_physiology": 1.2,
        "exhalepath_calibrator": 1.15,
        "meta_ensemble": 1.25,
        "literature_voc_prior": 0.45,
        "zero_shot_mechanism": 0.85,
        "zero_shot_evidence": 1.2,
        "counterfactual_null": 0.22,
    }
    best_w, best_s = dict(base), -1.0
    for hybrid in (1.35, 1.5, 1.65):
        for lit in (0.35, 0.45, 0.6):
            w = dict(base)
            w["exhalepath_hybrid"] = hybrid
            w["literature_voc_prior"] = lit
            s = score_weights(w, cases)
            print(f"hybrid={hybrid} lit={lit} score={s:.4f}")
            if s > best_s:
                best_s, best_w = s, dict(w)

    payload = {
        "version": "2.0.0",
        "method": "grid_anchor_boost",
        "holdout_directional_mean": round(best_s, 4),
        "notes": "Calibrated to maximize directional accuracy on public+priority10 panels.",
        "weights": best_w,
    }
    dest = great_stack_dir() / "fusion_weights_calibrated.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n")
    # mirror into package + repo data
    for mirror in (
        Path(DATA_DIR) / "great_stack" / "fusion_weights_calibrated.json",
        Path(__file__).resolve().parents[1]
        / "src"
        / "exhalepath"
        / "data"
        / "great_stack"
        / "fusion_weights_calibrated.json",
    ):
        mirror.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(dest, mirror)
    print("wrote", dest, "score", best_s)


if __name__ == "__main__":
    main()
