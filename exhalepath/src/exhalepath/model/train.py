from __future__ import annotations

import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split

from ..config import MODELS_DIR, PROCESSED_DIR
from ..knowledge.loader import KnowledgeBase
from .features import feature_frame_for_training


def train_calibrator(
    *,
    case_features_path: Path | None = None,
    voc_targets_path: Path | None = None,
    out_dir: Path | None = None,
    test_size: float = 0.2,
    random_state: int = 7,
) -> dict:
    """
    Train a pathway-feature → log2 VOC fold-change calibrator.

    Uses HistGradientBoosting when sample counts are large; falls back to GBR.
    """
    out_dir = Path(out_dir or MODELS_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    case_features_path = Path(case_features_path or PROCESSED_DIR / "case_pathway_features.csv")
    voc_targets_path = Path(voc_targets_path or PROCESSED_DIR / "voc_training_targets.csv")

    if not case_features_path.exists() or not voc_targets_path.exists():
        raise FileNotFoundError(
            "Training corpus missing. Run: python -m exhalepath build-corpus"
        )

    kb = KnowledgeBase()
    case_features = pd.read_csv(case_features_path)
    voc_targets = pd.read_csv(voc_targets_path)

    # Subsample extremely large matrices for tractable training while retaining diversity
    if len(voc_targets) > 400_000:
        voc_targets = voc_targets.sample(400_000, random_state=random_state)

    X, y, voc_ids = feature_frame_for_training(case_features, voc_targets, kb)

    # One model per VOC (specialized heads) for better biochemistry fidelity
    models = {}
    metrics = {}
    for voc_id in sorted(voc_ids.unique()):
        mask = voc_ids == voc_id
        Xi = X.loc[mask]
        yi = y.loc[mask]
        if len(Xi) < 30:
            continue
        Xtr, Xte, ytr, yte = train_test_split(
            Xi, yi, test_size=test_size, random_state=random_state
        )
        if len(Xtr) >= 500:
            model = HistGradientBoostingRegressor(
                max_depth=6,
                learning_rate=0.08,
                max_iter=200,
                random_state=random_state,
            )
        else:
            model = GradientBoostingRegressor(
                max_depth=3,
                learning_rate=0.08,
                n_estimators=120,
                random_state=random_state,
            )
        model.fit(Xtr, ytr)
        pred = model.predict(Xte)
        metrics[voc_id] = {
            "mae_log2fc": float(mean_absolute_error(yte, pred)),
            "r2": float(r2_score(yte, pred)) if len(np.unique(yte)) > 1 else None,
            "n_train": int(len(Xtr)),
            "n_test": int(len(Xte)),
        }
        models[voc_id] = model

    bundle = {
        "models": models,
        "feature_columns": list(X.columns),
        "metrics": metrics,
        "version": "exhalepath-calibrator-0.1",
    }
    model_path = out_dir / "voc_calibrator.joblib"
    metrics_path = out_dir / "voc_calibrator_metrics.json"
    joblib.dump(bundle, model_path)
    metrics_path.write_text(json.dumps(metrics, indent=2))
    print(f"Wrote {model_path}")
    print(json.dumps({"n_voc_models": len(models), "mean_mae": float(np.mean([m['mae_log2fc'] for m in metrics.values()]))}, indent=2))
    return {"model_path": model_path, "metrics_path": metrics_path, "metrics": metrics}
