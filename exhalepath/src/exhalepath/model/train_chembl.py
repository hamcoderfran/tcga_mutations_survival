"""
Train a chemogenomic auxiliary model on ChEMBL activity rows.

Supervised task
---------------
  X: smiles_len, standard_type encodings, pathway multi-hot, gene frequency
  y: pchembl_value

This learns which chemical / assay contexts potently modulate VOC-pathway
enzymes. Distilled pathway ligandability is written back for ExhalePath
predict-time features. It does **not** replace breath-ppb calibration.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split

from ..config import CHEMBL_DIR, KNOWLEDGE_DIR, MODELS_DIR
from ..ingest.chembl_harvest import distill_chembl_priors
from ..knowledge.loader import KnowledgeBase, clear_knowledge_cache


def _design_matrix(df: pd.DataFrame, pathway_ids: list[str]) -> tuple[pd.DataFrame, pd.Series]:
    d = df.dropna(subset=["pchembl_value"]).copy()
    d["pchembl_value"] = d["pchembl_value"].astype(float)
    d["smiles_len"] = d.get("smiles_len", pd.Series(0, index=d.index)).fillna(0).astype(float)
    # standard type one-hots
    for st in ["IC50", "Ki", "Kd", "EC50", "AC50", "Potency"]:
        d[f"std_{st}"] = (d.get("standard_type", pd.Series("", index=d.index)).astype(str) == st).astype(
            float
        )
    # pathway multi-hot
    pw_series = d.get("pathways", pd.Series("", index=d.index)).astype(str)
    for pid in pathway_ids:
        pat = rf"(?:^|\|){re.escape(pid)}(?:\||$)"
        d[f"pw_{pid}"] = pw_series.str.contains(pat, regex=True).astype(float)
    # gene frequency encoding
    gene_freq = d["gene"].astype(str).map(d["gene"].astype(str).value_counts()).astype(float)
    d["gene_freq"] = np.log1p(gene_freq)
    feat_cols = (
        ["smiles_len", "gene_freq"]
        + [f"std_{st}" for st in ["IC50", "Ki", "Kd", "EC50", "AC50", "Potency"]]
        + [f"pw_{pid}" for pid in pathway_ids]
    )
    X = d[feat_cols].fillna(0.0)
    y = d["pchembl_value"]
    return X, y


def train_chembl_aux_model(
    *,
    activities_path: Path | None = None,
    out_dir: Path | None = None,
    max_rows: int | None = 2_000_000,
    test_size: float = 0.1,
    random_state: int = 7,
    update_knowledge_priors: bool = True,
) -> dict:
    out_dir = Path(out_dir or MODELS_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    activities_path = Path(activities_path or CHEMBL_DIR / "chembl_pathway_activities.csv")
    if not activities_path.exists():
        raise FileNotFoundError(
            f"Missing {activities_path}. Run: python -m exhalepath harvest-chembl"
        )

    kb = KnowledgeBase()
    df = pd.read_csv(activities_path)
    if max_rows is not None and len(df) > max_rows:
        df = df.sample(max_rows, random_state=random_state)

    pathway_ids = list(kb.pathways.keys())
    X, y = _design_matrix(df, pathway_ids)
    if len(X) < 100:
        raise RuntimeError(f"Need ≥100 ChEMBL rows with pChEMBL; got {len(X)}")

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=test_size, random_state=random_state)
    model = HistGradientBoostingRegressor(
        max_depth=8,
        learning_rate=0.06,
        max_iter=250,
        random_state=random_state,
    )
    model.fit(Xtr, ytr)
    pred = model.predict(Xte)
    metrics = {
        "n_rows": int(len(X)),
        "n_train": int(len(Xtr)),
        "n_test": int(len(Xte)),
        "mae_pchembl": float(mean_absolute_error(yte, pred)),
        "r2_pchembl": float(r2_score(yte, pred)),
        "feature_columns": list(X.columns),
    }

    bundle = {
        "version": "chembl-aux-1.0",
        "task": "pchembl_regression",
        "model": model,
        "feature_columns": list(X.columns),
        "metrics": metrics,
        "note": (
            "Chemogenomic auxiliary model. Improves pathway ligandability priors; "
            "does not predict exhaled ppb by itself."
        ),
    }
    model_path = out_dir / "chembl_aux_model.joblib"
    metrics_path = out_dir / "chembl_aux_metrics.json"
    joblib.dump(bundle, model_path)
    metrics_path.write_text(json.dumps(metrics, indent=2))

    if update_knowledge_priors:
        voc_props_path = CHEMBL_DIR / "chembl_voc_properties.csv"
        voc_props = pd.read_csv(voc_props_path) if voc_props_path.exists() else None
        priors = distill_chembl_priors(df, kb=kb, voc_props=voc_props)
        priors["aux_model_metrics"] = metrics
        for p in (
            CHEMBL_DIR / "chembl_pathway_priors.json",
            KNOWLEDGE_DIR / "chembl_pathway_priors.json",
        ):
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(priors, indent=2))
        clear_knowledge_cache()

    return {
        "model_path": str(model_path),
        "metrics_path": str(metrics_path),
        "metrics": metrics,
    }
