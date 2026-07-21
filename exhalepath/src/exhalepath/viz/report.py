from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..model.predict import PredictionResult


def save_prediction_report(result: PredictionResult, out_dir: Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = result.to_dataframe()
    csv_path = out_dir / "voc_predictions.csv"
    json_path = out_dir / "prediction_bundle.json"
    png_path = out_dir / "voc_fold_changes.png"
    pathway_path = out_dir / "pathway_scores.csv"
    cell_path = out_dir / "cell_states.csv"
    physio_path = out_dir / "physiology_traces.csv"

    df.to_csv(csv_path, index=False)
    json_path.write_text(result.bundle.model_dump_json(indent=2))
    pd.DataFrame([p.model_dump() for p in result.bundle.pathway_scores]).to_csv(
        pathway_path, index=False
    )
    if result.bundle.cell_states:
        pd.DataFrame([c.model_dump() for c in result.bundle.cell_states]).to_csv(
            cell_path, index=False
        )
    phys_rows = [
        p.physiology.model_dump()
        for p in result.bundle.predictions
        if p.physiology is not None
    ]
    if phys_rows:
        pd.DataFrame(phys_rows).to_csv(physio_path, index=False)

    top = df.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = ["#E45756" if v >= 0 else "#4C78A8" for v in top["log2_fold_change"]]
    ax.barh(top["name"], top["log2_fold_change"], color=colors)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("Predicted log2 fold-change vs healthy breath")
    ax.set_title(f"ExhalePath — {result.bundle.disease_name}")
    fig.tight_layout()
    fig.savefig(png_path, dpi=160)
    plt.close(fig)

    out = {
        "csv": csv_path,
        "json": json_path,
        "plot": png_path,
        "pathways": pathway_path,
    }
    if cell_path.exists():
        out["cell_states"] = cell_path
    if physio_path.exists():
        out["physiology"] = physio_path
    return out
