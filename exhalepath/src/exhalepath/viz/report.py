from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from ..model.predict import PredictionResult


def save_biomarker_report(report, out_dir: Path) -> dict[str, Path]:
    """Write top-N exhaled biomarker panel (CSV/JSON/plot) for any disease × location."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = save_prediction_report(report.result, out_dir)
    panel = report.to_dataframe()
    panel_path = out_dir / "top_voc_biomarkers.csv"
    panel.to_csv(panel_path, index=False)
    summary_path = out_dir / "biomarker_summary.txt"
    summary_path.write_text(report.summary() + "\n\nNotes:\n" + "\n".join(f"- {n}" for n in report.notes))
    meta = {
        "disease_query": report.disease_query,
        "disease_id": report.disease_id,
        "disease_name": report.disease_name,
        "location_query": report.location_query,
        "location": report.location,
        "n_vocs_modeled": report.n_vocs_modeled,
        "top_n": len(report.top_vocs),
        "model_version": report.model_version,
    }
    meta_path = out_dir / "biomarker_meta.json"
    meta_path.write_text(__import__("json").dumps(meta, indent=2))
    # Ranked |Δppb| plot (quantity-first biomarker view)
    import matplotlib.pyplot as plt

    top = panel.head(min(20, len(panel))).iloc[::-1]
    fig, ax = plt.subplots(figsize=(9, 7))
    colors = ["#E45756" if v >= 0 else "#4C78A8" for v in top["delta_ppb"]]
    ax.barh(top["name"], top["delta_ppb"], color=colors)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("Predicted Δ exhaled VOC (ppb vs healthy)")
    ax.set_title(
        f"ExhalePath biomarkers — {report.disease_name}\n"
        f"@ {report.location.get('name') or report.location_query}"
    )
    fig.tight_layout()
    delta_plot = out_dir / "voc_delta_ppb.png"
    fig.savefig(delta_plot, dpi=160)
    plt.close(fig)
    paths.update(
        {
            "top_panel": panel_path,
            "summary": summary_path,
            "meta": meta_path,
            "delta_plot": delta_plot,
        }
    )
    return paths


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
