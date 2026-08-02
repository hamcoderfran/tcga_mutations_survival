"""Comprehensive visual + JSON/Markdown report for the fused stack."""

from __future__ import annotations

import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .types import StackResult


def _slug(text: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", str(text or "query").lower()).strip("_")
    return (s or "query")[:60]


def default_stack_dir(disease: str, location: str | None = None) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return Path("runs") / f"stack_{_slug(disease)}_{_slug(location or 'auto')}_{stamp}"


def result_to_dict(result: StackResult) -> dict[str, Any]:
    return {
        "disease_id": result.disease_id,
        "disease_name": result.disease_name,
        "location": result.location,
        "query": {
            "disease": result.query.disease,
            "location": result.query.location,
            "comorbidities": result.query.comorbidities,
            "age": result.query.age,
            "sex": result.query.sex,
            "genes": result.query.genes,
            "description": result.query.description,
            "smoking": result.query.smoking,
        },
        "summary": result.summary,
        "fusion_weights": result.fusion_weights,
        "fused_vocs": [
            {
                "voc_id": v.voc_id,
                "name": v.name,
                "fused_log2fc": v.fused_log2fc,
                "fused_delta_ppb": v.fused_delta_ppb,
                "fused_confidence": v.fused_confidence,
                "rrf_score": v.rrf_score,
                "n_models_agreeing": v.n_models_agreeing,
                "model_votes": v.model_votes,
                "model_confidences": v.model_confidences,
                "evidence": v.evidence,
                "epistemic_std": v.epistemic_std,
                "ci_low_log2fc": v.ci_low_log2fc,
                "ci_high_log2fc": v.ci_high_log2fc,
            }
            for v in result.fused_vocs
        ],
        "fused_aspects": {
            kind: [
                {
                    "id": h.id,
                    "name": h.name,
                    "score": h.score,
                    "kind": h.kind,
                    "evidence": h.evidence[:6],
                }
                for h in hits
            ]
            for kind, hits in result.fused_aspects.items()
        },
        "models": [
            {
                "model_id": m.model_id,
                "family": m.family,
                "aspect": m.aspect,
                "weight": m.weight,
                "status": m.status,
                "n_voc_signals": len(m.voc_signals),
                "n_aspects": len(m.aspects),
                "metadata": m.metadata,
                "notes": m.notes,
            }
            for m in result.model_outputs
        ],
        "notes": result.notes,
    }


def save_stack_report(result: StackResult, out_dir: Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    payload = result_to_dict(result)
    jp = out_dir / "STACK_RESULT.json"
    jp.write_text(json.dumps(payload, indent=2, default=str))
    paths["json"] = jp

    # CSV of fused VOCs
    import csv

    cp = out_dir / "fused_vocs.csv"
    with cp.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "rank",
                "voc_id",
                "name",
                "fused_log2fc",
                "fused_delta_ppb",
                "fused_confidence",
                "rrf_score",
                "n_models_agreeing",
                "epistemic_std",
                "ci_low_log2fc",
                "ci_high_log2fc",
            ]
        )
        for i, v in enumerate(result.fused_vocs, 1):
            w.writerow(
                [
                    i,
                    v.voc_id,
                    v.name,
                    f"{v.fused_log2fc:.4f}",
                    f"{v.fused_delta_ppb:.4f}",
                    f"{v.fused_confidence:.4f}",
                    f"{v.rrf_score:.4f}",
                    v.n_models_agreeing,
                    f"{v.epistemic_std:.4f}",
                    f"{v.ci_low_log2fc:.4f}",
                    f"{v.ci_high_log2fc:.4f}",
                ]
            )
    paths["csv"] = cp

    # model status table
    mp = out_dir / "model_status.csv"
    with mp.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model_id", "family", "aspect", "weight", "status", "n_voc", "n_aspects"])
        for m in result.model_outputs:
            w.writerow(
                [
                    m.model_id,
                    m.family,
                    m.aspect,
                    m.weight,
                    m.status,
                    len(m.voc_signals),
                    len(m.aspects),
                ]
            )
    paths["model_csv"] = mp

    # ---- dashboard figure ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Great Disease Stack · {result.disease_name}\n"
        f"{result.summary.get('n_models_ok')}/{result.summary.get('n_models_total')} models fused",
        fontsize=13,
        fontweight="bold",
    )
    top = list(reversed(result.fused_vocs[:12]))
    ax = axes[0, 0]
    if top:
        colors = ["#C44E52" if v.fused_log2fc >= 0 else "#4C78A8" for v in top]
        ax.barh([v.name for v in top], [v.fused_log2fc for v in top], color=colors)
        ax.axvline(0, color="#222", lw=0.8)
        ax.set_xlabel("Fused log2 fold-change")
        ax.set_title("Consensus VOC biomarkers")
    ax = axes[0, 1]
    if top:
        ax.barh([v.name for v in top], [v.n_models_agreeing for v in top], color="#54A24B")
        ax.set_xlabel("# models with |signal|≥0.05")
        ax.set_title("Cross-model agreement")
    ax = axes[1, 0]
    statuses = {"ok": 0, "degraded": 0, "skipped": 0}
    for m in result.model_outputs:
        statuses[m.status] = statuses.get(m.status, 0) + 1
    ax.bar(list(statuses.keys()), list(statuses.values()), color=["#54A24B", "#EECA3B", "#B8B8B8"])
    ax.set_title("Model health")
    ax.set_ylabel("count")
    ax = axes[1, 1]
    pw = (result.fused_aspects.get("pathway") or [])[:10]
    if pw:
        ax.barh([p.name[:28] for p in reversed(pw)], [p.score for p in reversed(pw)], color="#F58518")
        ax.set_title("Fused pathway aspects")
        ax.set_xlabel("score")
    else:
        ax.text(0.5, 0.5, "No pathway aspects", ha="center", va="center")
        ax.set_axis_off()
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    dash = out_dir / "stack_dashboard.png"
    fig.savefig(dash, dpi=170)
    plt.close(fig)
    paths["dashboard"] = dash

    # heatmap of model × top VOC votes
    top5 = result.fused_vocs[:10]
    model_ids = [m.model_id for m in result.model_outputs if m.status != "skipped"]
    if top5 and model_ids:
        mat = np.zeros((len(top5), len(model_ids)))
        for i, v in enumerate(top5):
            for j, mid in enumerate(model_ids):
                mat[i, j] = float(v.model_votes.get(mid, 0.0))
        fig, ax = plt.subplots(figsize=(max(8, len(model_ids) * 0.55), 5))
        im = ax.imshow(mat, aspect="auto", cmap="coolwarm")
        ax.set_yticks(range(len(top5)))
        ax.set_yticklabels([v.name for v in top5])
        ax.set_xticks(range(len(model_ids)))
        ax.set_xticklabels(model_ids, rotation=55, ha="right", fontsize=8)
        ax.set_title("Per-model log2fc votes (top fused VOCs)")
        fig.colorbar(im, ax=ax, fraction=0.03)
        fig.tight_layout()
        hm = out_dir / "model_vote_heatmap.png"
        fig.savefig(hm, dpi=160)
        plt.close(fig)
        paths["heatmap"] = hm

    # Markdown
    md = [
        f"# Great Disease Stack — {result.disease_name}",
        "",
        f"- Disease ID: `{result.disease_id}`",
        f"- Location: {result.location.get('name')}",
        f"- Models fused: **{result.summary.get('n_models_ok')}/{result.summary.get('n_models_total')}**",
        f"- Mean model agreement: {result.summary.get('mean_model_agreement')}",
        f"- Fusion: `{result.summary.get('fusion_method')}`",
        "",
        "## Fused VOC biomarkers",
        "",
        "| Rank | VOC | log2fc | Δppb | Conf | Agree |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for i, v in enumerate(result.fused_vocs, 1):
        md.append(
            f"| {i} | {v.name} | {v.fused_log2fc:+.3f} | {v.fused_delta_ppb:+.2f} | "
            f"{v.fused_confidence:.2f} | {v.n_models_agreeing} |"
        )
    md += ["", "## Models", ""]
    for m in result.model_outputs:
        md.append(
            f"- `{m.model_id}` ({m.family}/{m.aspect}) weight={m.weight:.2f} "
            f"status={m.status} vocs={len(m.voc_signals)}"
        )
    for kind, hits in result.fused_aspects.items():
        md += ["", f"## Aspect: {kind}", ""]
        for h in hits[:10]:
            md.append(f"- **{h.name}** (`{h.id}`) score={h.score:.3f}")
    md += ["", "![Dashboard](stack_dashboard.png)", ""]
    if "heatmap" in paths:
        md += ["![Votes](model_vote_heatmap.png)", ""]
    md += ["", "> Research / hypothesis-generation only — not a medical device.", ""]
    md_path = out_dir / "STACK_REPORT.md"
    md_path.write_text("\n".join(md))
    paths["md"] = md_path

    # HTML
    rows = []
    for i, v in enumerate(result.fused_vocs[:25], 1):
        tone = "#C44E52" if v.fused_log2fc >= 0 else "#4C78A8"
        rows.append(
            "<tr>"
            f"<td>{i}</td><td>{html.escape(v.name)}</td>"
            f"<td style='color:{tone}'>{v.fused_log2fc:+.3f}</td>"
            f"<td>{v.fused_delta_ppb:+.2f}</td><td>{v.fused_confidence:.2f}</td>"
            f"<td>{v.n_models_agreeing}</td></tr>"
        )
    model_rows = "".join(
        f"<tr><td><code>{html.escape(m.model_id)}</code></td><td>{html.escape(m.family)}</td>"
        f"<td>{html.escape(m.aspect)}</td><td>{m.weight:.2f}</td><td>{m.status}</td>"
        f"<td>{len(m.voc_signals)}</td></tr>"
        for m in result.model_outputs
    )
    import base64

    dash_b64 = base64.b64encode(dash.read_bytes()).decode("ascii")
    html_doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<title>Great Stack — {html.escape(result.disease_name)}</title>
<style>
body {{ font-family: 'IBM Plex Sans', Segoe UI, sans-serif; margin:0; background:
 linear-gradient(160deg,#eef3f7 0%,#f7f1e8 55%,#e8efe6 100%); color:#1b1f1c; }}
main {{ max-width: 1120px; margin:0 auto; padding:2rem 1.2rem 3rem; }}
h1 {{ font-family: 'Libre Baskerville', Georgia, serif; margin:0 0 .4rem; }}
.card {{ background:rgba(255,255,255,.85); border:1px solid #d5dcd6; border-radius:10px;
 padding:1rem 1.1rem; margin:1rem 0; box-shadow:0 8px 22px rgba(20,30,20,.06); }}
table {{ width:100%; border-collapse:collapse; font-size:.92rem; }}
th,td {{ padding:.4rem .5rem; border-bottom:1px solid #e3e8e3; }}
th {{ text-align:left; color:#5c665e; font-size:.78rem; text-transform:uppercase; }}
img {{ width:100%; border-radius:8px; border:1px solid #d5dcd6; }}
.pill {{ display:inline-block; background:#e7f3ea; padding:.2rem .55rem; border-radius:999px; margin-right:.35rem; }}
</style></head><body><main>
<h1>{html.escape(result.disease_name)}</h1>
<div>
 <span class="pill">{result.summary.get('n_models_ok')}/{result.summary.get('n_models_total')} models</span>
 <span class="pill">agree={result.summary.get('mean_model_agreement')}</span>
 <span class="pill">{html.escape(str(result.location.get('name') or ''))}</span>
</div>
<div class="card"><img alt="dashboard" src="data:image/png;base64,{dash_b64}"/></div>
<div class="card"><h2>Fused VOC biomarkers</h2>
<table><thead><tr><th>#</th><th>VOC</th><th>log2fc</th><th>Δppb</th><th>Conf</th><th>Agree</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
<div class="card"><h2>Model stack</h2>
<table><thead><tr><th>Model</th><th>Family</th><th>Aspect</th><th>Weight</th><th>Status</th><th>VOCs</th></tr></thead>
<tbody>{model_rows}</tbody></table></div>
<p style="color:#667066;font-size:.88rem">Research / hypothesis-generation only — not a medical device.</p>
</main></body></html>"""
    hp = out_dir / "STACK_REPORT.html"
    hp.write_text(html_doc)
    paths["html"] = hp

    (out_dir / "manifest.json").write_text(
        json.dumps({"artifacts": {k: str(v) for k, v in paths.items()}, "summary": result.summary}, indent=2)
    )
    return paths


def print_stack_console(result: StackResult, *, top_display: int = 15) -> None:
    from rich import box
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table

    console = Console()
    s = result.summary
    console.print(
        Panel(
            f"[bold]{result.disease_name}[/bold]\n"
            f"[dim]{result.disease_id}[/dim] · loc={result.location.get('name')}\n"
            f"models [green]{s.get('n_models_ok')}[/green]/{s.get('n_models_total')} "
            f"· agree={s.get('mean_model_agreement')}",
            title="Great Disease Stack",
            border_style="dark_sea_green4",
        )
    )
    table = Table(title="Fused VOC consensus", box=box.SIMPLE_HEAVY)
    for col in ("#", "VOC", "log2fc", "Δppb", "Conf", "Agree"):
        table.add_column(col, justify="right" if col != "VOC" else "left")
    for i, v in enumerate(result.fused_vocs[:top_display], 1):
        style = "red" if v.fused_log2fc >= 0 else "blue"
        table.add_row(
            str(i),
            v.name,
            f"[{style}]{v.fused_log2fc:+.3f}[/{style}]",
            f"{v.fused_delta_ppb:+.2f}",
            f"{v.fused_confidence:.2f}",
            str(v.n_models_agreeing),
        )
    console.print(table)

    mt = Table(title="Model stack", box=box.SIMPLE)
    mt.add_column("Model")
    mt.add_column("Aspect")
    mt.add_column("W", justify="right")
    mt.add_column("Status")
    mt.add_column("VOCs", justify="right")
    for m in result.model_outputs:
        st = {"ok": "green", "degraded": "yellow", "skipped": "dim"}.get(m.status, "white")
        mt.add_row(m.model_id, m.aspect, f"{m.weight:.2f}", f"[{st}]{m.status}[/{st}]", str(len(m.voc_signals)))
    console.print(mt)

    for kind in ("gene", "pathway", "microbe", "cell_state"):
        hits = result.fused_aspects.get(kind) or []
        if not hits:
            continue
        console.print(f"[bold cyan]{kind}[/bold cyan]: " + ", ".join(f"{h.name}({h.score:.2f})" for h in hits[:8]))
