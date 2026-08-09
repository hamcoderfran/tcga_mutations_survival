"""Rich one-line biomarker dashboard: multi-panel plots + HTML report."""

from __future__ import annotations

import base64
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

from .report import save_biomarker_report


def _slug(text: str) -> str:
    s = re.sub(r"[^a-z0-9]+", "_", str(text or "query").lower()).strip("_")
    return (s or "query")[:60]


def default_report_dir(disease: str, location: str | None = None) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    loc = _slug(location or "auto")
    return Path("runs") / f"voc_{_slug(disease)}_{loc}_{stamp}"


def _b64_png(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


def save_visual_dashboard(report, out_dir: Path) -> dict[str, Path]:
    """
    Write a full visual pack:
      - standard CSV/JSON/plots (save_biomarker_report)
      - dashboard.png (4-panel overview)
      - REPORT.html (comprehensive one-page view)
      - REPORT.md
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = save_biomarker_report(report, out_dir)

    panel = report.to_dataframe()
    meta = report.result.bundle.metadata or {}
    pathways = list(report.result.bundle.pathway_scores or [])
    cells = list(report.result.bundle.cell_states or [])

    # ---- multi-panel dashboard ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"ExhalePath · {report.disease_name}\n"
        f"@ {report.location.get('name') or report.location_query} · {report.model_version}",
        fontsize=13,
        fontweight="bold",
    )

    # 1) Δppb top 15
    ax = axes[0, 0]
    top = panel.head(min(15, len(panel))).iloc[::-1]
    colors = ["#C44E52" if v >= 0 else "#4C78A8" for v in top["delta_ppb"]]
    ax.barh(top["name"], top["delta_ppb"], color=colors)
    ax.axvline(0, color="#222", lw=0.8)
    ax.set_xlabel("Δ exhaled VOC (ppb vs healthy)")
    ax.set_title("Top quantity shifts")

    # 2) log2fc top 15
    ax = axes[0, 1]
    colors = ["#C44E52" if v >= 0 else "#4C78A8" for v in top["log2_fold_change"]]
    ax.barh(top["name"], top["log2_fold_change"], color=colors)
    ax.axvline(0, color="#222", lw=0.8)
    ax.set_xlabel("log2 fold-change vs healthy")
    ax.set_title("Directional fold changes")

    # 3) pathway scores
    ax = axes[1, 0]
    pw = sorted(pathways, key=lambda p: abs(float(p.score)), reverse=True)[:12]
    if pw:
        names = [p.name[:28] for p in pw][::-1]
        vals = [float(p.score) for p in pw][::-1]
        ax.barh(names, vals, color="#54A24B")
        ax.set_xlabel("Pathway activity score")
        ax.set_title("Pathway drivers")
    else:
        ax.text(0.5, 0.5, "No pathway scores", ha="center", va="center")
        ax.set_axis_off()

    # 4) cell-state effective sources
    ax = axes[1, 1]
    cs = sorted(cells, key=lambda c: float(c.effective_source), reverse=True)[:10]
    if cs:
        names = [c.name[:28] for c in cs][::-1]
        vals = [float(c.effective_source) for c in cs][::-1]
        ax.barh(names, vals, color="#EECA3B")
        ax.set_xlabel("Density × activity")
        ax.set_title("Affected cell states")
    else:
        ax.text(0.5, 0.5, "No cell-state activity", ha="center", va="center")
        ax.set_axis_off()

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    dash_path = out_dir / "dashboard.png"
    fig.savefig(dash_path, dpi=170)
    plt.close(fig)
    paths["dashboard"] = dash_path

    # ---- markdown report ----
    md_lines = [
        f"# ExhalePath report — {report.disease_name}",
        "",
        f"- **Disease ID:** `{report.disease_id}`",
        f"- **Location:** {report.location.get('name') or report.location_query}",
        f"- **Model:** {report.model_version}",
        f"- **VOCs modeled:** {report.n_vocs_modeled} · showing top {len(report.top_vocs)}",
    ]
    if meta.get("zero_shot"):
        md_lines.append(
            f"- **Zero-shot:** mode=`{meta.get('zero_shot_mode')}` "
            f"mechanism_confidence=`{meta.get('mechanism_confidence')}`"
        )
    md_lines.extend(["", "## Top VOC biomarkers", "", "| Rank | VOC | Healthy | Predicted | Δppb | Fold | Conf |", "|---|---|---:|---:|---:|---:|---:|"])
    for i, p in enumerate(report.top_vocs[:20], 1):
        md_lines.append(
            f"| {i} | {p.name} | {p.healthy_ppb:.2f} | {p.predicted_ppb:.2f} | "
            f"{p.delta_ppb:+.2f} | {p.fold_change:.2f}x | {p.confidence:.2f} |"
        )
    md_lines.extend(["", "## Pathway drivers", ""])
    for p in pw[:10]:
        md_lines.append(f"- **{p.name}** (`{p.pathway_id}`) score={p.score:.3f}")
    if report.mechanisms:
        md_lines.extend(["", "## Why (mechanisms)", ""])
        for m in report.mechanisms[:8]:
            md_lines.append(f"- {m.get('why')}")
    md_lines.extend(["", "## Notes", ""])
    for n in report.notes:
        md_lines.append(f"- {n}")
    md_lines.extend(
        [
            "",
            "## Figures",
            "",
            "![Dashboard](dashboard.png)",
            "",
            "![Δppb](voc_delta_ppb.png)",
            "",
            "![log2fc](voc_fold_changes.png)",
            "",
        ]
    )
    md_path = out_dir / "REPORT.md"
    md_path.write_text("\n".join(md_lines) + "\n")
    paths["report_md"] = md_path

    # ---- HTML report ----
    rows_html = []
    for i, p in enumerate(report.top_vocs[:25], 1):
        tone = "#C44E52" if p.delta_ppb >= 0 else "#4C78A8"
        rows_html.append(
            "<tr>"
            f"<td>{i}</td><td>{html.escape(p.name)}</td>"
            f"<td style='text-align:right'>{p.healthy_ppb:.2f}</td>"
            f"<td style='text-align:right'>{p.predicted_ppb:.2f}</td>"
            f"<td style='text-align:right;color:{tone}'>{p.delta_ppb:+.2f}</td>"
            f"<td style='text-align:right'>{p.fold_change:.2f}x</td>"
            f"<td style='text-align:right'>{p.confidence:.2f}</td>"
            "</tr>"
        )
    pw_html = "".join(
        f"<li><b>{html.escape(p.name)}</b> "
        f"<code>{html.escape(p.pathway_id)}</code> — {p.score:.3f}</li>"
        for p in pw[:12]
    )
    mech_html = "".join(
        f"<li>{html.escape(str(m.get('why') or ''))}</li>"
        for m in (report.mechanisms or [])[:8]
    )
    notes_html = "".join(f"<li>{html.escape(n)}</li>" for n in report.notes)
    zs_banner = ""
    if meta.get("zero_shot"):
        zs_banner = (
            f"<div class='zs'>Zero-shot · mode=<b>{html.escape(str(meta.get('zero_shot_mode')))}</b> · "
            f"mechanism confidence=<b>{meta.get('mechanism_confidence')}</b></div>"
        )
    dash_b64 = _b64_png(dash_path)
    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>ExhalePath — {html.escape(report.disease_name)}</title>
<style>
body {{ font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif; margin: 0; background:
  linear-gradient(165deg, #f7f3ea 0%, #e8eef2 55%, #f3ebe3 100%); color: #1c1a17; }}
main {{ max-width: 1100px; margin: 0 auto; padding: 2rem 1.25rem 3rem; }}
h1 {{ font-family: 'Libre Baskerville', Georgia, serif; font-size: 2rem; margin: 0 0 .35rem; }}
.sub {{ color: #5a5348; margin-bottom: 1.25rem; }}
.zs {{ background: #fff4d6; border-left: 4px solid #d4a017; padding: .75rem 1rem; margin: 1rem 0; }}
.card {{ background: rgba(255,255,255,.82); border: 1px solid #d9d0c3; border-radius: 10px;
  padding: 1rem 1.1rem; margin: 1rem 0; box-shadow: 0 8px 24px rgba(40,30,10,.06); }}
table {{ width: 100%; border-collapse: collapse; font-size: .92rem; }}
th, td {{ padding: .45rem .5rem; border-bottom: 1px solid #e4ddd2; }}
th {{ text-align: left; font-size: .78rem; text-transform: uppercase; letter-spacing: .04em; color: #6a6258; }}
img.dash {{ width: 100%; border-radius: 8px; border: 1px solid #d9d0c3; }}
.disclaimer {{ font-size: .85rem; color: #6a6258; margin-top: 1.5rem; }}
code {{ background: #f0ebe3; padding: .1rem .35rem; border-radius: 4px; }}
</style>
</head>
<body>
<main>
  <h1>{html.escape(report.disease_name)}</h1>
  <div class="sub">
    Location: <b>{html.escape(str(report.location.get('name') or report.location_query))}</b>
    · Disease ID: <code>{html.escape(report.disease_id)}</code>
    · Model: {html.escape(report.model_version)}
    · {report.n_vocs_modeled} VOCs modeled
  </div>
  {zs_banner}
  <div class="card">
    <img class="dash" alt="dashboard" src="data:image/png;base64,{dash_b64}"/>
  </div>
  <div class="card">
    <h2>Top VOC biomarkers</h2>
    <table>
      <thead><tr><th>#</th><th>VOC</th><th>Healthy</th><th>Predicted</th><th>Δppb</th><th>Fold</th><th>Conf</th></tr></thead>
      <tbody>
      {''.join(rows_html)}
      </tbody>
    </table>
  </div>
  <div class="card">
    <h2>Pathway drivers</h2>
    <ul>{pw_html or '<li>None</li>'}</ul>
  </div>
  <div class="card">
    <h2>Why (mechanisms)</h2>
    <ul>{mech_html or '<li>No mechanism narrative</li>'}</ul>
  </div>
  <div class="card">
    <h2>Notes</h2>
    <ul>{notes_html}</ul>
  </div>
  <p class="disclaimer">Research / hypothesis-generation only — not a medical device or clinical diagnosis.
  Validate against breath GC-MS / PTR-MS cohorts.</p>
</main>
</body>
</html>
"""
    html_path = out_dir / "REPORT.html"
    html_path.write_text(html_doc)
    paths["report_html"] = html_path

    # index pointer for latest run
    latest = Path("runs") / "latest"
    try:
        if latest.is_symlink() or latest.exists():
            if latest.is_symlink() or latest.is_file():
                latest.unlink()
            # if directory, leave alone and write pointer file instead
        if not latest.exists():
            latest.symlink_to(out_dir.resolve())
    except OSError:
        (Path("runs") / "LATEST.txt").write_text(str(out_dir.resolve()) + "\n")

    # Research-facing exports (literature overlay, Methods, Prism/GraphPad long CSV)
    from .analysis_export import (
        literature_overlay,
        methods_markdown,
        next_experiments,
        write_graphpad_long_csv,
        write_literature_overlay_csv,
    )

    pred_map = {
        str(row.get("voc_id") or row.get("name")): float(row["log2_fold_change"])
        for _, row in panel.iterrows()
        if "log2_fold_change" in row and row.get("voc_id") is not None
    }
    # panel may use index; also pull from top_vocs
    for p in report.top_vocs:
        pred_map[p.voc_id] = float(p.log2_fold_change)
    overlay = literature_overlay(report.disease_id, pred_map)
    if not overlay.get("available"):
        overlay = literature_overlay(report.disease_name, pred_map)
    paths["literature_csv"] = write_literature_overlay_csv(
        overlay, out_dir / "literature_overlay.csv"
    )
    (out_dir / "literature_overlay.json").write_text(json.dumps(overlay, indent=2))
    paths["literature_json"] = out_dir / "literature_overlay.json"

    gp_rows = []
    for p in report.top_vocs:
        gp_rows.append(
            {
                "group": report.disease_id,
                "voc": p.name,
                "metric": "log2_fold_change",
                "estimate": f"{p.log2_fold_change:.6f}",
                "ci_low": f"{getattr(p, 'ci_low_log2fc', getattr(p, 'ci_low_ppb', ''))}",
                "ci_high": f"{getattr(p, 'ci_high_log2fc', getattr(p, 'ci_high_ppb', ''))}",
                "source": "exhalepath_biomarker",
                "confidence": f"{p.confidence:.4f}",
            }
        )
    paths["graphpad_csv"] = write_graphpad_long_csv(
        gp_rows, out_dir / "graphpad_voc_long.csv"
    )
    tips = next_experiments(
        top_vocs=[p.voc_id for p in report.top_vocs[:8]],
        overlay=overlay,
    )
    (out_dir / "METHODS.md").write_text(
        methods_markdown(
            tool="voc / ExhaleBiomarkerEngine",
            disease=report.disease_name,
            location=str(report.location.get("name") or report.location_query or ""),
            model_version=report.model_version,
            extra_bullets=[
                f"Zero-shot: {bool(meta.get('zero_shot'))} mode={meta.get('zero_shot_mode')}.",
                f"VOCs modeled: {report.n_vocs_modeled}.",
            ],
        )
    )
    paths["methods"] = out_dir / "METHODS.md"
    (out_dir / "NEXT_EXPERIMENTS.md").write_text(
        "# Next experiments\n\n" + "\n".join(f"- {t}" for t in tips) + "\n"
    )
    paths["next_experiments"] = out_dir / "NEXT_EXPERIMENTS.md"

    # Append literature + next-experiment sections to REPORT.md if present
    md_path = out_dir / "REPORT.md"
    if md_path.exists():
        extra = ["", "## Literature overlay", ""]
        if overlay.get("available"):
            extra.append(
                f"Directional agreement: **{overlay.get('n_agree')}/{overlay.get('n_compared')}**"
            )
            extra += [
                "",
                "| VOC | Literature | Predicted | Agree |",
                "|---|---:|---:|:---:|",
            ]
            for r in overlay.get("rows") or []:
                pred = "—" if r["predicted_log2fc"] is None else f"{r['predicted_log2fc']:+.2f}"
                ag = "—" if r["agree"] is None else ("yes" if r["agree"] else "no")
                extra.append(
                    f"| {r['voc_id']} | {r['literature_log2fc']:+.2f} | {pred} | {ag} |"
                )
            if overlay.get("refs"):
                extra += ["", "### References", ""]
                for ref in overlay["refs"]:
                    extra.append(
                        f"- {ref.get('title')} ({ref.get('year')}) — `{ref.get('doi')}`"
                    )
        else:
            extra.append("_No curated literature panel matched this disease id._")
        extra += ["", "## Next experiments", ""] + [f"- {t}" for t in tips] + ["",
            "Also see `METHODS.md`, `graphpad_voc_long.csv`, `literature_overlay.csv`.", ""]
        md_path.write_text(md_path.read_text() + "\n".join(extra))

    manifest = {
        "disease": report.disease_name,
        "disease_id": report.disease_id,
        "location": report.location,
        "out_dir": str(out_dir),
        "artifacts": {k: str(v) for k, v in paths.items()},
        "zero_shot": bool(meta.get("zero_shot")),
        "zero_shot_mode": meta.get("zero_shot_mode"),
        "mechanism_confidence": meta.get("mechanism_confidence"),
        "literature_overlay_available": bool(overlay.get("available")),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    man_path = out_dir / "manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2))
    paths["manifest"] = man_path
    return paths


def print_comprehensive_console(report, *, top_display: int = 15) -> None:
    """Rich terminal view used by one-line voc invocations."""
    from rich import box
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table

    console = Console()
    meta = report.result.bundle.metadata or {}
    loc = report.location.get("name") or report.location_query
    header = (
        f"[bold]{report.disease_name}[/bold]\n"
        f"[dim]{report.disease_id}[/dim]  ·  location: [cyan]{loc}[/cyan]\n"
        f"model: {report.model_version}  ·  {report.n_vocs_modeled} VOCs"
    )
    if meta.get("zero_shot"):
        header += (
            f"\n[yellow]zero-shot[/yellow] mode={meta.get('zero_shot_mode')} "
            f"mech_conf={meta.get('mechanism_confidence')}"
        )
    console.print(Panel(header, title="ExhalePath", border_style="sea_green3"))

    table = Table(title="Top exhaled VOC biomarkers", box=box.SIMPLE_HEAVY)
    table.add_column("#", justify="right")
    table.add_column("VOC")
    table.add_column("Healthy", justify="right")
    table.add_column("Predicted", justify="right")
    table.add_column("Δ ppb", justify="right")
    table.add_column("Fold", justify="right")
    table.add_column("Conf", justify="right")
    for i, p in enumerate(report.top_vocs[:top_display], 1):
        style = "red" if p.delta_ppb >= 0 else "blue"
        table.add_row(
            str(i),
            p.name,
            f"{p.healthy_ppb:.2f}",
            f"{p.predicted_ppb:.2f}",
            f"[{style}]{p.delta_ppb:+.2f}[/{style}]",
            f"{p.fold_change:.2f}x",
            f"{p.confidence:.2f}",
        )
    console.print(table)

    pw = sorted(
        report.result.bundle.pathway_scores or [],
        key=lambda p: abs(float(p.score)),
        reverse=True,
    )[:8]
    if pw:
        pt = Table(title="Pathway drivers", box=box.SIMPLE)
        pt.add_column("Pathway")
        pt.add_column("Score", justify="right")
        for p in pw:
            pt.add_row(p.name, f"{p.score:.3f}")
        console.print(pt)

    if report.mechanisms:
        console.print("[bold cyan]Why[/bold cyan]")
        for m in report.mechanisms[:5]:
            console.print(f"  • {m.get('why')}")

    for note in report.notes[:6]:
        console.print(f"[dim]• {note}[/dim]")
