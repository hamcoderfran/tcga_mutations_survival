from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import typer
from rich import print as rprint
from rich.table import Table

from .config import DEFAULT_GDC_PROJECTS, MODELS_DIR, PROCESSED_DIR
from .ingest.build_corpus import build_training_corpus
from .ingest.census_harvest import CENSUS_DIR, calibrate_atlas_from_census, harvest_census_compositions
from .model.predict import ExhalePathPredictor
from .model.train import train_calibrator
from .schemas import DiseaseQuery, TumorContext
from .viz.report import save_prediction_report

app = typer.Typer(
    name="exhalepath",
    help="ExhalePath: pathway-informed exhaled VOC (ppb) prediction for any disease.",
    add_completion=False,
)


@app.command("build-corpus")
def build_corpus_cmd(
    max_projects: Optional[int] = typer.Option(
        8,
        help="Limit GDC projects (default 8 for a fast rich corpus; None for all).",
    ),
    no_reactome: bool = typer.Option(False, help="Skip Reactome gene-set expansion"),
    offline_demo: bool = typer.Option(
        False,
        help="Build a large synthetic multi-project corpus without network calls",
    ),
    demo_cases: int = typer.Option(400, help="Cases per project in offline-demo mode"),
    out_dir: Path = typer.Option(PROCESSED_DIR, help="Output directory"),
):
    """Pull multi-cohort TCGA/GDC + pathway data and materialize a training corpus."""
    paths = build_training_corpus(
        projects=DEFAULT_GDC_PROJECTS[:max_projects] if max_projects else None,
        out_dir=out_dir,
        expand_reactome=not no_reactome,
        max_projects=max_projects,
        offline_demo=offline_demo,
        demo_cases_per_project=demo_cases,
    )
    rprint("[green]Corpus ready[/green]")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")


@app.command("train")
def train_cmd(
    case_features: Path = typer.Option(PROCESSED_DIR / "case_pathway_features.csv"),
    voc_targets: Path = typer.Option(PROCESSED_DIR / "voc_training_targets.csv"),
    out_dir: Path = typer.Option(MODELS_DIR),
):
    """Train per-VOC ppb/log2fc calibrators on the built corpus."""
    result = train_calibrator(
        case_features_path=case_features,
        voc_targets_path=voc_targets,
        out_dir=out_dir,
    )
    rprint("[green]Training complete[/green]", result["model_path"])


@app.command("predict")
def predict_cmd(
    disease: str = typer.Argument(..., help="Disease name, alias, or TCGA project"),
    stage: Optional[str] = typer.Option(None, help="Tumor stage, e.g. II, IIIA, IV"),
    site: Optional[str] = typer.Option(None, help="Primary site / position"),
    histology: Optional[str] = typer.Option(None, help="Histologic type"),
    tumor_type: Optional[str] = typer.Option(None, help="solid / hematologic / metastatic"),
    metastatic: bool = typer.Option(False, help="Mark as metastatic disease"),
    genes: Optional[str] = typer.Option(
        None, help="Comma-separated mutated genes, e.g. KRAS,TP53,CDKN2A"
    ),
    age: Optional[float] = typer.Option(None),
    sex: Optional[str] = typer.Option(None),
    smoking: Optional[str] = typer.Option(None, help="never|former|current"),
    out_dir: Optional[Path] = typer.Option(None, help="Write CSV/JSON/plot report here"),
    no_opentargets: bool = typer.Option(False, help="Disable Open Targets enrichment"),
    top: int = typer.Option(12, help="Rows to display"),
    mode: str = typer.Option(
        "hybrid",
        help="physiology | hybrid | legacy — cell/blood/alveolar physio vs prior-only",
    ),
    va: Optional[float] = typer.Option(None, help="Alveolar ventilation VA (L/min)"),
    q: Optional[float] = typer.Option(None, help="Cardiac output Q (L/min)"),
    cell_fractions: Optional[str] = typer.Option(
        None,
        help="scRNA cell-state fractions as state=frac pairs, e.g. tumor_epithelial_warburg=0.4,hepatocyte_ketogenic=0.5",
    ),
):
    """Predict exhaled VOC concentration shifts (ppb) for a disease context."""
    tumor = None
    if any([stage, site, histology, tumor_type, metastatic]):
        tumor = TumorContext(
            stage=stage,
            primary_site=site,
            histology=histology,
            tumor_type=tumor_type,
            metastatic=metastatic,
        )
    frac_map = {}
    if cell_fractions:
        for part in cell_fractions.split(","):
            if "=" not in part:
                continue
            k, v = part.split("=", 1)
            frac_map[k.strip()] = float(v)
    mode_norm = mode if mode in {"physiology", "hybrid", "legacy"} else "hybrid"
    query = DiseaseQuery(
        disease=disease,
        tumor=tumor,
        mutated_genes=[g.strip().upper() for g in genes.split(",")] if genes else [],
        age_years=age,
        sex=sex if sex in {"female", "male", "other"} else None,
        smoking_status=smoking if smoking in {"never", "former", "current"} else None,
        mode=mode_norm,  # type: ignore[arg-type]
        alveolar_ventilation_l_per_min=va,
        cardiac_output_l_per_min=q,
        cell_state_fractions=frac_map,
    )
    predictor = ExhalePathPredictor(use_opentargets=not no_opentargets)
    result = predictor.predict(query)

    table = Table(title=f"ExhalePath · {result.bundle.disease_name} · {mode_norm}")
    table.add_column("VOC")
    table.add_column("Healthy ppb", justify="right")
    table.add_column("Predicted ppb", justify="right")
    table.add_column("Δ ppb", justify="right")
    table.add_column("Fold", justify="right")
    table.add_column("Blood→Alv", justify="right")
    table.add_column("Drivers")
    for p in result.bundle.predictions[:top]:
        alv = f"{p.physiology.alveolar_fraction:.3f}" if p.physiology else "—"
        table.add_row(
            p.name,
            f"{p.healthy_ppb:.2f}",
            f"{p.predicted_ppb:.2f}",
            f"{p.delta_ppb:+.2f}",
            f"{p.fold_change:.2f}x",
            alv,
            ", ".join(p.top_pathway_drivers[:3]) or "—",
        )
    rprint(table)
    if result.bundle.cell_states:
        cs = Table(title="Affected cell states (density × activity)")
        cs.add_column("State")
        cs.add_column("Tissue")
        cs.add_column("Density", justify="right")
        cs.add_column("Activity", justify="right")
        cs.add_column("Source", justify="right")
        for s in result.bundle.cell_states[:8]:
            cs.add_row(
                s.name,
                s.tissue,
                f"{s.density:.3f}",
                f"{s.activity:.3f}",
                f"{s.effective_source:.3f}",
            )
        rprint(cs)
    for note in result.bundle.notes:
        rprint(f"[dim]• {note}[/dim]")

    if out_dir:
        paths = save_prediction_report(result, out_dir)
        rprint("[green]Wrote report:[/green]", json.dumps({k: str(v) for k, v in paths.items()}, indent=2))


@app.command("list-diseases")
def list_diseases():
    """List curated disease atlas entries."""
    from .knowledge.loader import KnowledgeBase

    kb = KnowledgeBase()
    for d in kb.diseases.values():
        aliases = ", ".join(d.get("aliases", [])[:4])
        rprint(f"[bold]{d['disease_id']}[/bold] — {d['name']}  ({aliases})")


@app.command("list-vocs")
def list_vocs():
    """List VOC catalog with healthy breath baselines."""
    from .knowledge.loader import KnowledgeBase

    kb = KnowledgeBase()
    for v in kb.vocs.values():
        rprint(
            f"{v['voc_id']:18} {v['name']:22} healthy≈{v['healthy_ppb_median']} ppb  CAS {v.get('cas')}"
        )


@app.command("harvest-census")
def harvest_census_cmd(
    out_dir: Path = typer.Option(CENSUS_DIR, help="Output directory for Census summaries"),
    min_cells: int = typer.Option(5_000, help="Min Census cells for a disease to harvest"),
    top_n: int = typer.Option(100, help="Max diseases to harvest (re-ranked by cell abundance)"),
    calibrate: bool = typer.Option(True, help="Calibrate cell_state_atlas from harvested fractions"),
    no_opportunistic: bool = typer.Option(
        False, help="Only harvest top-100 US list matches (skip other data-rich Census diseases)"
    ),
):
    """
    Download CELLxGENE Census single-cell compositions for top US diseases +
    healthy cells across body tissues. Lean into diseases/tissues with the most cells.
    Stores summaries (not full raw matrices), then calibrates ExhalePath cell states.
    """
    paths = harvest_census_compositions(
        out_dir=out_dir,
        min_disease_cells=min_cells,
        top_n_diseases=top_n,
        include_opportunistic=not no_opportunistic,
    )
    rprint("[green]Census harvest complete[/green]")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")
    if calibrate:
        atlas = calibrate_atlas_from_census(census_dir=out_dir)
        rprint(f"[green]Atlas calibrated:[/green] {atlas}")


@app.command("audit")
def audit_cmd(
    out: Path = typer.Option(
        Path("runs/audit/audit_report.json"),
        help="Where to write the full JSON audit report",
    ),
):
    """Run literature accuracy, calibrator holdout, zero-shot, and stress audits."""
    from .eval.audit import run_full_audit

    report = run_full_audit(out_path=out)
    lit = report.literature_accuracy
    rprint(
        f"[bold]Audit {'PASSED' if report.passed else 'FAILED'}[/bold]\n"
        f"  literature case pass rate: {lit.get('case_pass_rate', 0):.1%}\n"
        f"  directional accuracy:      {lit.get('directional_accuracy')}\n"
        f"  min-fold accuracy:         {lit.get('min_fold_accuracy')}\n"
        f"  zero-shot pass rate:       {report.zero_shot.get('pass_rate')}\n"
        f"  stress failures:           {len(report.stress.get('failures', []))}\n"
        f"  calibrator:                {report.calibrator_holdout}\n"
        f"  report: {out}"
    )
    if not report.passed:
        for c in lit.get("cases", []):
            if not c.get("passed"):
                rprint(f"  [red]FAIL[/red] {c['case_id']}: {c.get('errors')}")
        for f in report.stress.get("failures", []):
            rprint(f"  [red]STRESS[/red] {f}")
        raise typer.Exit(code=1)


@app.command("eval-multisite")
def eval_multisite_cmd(
    out_dir: Path = typer.Option(Path("runs/multisite_eval"), help="Output directory"),
    mode: str = typer.Option("hybrid", help="physiology | hybrid | legacy"),
):
    """Test 20 diseases × 5 sites against literature VOC expectations; report accuracy."""
    from .eval.multisite import run_multisite_literature_eval

    report = run_multisite_literature_eval(out_dir=out_dir, mode=mode)
    o = report["overall"]
    rprint(
        f"[bold]Multisite literature eval[/bold]\n"
        f"  diseases × sites:     {o['n_diseases']} × {o['n_sites_per_disease']}\n"
        f"  predictions:          {o['n_predictions']}\n"
        f"  directional accuracy: {o['directional_accuracy']:.1%}\n"
        f"  min-fold accuracy:    {o['min_fold_accuracy']}\n"
        f"  site sensitivity:     {o['site_sensitivity']}\n"
        f"  disease pass rate:    {o['disease_pass_rate']:.1%} "
        f"({o['n_diseases_passed']}/{o['n_diseases']})\n"
        f"  composite accuracy:   {o['composite_accuracy']:.1%}\n"
        f"  grade: {o['grading']}\n"
        f"  report: {out_dir / 'multisite_accuracy_report.json'}"
    )
    fails = [d for d in report["diseases"] if not d["passed"]]
    if fails:
        rprint("[yellow]Diseases below pass threshold:[/yellow]")
        for d in fails:
            rprint(
                f"  • {d['disease']}: dir={d['directional_accuracy']:.0%} "
                f"fold={d['min_fold_accuracy']} site={d['site_sensitivity']}"
            )


def main():
    app()


if __name__ == "__main__":
    main()
