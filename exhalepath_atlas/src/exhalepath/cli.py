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


@app.command("eval-completion")
def eval_completion_cmd(
    out_dir: Path = typer.Option(Path("runs/completion"), help="Output directory"),
):
    """
    Full dataset completion gate: pytest + audit + multisite + public breath.

    Exits non-zero if any research-ready gate fails.
    """
    from .eval.completion import run_completion_suite

    report = run_completion_suite(out_dir=out_dir)
    g = report["gates"]
    rprint(f"[bold]Completion {'PASSED' if report['passed'] else 'FAILED'}[/bold]")
    for k, v in g.items():
        mark = "✓" if v else "✗"
        rprint(f"  {mark} {k}")
    pb = report["public_breath"]
    rprint(
        f"  lit recall@15={pb['literature'].get('mean_elevated_recall_at_k')} "
        f"dir={pb['literature'].get('mean_directional_accuracy')}"
    )
    rprint(
        f"  sci_data elev_dir={pb['scientific_data'].get('mean_elevated_directional_accuracy')} "
        f"recall={pb['scientific_data'].get('mean_elevated_recall_at_k')}"
    )
    rprint(f"  report: {out_dir / 'completion_report.json'}")
    if not report["passed"]:
        raise typer.Exit(code=1)


@app.command("integrate-datasources")
def integrate_datasources_cmd(
    offline: bool = typer.Option(
        False, help="Skip live HTTP; still write curated priority 1–14 tables"
    ),
    priorities: Optional[str] = typer.Option(
        None, help="Comma-separated priorities to run, e.g. 1,2,13,14 (default: all 1–14)"
    ),
):
    """
    Harvest + fuse priority datasources 1–14 into exhalepath_atlas knowledge/.

    1 metabolomics repos · 2 HMDB · 3 λ partition · 4 mVOC · 5 PubChem · 6 Reactome
    7 GTEx priors · 8 Open Targets/GWAS · 9 GDC/TCGA · 10 BindingDB · 11 blood proxy · 12 NIST RI
    13 HBDB/VOLATILOME · 14 KEGG VOC pathways
    """
    from pathlib import Path as P

    from .datasources import integrate_all_datasources

    root = P(__file__).resolve().parents[2]
    pri = None
    if priorities:
        pri = [int(x.strip()) for x in priorities.split(",") if x.strip()]
    manifest = integrate_all_datasources(root=root, offline=offline, priorities=pri)
    rprint(
        f"[green]Integrated {manifest['n_datasources']} datasources[/green] "
        f"→ {manifest.get('capability', {}).get('n_diseases')} diseases capable"
    )
    for d in manifest["datasources"]:
        rprint(f"  [{d['priority']:02d}] {d['key']}: {d['title']}")


@app.command("harvest-public-breath")
def harvest_public_breath_cmd(
    out_dir: Path = typer.Option(None, help="Download dir (default data/public_breath)"),
):
    """Download Scientific Data breathomics peak tables and build public validation cases."""
    from .config import DATA_DIR
    from .ingest.public_breath import build_public_breath_benchmark

    paths = build_public_breath_benchmark(out_dir=out_dir or (DATA_DIR / "public_breath"))
    rprint("[green]Public breath benchmark ready[/green]")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")


@app.command("eval-public-breath")
def eval_public_breath_cmd(
    out_dir: Path = typer.Option(Path("runs/public_breath_eval")),
    top_k: int = typer.Option(15, help="Top-k predicted VOCs for recall"),
    mode: str = typer.Option("hybrid"),
):
    """Evaluate ExhalePath against public breath VOC benchmarks (top-k + direction)."""
    from .eval.public_breath import evaluate_public_breath

    report = evaluate_public_breath(top_k=top_k, mode=mode, out_dir=out_dir)
    o = report["overall"]
    rprint(
        f"[bold]Public breath eval[/bold]\n"
        f"  cases: {o['n_cases']}\n"
        f"  mean elevated recall@{o['top_k']}: {o['mean_elevated_recall_at_k']}\n"
        f"  mean directional accuracy: {o['mean_directional_accuracy']}\n"
        f"  report: {out_dir / 'public_breath_eval.json'}"
    )


@app.command("harvest-clinical-comorbidity")
def harvest_clinical_comorbidity_cmd():
    """
    Harvest Magdeburg SZ breath + ST003181 depression plasma (n=401) and build
    comorbidity clinical benchmarks (depression±obesity, schizophrenia±heart).
    """
    from .ingest.clinical_comorbidity import build_comorbidity_clinical_benchmarks

    path = build_comorbidity_clinical_benchmarks()
    rprint(f"[green]Clinical comorbidity benchmarks:[/green] {path}")


@app.command("eval-comorbidity-clinical")
def eval_comorbidity_clinical_cmd(
    out_dir: Path = typer.Option(Path("runs/comorbidity_clinical_eval")),
    top_k: int = typer.Option(15, help="Top-k predicted VOCs for recall"),
    mode: str = typer.Option("hybrid"),
    comorbidity_weight: float = typer.Option(0.65),
):
    """Evaluate comorbidity-aware predictions against clinical breath/metabolome cases."""
    from .eval.comorbidity_clinical import evaluate_comorbidity_clinical

    report = evaluate_comorbidity_clinical(
        top_k=top_k,
        mode=mode,
        comorbidity_weight=comorbidity_weight,
        out_dir=out_dir,
    )
    rprint("[bold]Comorbidity clinical eval[/bold]")
    rprint(f"  cases: {report['n_cases']}")
    rprint(f"  mean elevated recall@{top_k}: {report['mean_elevated_recall_at_k']}")
    rprint(f"  mean elevated directional: {report['mean_elevated_directional_accuracy']}")
    rprint(
        f"  mean suppressed directional: {report['mean_suppressed_directional_accuracy']}"
    )
    rprint(f"  mean directional accuracy: {report['mean_directional_accuracy']}")
    for c in report["cases"]:
        rprint(
            f"  • {c['case_id']}: dir={c['directional_accuracy']} "
            f"elev@k={c['elevated_recall_at_k']} "
            f"comorbid={c['comorbidities']}"
        )
    rprint(f"  report: {out_dir / 'comorbidity_clinical_eval.json'}")


@app.command("build-mechanism-packs")
def build_mechanism_packs_cmd(
    out: Path = typer.Option(
        Path("data/knowledge/disease_mechanism_packs.json"),
        help="Output JSON for all atlas diseases",
    ),
    top_vocs: int = typer.Option(15, help="VOCs to explain per disease"),
    mode: str = typer.Option("hybrid"),
):
    """
    Build WHY-mechanism packs for ~100 atlas diseases:

    VOC Δppb → pathways → Census cell populations → driver genes / alterations.
    """
    from .explain.mechanisms import build_all_mechanism_packs

    path = build_all_mechanism_packs(out_path=out, top_n_vocs=top_vocs, mode=mode)
    rprint(f"[green]Mechanism packs written:[/green] {path}")


@app.command("explain")
def explain_cmd(
    disease: str = typer.Argument(...),
    location: Optional[str] = typer.Option(None, "--location", "-l"),
    genes: Optional[str] = typer.Option(None),
    top: int = typer.Option(10, help="Top VOCs to explain"),
    out: Optional[Path] = typer.Option(None, help="Write mechanism pack JSON"),
):
    """Explain why top VOCs change: pathways, cell states, Census populations, genes."""
    from .explain.mechanisms import MechanismExplainer

    explainer = MechanismExplainer(use_opentargets=False)
    pack = explainer.build_disease_pack(
        disease,
        location=location,
        genes=[g.strip().upper() for g in genes.split(",")] if genes else None,
        top_n=top,
    )
    rprint(f"[bold]{pack['disease_name']}[/bold] @ {pack['location'].get('name')}")
    if pack.get("census", {}).get("top_cell_types"):
        rprint("[cyan]Top Census cell populations[/cyan]")
        for ct in pack["census"]["top_cell_types"][:5]:
            rprint(
                f"  • {ct['cell_type']} ({ct['tissue']}) n={ct['n_cells']:,} "
                f"({ct['fraction']:.1%})"
            )
    for m in pack["voc_mechanisms"][:top]:
        rprint(f"\n[green]{m['name']}[/green] Δppb={m['delta_ppb']:+.2f} ({m['direction']})")
        rprint(f"  {m['why']}")
        if m.get("driver_genes"):
            rprint(f"  genes: {', '.join(m['driver_genes'][:8])}")
    if out:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(pack, indent=2))
        rprint(f"[green]Wrote[/green] {out}")


@app.command("harvest-chembl")
def harvest_chembl_cmd(
    out_dir: Path = typer.Option(None, help="Output dir (default data/chembl)"),
    max_per_target: int = typer.Option(5_000, help="Max activities per gene target"),
    max_rows: Optional[int] = typer.Option(
        500_000, help="Cap total harvested activity rows (None = no cap)"
    ),
    max_genes: Optional[int] = typer.Option(None, help="Limit seed genes (debug)"),
    offline_demo: bool = typer.Option(
        False,
        help="Synthesize multi-million ChEMBL-like rows without API (training drill)",
    ),
    demo_rows: int = typer.Option(1_000_000, help="Rows to synthesize in offline-demo"),
    apply_lambda_hints: bool = typer.Option(
        True, help="Merge VOC AlogP→λ hints into physio_constants (setdefault)"
    ),
):
    """
    Harvest ChEMBL bioactivities for VOC-pathway seed genes (+ VOC physchem).

    ChEMBL (~24M activities) trains the chemogenomic middle layer — enzyme/pathway
    modulation and VOC physicochemical priors — NOT exhaled ppb labels.
    """
    from .config import CHEMBL_DIR
    from .ingest.chembl_harvest import (
        apply_chembl_lambda_hints_to_physio,
        harvest_chembl_for_pathways,
    )

    paths = harvest_chembl_for_pathways(
        out_dir=out_dir or CHEMBL_DIR,
        max_per_target=max_per_target,
        max_rows=None if offline_demo else max_rows,
        max_genes=max_genes,
        offline_demo=offline_demo,
        demo_rows=demo_rows,
    )
    rprint("[green]ChEMBL harvest complete[/green]")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")
    if apply_lambda_hints and paths.get("priors"):
        phys = apply_chembl_lambda_hints_to_physio(priors_path=paths["priors"])
        rprint(f"[green]Physio λ hints updated:[/green] {phys}")


@app.command("train-chembl")
def train_chembl_cmd(
    activities: Path = typer.Option(
        None, help="chembl_pathway_activities.csv (default data/chembl/...)"
    ),
    out_dir: Path = typer.Option(MODELS_DIR),
    max_rows: Optional[int] = typer.Option(
        2_000_000, help="Subsample cap for fitting (None = all rows)"
    ),
):
    """Train chemogenomic auxiliary model on ChEMBL activity rows (pChEMBL)."""
    from .config import CHEMBL_DIR
    from .model.train_chembl import train_chembl_aux_model

    result = train_chembl_aux_model(
        activities_path=activities or (CHEMBL_DIR / "chembl_pathway_activities.csv"),
        out_dir=out_dir,
        max_rows=max_rows,
    )
    m = result["metrics"]
    rprint(
        f"[green]ChEMBL aux model trained[/green] {result['model_path']}\n"
        f"  rows={m['n_rows']:,}  MAE(pChEMBL)={m['mae_pchembl']:.3f}  "
        f"R²={m['r2_pchembl']:.3f}"
    )


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


@app.command("biomarker")
def biomarker_cmd(
    disease: str = typer.Argument(..., help="Disease name / atlas id / alias"),
    location: Optional[str] = typer.Option(
        None,
        "--location",
        "-l",
        help="Any anatomic site of affected cells (lung, brain, left breast, gut, …)",
    ),
    top: int = typer.Option(50, help="Top-N VOCs by |Δppb| (default 50)"),
    stage: Optional[str] = typer.Option(None, help="Tumor stage, e.g. II, IIIA, IV"),
    genes: Optional[str] = typer.Option(
        None, help="Comma-separated mutated genes, e.g. KRAS,TP53"
    ),
    affected_fraction: Optional[float] = typer.Option(
        None, help="Density of affected cells at location (0–1)"
    ),
    affected_activity: Optional[float] = typer.Option(
        None, help="Metabolic activity of affected cells (≥0)"
    ),
    mode: str = typer.Option("hybrid", help="physiology | hybrid | legacy"),
    age: Optional[float] = typer.Option(None),
    sex: Optional[str] = typer.Option(None),
    smoking: Optional[str] = typer.Option(None, help="never|former|current"),
    comorbidities: Optional[str] = typer.Option(
        None,
        "--comorbidities",
        "-c",
        help="Comma-separated comorbidities, e.g. obesity,heart_disease",
    ),
    comorbidity_weight: float = typer.Option(
        0.65, help="Relative weight of comorbidity priors (0–1.5)"
    ),
    metastatic: bool = typer.Option(False),
    out_dir: Optional[Path] = typer.Option(
        None, help="Write top-50 biomarker CSV/JSON/plots here"
    ),
    no_opentargets: bool = typer.Option(False, help="Disable Open Targets enrichment"),
    no_explain: bool = typer.Option(
        False, help="Skip VOC mechanism explanations (pathways/cells/genes)"
    ),
):
    """
    Whole-body exhaled biomarker prediction: disease + any cell location → top VOCs (ppb).

    Models ~100 atlas diseases × full-body tissue map × 50-VOC panel.
    Comorbidities fuse pathway bias, VOC priors, and cell-state modulation.
    """
    from .biomarker import ExhaleBiomarkerEngine
    from .viz.report import save_biomarker_report

    engine = ExhaleBiomarkerEngine(
        use_opentargets=not no_opentargets, reload_knowledge=True
    )
    mode_norm = mode if mode in {"physiology", "hybrid", "legacy"} else "hybrid"
    comorb = (
        [x.strip() for x in comorbidities.split(",") if x.strip()]
        if comorbidities
        else None
    )
    report = engine.predict(
        disease,
        location=location,
        top_n=top,
        stage=stage,
        genes=[g.strip().upper() for g in genes.split(",")] if genes else None,
        affected_fraction=affected_fraction,
        affected_activity=affected_activity,
        mode=mode_norm,
        sex=sex,
        age_years=age,
        smoking_status=smoking,
        metastatic=metastatic,
        explain=not no_explain,
        comorbidities=comorb,
        comorbidity_weight=comorbidity_weight,
    )
    comorb_label = ""
    meta = report.result.bundle.metadata or {}
    if meta.get("comorbidities"):
        comorb_label = " + " + "+".join(
            str(c.get("name") or c.get("disease_id")) for c in meta["comorbidities"]
        )
    table = Table(
        title=(
            f"ExhalePath Biomarker · {report.disease_name}{comorb_label} @ "
            f"{report.location.get('name') or location}"
        )
    )
    table.add_column("#", justify="right")
    table.add_column("VOC")
    table.add_column("Healthy ppb", justify="right")
    table.add_column("Predicted ppb", justify="right")
    table.add_column("Δ ppb", justify="right")
    table.add_column("Fold", justify="right")
    table.add_column("Conf", justify="right")
    for i, p in enumerate(report.top_vocs, 1):
        table.add_row(
            str(i),
            p.name,
            f"{p.healthy_ppb:.2f}",
            f"{p.predicted_ppb:.2f}",
            f"{p.delta_ppb:+.2f}",
            f"{p.fold_change:.2f}x",
            f"{p.confidence:.2f}",
        )
    rprint(table)
    rprint(
        f"[dim]Modeled {report.n_vocs_modeled} VOCs · showing top {len(report.top_vocs)} "
        f"by |Δppb| · {report.model_version}[/dim]"
    )
    for note in report.notes[:6]:
        rprint(f"[dim]• {note}[/dim]")
    if report.mechanisms:
        rprint("[cyan]Why (mechanisms)[/cyan]")
        for m in report.mechanisms[:5]:
            rprint(f"  • {m.get('why')}")
    if out_dir:
        paths = save_biomarker_report(report, out_dir)
        rprint(
            "[green]Wrote biomarker report:[/green]",
            json.dumps({k: str(v) for k, v in paths.items()}, indent=2),
        )


@app.command("list-diseases")
def list_diseases():
    """List curated disease atlas entries."""
    from .knowledge.loader import KnowledgeBase

    kb = KnowledgeBase()
    for d in kb.diseases.values():
        aliases = ", ".join(d.get("aliases", [])[:4])
        rprint(f"[bold]{d['disease_id']}[/bold] — {d['name']}  ({aliases})")


@app.command("list-locations")
def list_locations():
    """List whole-body anatomic locations for affected-cell placement."""
    from .body.tissues import WholeBodyMap

    body = WholeBodyMap()
    for t in body.list_locations():
        n = t.get("n_census_healthy_cells") or 0
        rprint(
            f"[bold]{t['tissue_id']}[/bold] — {t['name']}  "
            f"(Census healthy cells: {n:,})"
        )


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
