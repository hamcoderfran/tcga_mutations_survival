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
    name="voc",
    help=(
        "voc — exhaled VOC biomarker prediction.\n\n"
        "Quick start:\n"
        '  voc "depression" -l brain -c obesity --age 24 --sex male\n'
        '  voc stack "depression" -l brain -c obesity --age 24 --sex male\n'
        '  voc patient "35M schizophrenia, smokes, on olanzapine, BMI 32"\n'
        '  voc ask                    # interactive question fields\n'
        '  voc nl "24yo obese male with depression"   # optional light LLM\n'
        '  voc "lung adenocarcinoma" -l lung --stage II --genes KRAS,TP53\n'
        "  voc list-diseases\n"
        "  voc --help"
    ),
    add_completion=False,
    no_args_is_help=True,
)


def _disease_catalog() -> list[dict]:
    from .knowledge.loader import KnowledgeBase

    return list(KnowledgeBase().diseases.values())


def _run_biomarker_from_slots(
    slots,
    *,
    out_dir: Optional[Path] = None,
    no_opentargets: bool = False,
    no_explain: bool = False,
    no_save: bool = False,
    top_display: int = 15,
):
    """Shared path: QuerySlots → ExhaleBiomarkerEngine → visual + comprehensive report."""
    from .biomarker import ExhaleBiomarkerEngine
    from .viz.dashboard import (
        default_report_dir,
        print_comprehensive_console,
        save_visual_dashboard,
    )

    miss = slots.missing_required()
    if miss:
        raise typer.BadParameter(f"Missing required fields: {', '.join(miss)}")

    engine = ExhaleBiomarkerEngine(
        use_opentargets=not no_opentargets, reload_knowledge=True
    )
    kwargs = slots.to_biomarker_kwargs()
    report = engine.predict(**kwargs, explain=not no_explain)
    print_comprehensive_console(report, top_display=top_display)

    if not no_save:
        dest = Path(out_dir) if out_dir else default_report_dir(
            report.disease_name, slots.location
        )
        paths = save_visual_dashboard(report, dest)
        rprint(f"[green]Full visual report:[/green] {paths.get('report_html')}")
        rprint(f"[green]Dashboard:[/green] {paths.get('dashboard')}")
        rprint(f"[dim]All artifacts → {dest}[/dim]")
    return report


@app.command("ask")
def ask_cmd(
    nl: Optional[str] = typer.Option(
        None,
        "--nl",
        help="Optional natural-language vignette to seed the questionnaire",
    ),
    llm: str = typer.Option(
        "auto",
        "--llm",
        help="NL backend: auto | rules | ollama | openai (very light; falls back to rules)",
    ),
    yes: bool = typer.Option(
        False, "--yes", "-y", help="Skip confirmation prompt"
    ),
    no_prompt: bool = typer.Option(
        False,
        "--no-prompt",
        help="Do not prompt for missing fields (NL/flags only; fails if disease missing)",
    ),
    out_dir: Optional[Path] = typer.Option(None, help="Write biomarker report here"),
    no_opentargets: bool = typer.Option(False),
    no_explain: bool = typer.Option(False),
):
    """
    Interactive question fields: disease, comorbidities, location, age, sex, …

    Optionally seed from natural language (``--nl``). Use ``--llm ollama`` for a
    tiny local model (default ``qwen2.5:0.5b``), ``--llm openai`` for an API, or
    ``--llm rules`` for zero-weight pattern parsing. Default ``auto`` tries
    Ollama → OpenAI key → rules.
    """
    from .nl import (
        QuerySlots,
        confirm_slots,
        fill_slots_interactively,
        parse_with_optional_llm,
    )

    catalog = _disease_catalog()
    if nl:
        slots = parse_with_optional_llm(nl, llm=llm, disease_catalog=catalog)
        rprint("[cyan]Seeded from natural language[/cyan]")
        for line in slots.summary_lines():
            rprint(f"  {line}")
    else:
        slots = QuerySlots()

    if not no_prompt:
        names = [d.get("name") or d.get("disease_id") or "" for d in catalog]
        slots = fill_slots_interactively(slots, disease_catalog=names)
    elif slots.missing_required():
        raise typer.BadParameter(
            "disease is required; pass --nl that names a disease or omit --no-prompt"
        )

    if not yes and not confirm_slots(slots):
        rprint("[yellow]Cancelled[/yellow]")
        raise typer.Exit(code=0)

    _run_biomarker_from_slots(
        slots,
        out_dir=out_dir,
        no_opentargets=no_opentargets,
        no_explain=no_explain,
    )


@app.command("nl")
def nl_cmd(
    text: str = typer.Argument(..., help="Natural-language clinical vignette"),
    llm: str = typer.Option(
        "auto",
        "--llm",
        help="auto | rules | ollama | openai — very light extractors; rules always available",
    ),
    prompt: bool = typer.Option(
        False,
        "--prompt/--no-prompt",
        help="Fill missing fields interactively after parsing (default: only if disease missing)",
    ),
    force_prompt: bool = typer.Option(
        False,
        "--ask",
        help="Always open the questionnaire after NL parse (confirm/edit fields)",
    ),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation"),
    out_dir: Optional[Path] = typer.Option(None),
    no_opentargets: bool = typer.Option(False),
    no_explain: bool = typer.Option(False),
    show_slots: bool = typer.Option(
        False, "--show-slots", help="Print parsed JSON slots and exit"
    ),
):
    """
    Natural-language input → disease / comorbidities / location / … → biomarker.

    Default path uses a **very light** optional LLM (local Ollama ``qwen2.5:0.5b``
    or OpenAI-compatible) when available; otherwise a zero-weight rule parser.
    If disease is still missing, falls back to the interactive questionnaire.
    """
    from .nl import (
        confirm_slots,
        fill_slots_interactively,
        parse_with_optional_llm,
    )

    catalog = _disease_catalog()
    slots = parse_with_optional_llm(text, llm=llm, disease_catalog=catalog)

    if show_slots:
        rprint(slots.model_dump_json(indent=2))
        raise typer.Exit(code=0)

    need_ask = force_prompt or prompt or bool(slots.missing_required())
    if need_ask:
        names = [d.get("name") or d.get("disease_id") or "" for d in catalog]
        slots = fill_slots_interactively(slots, disease_catalog=names)

    if slots.missing_required():
        raise typer.BadParameter(
            "Could not resolve disease from text; re-run with --ask or use `voc ask`"
        )

    if not yes:
        if not confirm_slots(slots):
            rprint("[yellow]Cancelled[/yellow]")
            raise typer.Exit(code=0)

    _run_biomarker_from_slots(
        slots,
        out_dir=out_dir,
        no_opentargets=no_opentargets,
        no_explain=no_explain,
    )


@app.command("build-corpus")
def build_corpus_cmd(
    max_projects: Optional[int] = typer.Option(
        8,
        help="(legacy) Limit GDC projects — ignored; real corpus is used.",
    ),
    no_reactome: bool = typer.Option(False, help="Unused (legacy flag)"),
    offline_demo: bool = typer.Option(
        False,
        help="DISABLED — synthetic corpora raise unless VOC_ALLOW_SYNTHETIC=1",
    ),
    demo_cases: int = typer.Option(400, help="Unused (legacy flag)"),
    out_dir: Path = typer.Option(PROCESSED_DIR, help="Output directory"),
):
    """Build REAL breath VOC training corpus (redirects to build-real-corpus)."""
    paths = build_training_corpus(
        projects=DEFAULT_GDC_PROJECTS[:max_projects] if max_projects else None,
        out_dir=out_dir,
        expand_reactome=not no_reactome,
        max_projects=max_projects,
        offline_demo=offline_demo,
        demo_cases_per_project=demo_cases,
    )
    rprint("[green]Real corpus ready[/green]")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")


@app.command("build-real-corpus")
def build_real_corpus_cmd(
    out_dir: Path = typer.Option(PROCESSED_DIR, help="Output directory"),
    keep_synthetic: bool = typer.Option(
        False, help="Do not delete leftover synthetic mutation_edges.csv"
    ),
):
    """
    Build training tables from measured MW/Sci Data cohorts + literature panels only.

    Completely excludes synthetic mechanistic+noise VOC labels.
    """
    from .ingest.real_breath_corpus import build_real_training_corpus

    paths = build_real_training_corpus(
        out_dir=out_dir, remove_synthetic=not keep_synthetic
    )
    man = json.loads(Path(paths["manifest"]).read_text())
    rprint("[green]Real-breath corpus ready[/green]")
    rprint(f"  diseases: {man.get('diseases')}")
    rprint(f"  n_targets: {man.get('n_voc_targets')}  synthetic={man.get('synthetic_labels')}")
    rprint(f"  origins: {man.get('label_origins')}")
    for k, v in paths.items():
        rprint(f"  {k}: {v}")


@app.command("eval-priority10")
def eval_priority10_cmd(
    out_dir: Path = typer.Option(Path("runs/priority10_eval")),
    top: int = typer.Option(15),
):
    """Evaluate priority pulmonary/infectious diseases against real literature panels."""
    from .ingest.real_breath_corpus import PRIORITY_DISEASES
    from .biomarker import ExhaleBiomarkerEngine

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_path = Path("data/real_breath/literature_panels/priority10_voc_panels.json")
    panels = {
        p["disease_id"]: p
        for p in json.loads(panel_path.read_text()).get("panels", [])
    }
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    rows = []
    for did in PRIORITY_DISEASES:
        panel = panels.get(did) or {}
        expect = panel.get("measured_log2fc") or {}
        if not expect:
            continue
        loc = panel.get("location") or None
        report = engine.predict(
            did, location=loc, top_n=50, mode="hybrid", explain=False
        )
        pred = {p.voc_id: p.log2_fold_change for p in report.result.bundle.predictions}
        hits = 0
        total = 0
        for voc, exp in expect.items():
            if voc not in pred:
                continue
            total += 1
            # directional agreement
            if (float(exp) >= 0 and pred[voc] >= 0) or (float(exp) < 0 and pred[voc] < 0):
                hits += 1
        acc = hits / total if total else None
        rows.append(
            {
                "disease_id": did,
                "n_vocs_scored": total,
                "directional_accuracy": acc,
                "n_panel_vocs": len(expect),
            }
        )
        mark = f"{acc:.0%}" if acc is not None else "n/a"
        rprint(f"  {did:22} dir_acc={mark}  n={total}")
    overall = [r["directional_accuracy"] for r in rows if r["directional_accuracy"] is not None]
    mean_acc = sum(overall) / len(overall) if overall else 0.0
    summary = {
        "n_diseases": len(rows),
        "mean_directional_accuracy": mean_acc,
        "cases": rows,
    }
    (out_dir / "priority10_eval.json").write_text(json.dumps(summary, indent=2))
    rprint(
        f"[bold]Priority-10 mean directional accuracy: {mean_acc:.1%}[/bold] "
        f"→ {out_dir / 'priority10_eval.json'}"
    )


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
    # also copy into packaged data models for runtime when distinct
    import shutil

    pkg_models = Path(__file__).resolve().parent / "data" / "models"
    pkg_models.mkdir(parents=True, exist_ok=True)
    for name in ("voc_calibrator.joblib", "voc_calibrator_metrics.json"):
        src = Path(result["model_path"]).parent / name
        dst = pkg_models / name
        if src.exists() and src.resolve() != dst.resolve():
            shutil.copy2(src, dst)
            rprint(f"[dim]copied {name} → {pkg_models}[/dim]")


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


@app.command("eval-stack-holdout")
def eval_stack_holdout_cmd(
    out_dir: Path = typer.Option(Path("runs/stack_holdout")),
    skip_slow: bool = typer.Option(
        False, help="Only adversarial + e2e smoke (skip public/lit/priority10)"
    ),
):
    """
    Holdout + adversarial evaluation of Great Disease Stack + PatientTemplate.

    Compares stack vs hybrid on public breath, literature, and priority-10 panels;
    runs a PatientTemplate break suite; writes STACK_HOLDOUT_REPORT.md with an
    honest SOTA assessment.
    """
    from .eval.stack_holdout import evaluate_stack_holdout

    report = evaluate_stack_holdout(out_dir=out_dir, skip_slow=skip_slow)
    sota = report.get("sota_assessment") or {}
    brk = report.get("patient_break") or {}
    pb = report.get("public_breath") or {}
    lit = report.get("literature") or {}
    p10 = report.get("priority10") or {}
    rprint("[bold]Stack holdout[/bold]")
    rprint(
        f"  patient hard pass: {brk.get('hard_pass')}/{brk.get('hard_total')} "
        f"({brk.get('hard_pass_rate')})"
    )
    if pb:
        rprint(
            f"  public_breath dir  hybrid={pb.get('hybrid_mean_directional')} "
            f"stack={pb.get('stack_mean_directional')}"
        )
    if lit:
        rprint(
            f"  literature dir     hybrid={lit.get('hybrid_mean_directional')} "
            f"stack={lit.get('stack_mean_directional')}"
        )
    if p10:
        rprint(
            f"  priority10 dir     hybrid={p10.get('hybrid_mean_directional')} "
            f"stack={p10.get('stack_mean_directional')}"
        )
    rprint(f"  clinical SOTA? {sota.get('is_sota_clinical_breathomics')}")
    rprint(f"  cutting-edge systems stack? {sota.get('is_cutting_edge_systems_stack')}")
    rprint(f"  report → {out_dir / 'STACK_HOLDOUT_REPORT.md'}")


@app.command("eval-patient-diagnostic")
def eval_patient_diagnostic_cmd(
    study: str = typer.Option(
        "ST000883",
        "--study",
        "-s",
        help="Bundled MW study id (ST000883 malaria GC-MS, ST000587 HF EBC)",
    ),
    all_studies: bool = typer.Option(
        False, "--all", help="Run all bundled diagnostic studies"
    ),
    benchmark: bool = typer.Option(
        False,
        "--benchmark",
        help="Compare hybrid/stack/literature × cosine/dot on one study",
    ),
    signature: str = typer.Option(
        "hybrid",
        "--signature",
        help="Disease signature source: hybrid | stack | literature",
    ),
    method: str = typer.Option("cosine", "--method", help="cosine | dot"),
    out_dir: Path = typer.Option(Path("runs/patient_diagnostic")),
    seed: int = typer.Option(42, "--seed"),
    n_splits: int = typer.Option(5, "--n-splits"),
):
    """
    Patient-level GC-MS diagnostic research eval (AUROC / sens / spec / locked splits).

    Uses bundled Metabolomics Workbench patient×VOC matrices. Writes a paper-ready
    pack (ROC, TRIPOD+AI checklist stub, split SHA256). Research enablement only —
    not a clinical validation claim.
    """
    from .eval.patient_diagnostic import (
        evaluate_patient_diagnostic,
        run_multi_study_diagnostic,
        run_signature_benchmark,
    )

    if signature not in {"hybrid", "stack", "literature"}:
        raise typer.BadParameter("signature must be hybrid|stack|literature")
    if method not in {"cosine", "dot"}:
        raise typer.BadParameter("method must be cosine|dot")

    if benchmark:
        payload = run_signature_benchmark(study_id=study, out_dir=out_dir / "benchmark")
        rprint("[bold]Signature benchmark[/bold]")
        for row in payload.get("rows") or []:
            if row.get("ok"):
                rprint(
                    f"  {row['signature']}/{row['method']}: "
                    f"AUROC={row.get('auroc'):.3f} nested={row.get('nested_auroc'):.3f}"
                )
            else:
                rprint(f"  {row['signature']}/{row['method']}: ERROR {row.get('error')}")
        rprint(f"  report → {out_dir / 'benchmark' / 'SIGNATURE_BENCHMARK.md'}")
        return

    if all_studies:
        summary = run_multi_study_diagnostic(
            signature_source=signature,  # type: ignore[arg-type]
            out_dir=out_dir,
        )
        rprint("[bold]Multi-study patient diagnostic[/bold]")
        for s in summary.get("studies") or []:
            if s.get("ok"):
                rprint(
                    f"  {s['study_id']} ({s['disease_id']}): "
                    f"AUROC={s.get('auroc')} nested={s.get('nested_auroc')} "
                    f"logistic_nested={s.get('logistic_nested_auroc')}"
                )
            else:
                rprint(f"  {s['study_id']}: ERROR {s.get('error')}")
        rprint(f"  summary → {out_dir / 'MULTI_STUDY_SUMMARY.md'}")
        return

    report = evaluate_patient_diagnostic(
        study_id=study,
        signature_source=signature,  # type: ignore[arg-type]
        score_method=method,  # type: ignore[arg-type]
        seed=seed,
        n_splits=n_splits,
        out_dir=out_dir / study,
    )
    m = report.get("metrics") or {}
    nested = report.get("nested") or {}
    rprint("[bold]Patient-level GC-MS diagnostic[/bold]")
    rprint(f"  study={report.get('study_id')} disease={report.get('disease_id')}")
    rprint(
        f"  n={m.get('n_subjects')} (pos={m.get('n_positive')} neg={m.get('n_negative')})"
    )
    rprint(
        f"  AUROC={m.get('auroc')}  CI95={m.get('auroc_ci95')}  AUPRC={m.get('auprc')}"
    )
    rprint(
        f"  sens/spec={m.get('sensitivity')}/{m.get('specificity')}  "
        f"confusion={m.get('confusion')}"
    )
    rprint(
        f"  nested sig AUROC={nested.get('mean_test_auroc')}  "
        f"optimism_gap={nested.get('optimism_gap')}"
    )
    log_b = (nested.get("logistic_baseline") or {}).get("mean_test_auroc")
    rprint(f"  nested logistic AUROC={log_b}")
    rprint(f"  split sha256={(nested.get('content_sha256') or '')[:20]}…")
    rprint(f"  report → {out_dir / study / 'PATIENT_DIAGNOSTIC_REPORT.md'}")


@app.command("lock-split")
def lock_split_cmd(
    study: str = typer.Option("ST000883", "--study", "-s"),
    strategy: str = typer.Option("stratified_kfold", "--strategy"),
    n_splits: int = typer.Option(5, "--n-splits"),
    seed: int = typer.Option(42, "--seed"),
    out: Path = typer.Option(
        Path("data/knowledge/gcms_splits/ST000883_split_v1.json"), "--out"
    ),
):
    """Write a preregistration-style locked patient-level split manifest (SHA256)."""
    from .gcms import load_mw_patient_matrix, lock_split, verify_split_manifest

    matrix = load_mw_patient_matrix(study)
    # default out path per study if user left the malaria default while changing study
    if out.name.startswith("ST000883") and study != "ST000883":
        out = Path("data/knowledge/gcms_splits") / f"{study}_split_v1.json"
    manifest = lock_split(
        matrix, strategy=strategy, n_splits=n_splits, seed=seed, out_path=out
    )
    problems = verify_split_manifest(manifest, matrix)
    rprint(f"[bold]Locked split[/bold] → {out}")
    rprint(f"  strategy={manifest.strategy} folds={manifest.n_splits} seed={manifest.seed}")
    rprint(f"  sha256={manifest.content_sha256}")
    if problems:
        rprint(f"  [yellow]integrity issues:[/yellow] {problems}")
    else:
        rprint("  integrity: OK")


@app.command("score-sample")
def score_sample_cmd(
    disease: str = typer.Argument(..., help="Disease id/name for signature"),
    voc_json: str = typer.Option(
        ...,
        "--vocs",
        help='JSON object of voc_id→value (log2fc or relative intensity), e.g. \'{"acetone":0.4}\'',
    ),
    signature: str = typer.Option("hybrid", "--signature"),
    method: str = typer.Option("cosine", "--method"),
):
    """Score an observed VOC vector against a disease signature (research template match)."""
    import json as _json

    from .gcms import disease_signature, score_observed_vector

    observed = _json.loads(voc_json)
    if not isinstance(observed, dict):
        raise typer.BadParameter("--vocs must be a JSON object")
    sig = disease_signature(disease, source=signature)  # type: ignore[arg-type]
    result = score_observed_vector(observed, sig, method=method)  # type: ignore[arg-type]
    rprint("[bold]Sample signature score[/bold]")
    rprint(f"  disease={disease} signature={signature} method={method}")
    rprint(f"  score={result.get('score')} overlap={result.get('n_overlap')}")
    if result.get("vocs_used"):
        rprint(f"  vocs={', '.join(result['vocs_used'][:12])}")


@app.command("eval-scidata-samples")
def eval_scidata_samples_cmd(
    cohort: str = typer.Option(
        "asthma",
        "--cohort",
        "-c",
        help="Positive cohort for one-vs-rest: asthma | copd | bronchiectasis",
    ),
    all_cohorts: bool = typer.Option(
        False, "--all", help="Run all Sci Data pulmonary cohorts"
    ),
    out_dir: Path = typer.Option(Path("runs/scidata_samples")),
    signature: str = typer.Option("hybrid", "--signature"),
    seed: int = typer.Option(42, "--seed"),
    n_splits: int = typer.Option(5, "--n-splits"),
    min_detect_frac: float = typer.Option(0.3, "--min-detect-frac"),
):
    """
    Per-sample Sci Data 2024 peak-table diagnostic (not cohort means).

    One-vs-rest across Asthma/COPD/Bronchiectasis (no healthy arm). Writes
    stratified AUCs (age/sex), blank/detection filter report, and a paper pack zip.
    """
    from .eval.scidata_samples import (
        evaluate_all_scidata_cohorts,
        evaluate_scidata_samples,
    )

    if cohort not in {"asthma", "copd", "bronchiectasis"}:
        raise typer.BadParameter("cohort must be asthma|copd|bronchiectasis")

    if all_cohorts:
        summary = evaluate_all_scidata_cohorts(
            out_dir=out_dir,
            signature_source=signature,
            seed=seed,
            n_splits=n_splits,
            min_detect_frac=min_detect_frac,
        )
        rprint("[bold]Sci Data per-sample (all cohorts)[/bold]")
        for row in summary.get("cohorts") or []:
            if row.get("ok"):
                rprint(
                    f"  {row['cohort']}: AUROC={row.get('auroc')} "
                    f"nested={row.get('nested_auroc')} n={row.get('n')}"
                )
            else:
                rprint(f"  {row['cohort']}: ERROR {row.get('error')}")
        rprint(f"  summary → {out_dir / 'SCIDATA_ALL_SUMMARY.md'}")
        return

    report = evaluate_scidata_samples(
        positive_cohort=cohort,
        out_dir=out_dir / cohort,
        signature_source=signature,
        seed=seed,
        n_splits=n_splits,
        min_detect_frac=min_detect_frac,
    )
    m = report.get("metrics") or {}
    nested = report.get("nested") or {}
    pack = report.get("paper_pack") or {}
    rprint("[bold]Sci Data per-sample diagnostic[/bold]")
    rprint(f"  cohort={cohort} study={report.get('study_id')}")
    rprint(
        f"  n={m.get('n_subjects')} (pos={m.get('n_positive')} neg={m.get('n_negative')})"
    )
    rprint(
        f"  logistic AUROC={m.get('auroc')}  nested={nested.get('mean_test_auroc')}"
    )
    strata = (report.get("stratified_auroc") or {}).get("strata") or {}
    rprint(f"  stratified keys={list(strata.keys()) or ['none']}")
    rprint(f"  filter={(report.get('filter') or {}).get('note')}")
    if pack.get("zip_path"):
        rprint(f"  paper pack → {pack['zip_path']}")
    rprint(f"  report → {out_dir / cohort / 'SCIDATA_SAMPLE_EVAL.md'}")


@app.command("export-metabolights")
def export_metabolights_cmd(
    study: str = typer.Option(
        "ST000883",
        "--study",
        "-s",
        help="Bundled MW study id, or Sci Data cohort via scidata:<asthma|copd|bronchiectasis>",
    ),
    out_dir: Path = typer.Option(Path("runs/metabolights_export")),
):
    """Export mzTab-M-like + ISA-Tab scaffold for MetaboLights-oriented deposit."""
    from .gcms import export_metabolights_bundle, load_mw_patient_matrix, load_scidata_ovr_matrix

    if study.startswith("scidata:"):
        cohort = study.split(":", 1)[1]
        matrix = load_scidata_ovr_matrix(cohort, mapped_vocs_only=True)
    else:
        matrix = load_mw_patient_matrix(study)
    paths = export_metabolights_bundle(matrix, out_dir / matrix.study_id)
    rprint("[bold]MetaboLights-oriented export[/bold]")
    rprint(f"  study={matrix.study_id} disease={matrix.disease_id}")
    for k, p in paths.items():
        rprint(f"  {k} → {p}")


@app.command("export-paper-pack")
def export_paper_pack_cmd(
    run_dir: Path = typer.Argument(..., help="Directory with diagnostic/report artifacts"),
    out: Path = typer.Option(
        None, "--out", help="Output zip path (default: <run_dir>/<name>_paper_pack.zip)"
    ),
):
    """Zip figures + Methods + literature overlay + split hash into one paper pack."""
    from .gcms import export_paper_pack

    info = export_paper_pack(run_dir, out_zip=out)
    rprint("[bold]Paper pack[/bold]")
    rprint(f"  zip → {info['zip_path']}")
    rprint(f"  files={info['n_files']}  zip_sha256={info['zip_sha256'][:24]}…")
    rprint(f"  split_sha256={info.get('split_content_sha256') or '—'}")


@app.command("eval-coverage")
def eval_coverage_cmd(
    out_dir: Path = typer.Option(Path("runs/coverage_audit")),
    offline: bool = typer.Option(True, help="Skip live network expansion (default True)"),
    no_expand: bool = typer.Option(False, help="Audit only; do not expand catalogs"),
):
    """Secure open-VOC coverage audit (VOLATILOME 99% target) + anti-poisoning checks."""
    from .eval.coverage_audit import run_coverage_audit

    report = run_coverage_audit(expand=not no_expand, offline=offline)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "COVERAGE_AUDIT.json").write_text(
        __import__("json").dumps(report, indent=2) + "\n"
    )
    src_md = Path("data/knowledge/COVERAGE_AUDIT.md")
    if src_md.exists():
        (out_dir / "COVERAGE_AUDIT.md").write_text(src_md.read_text())
    m = report["metrics"]
    rprint("[bold]VOC coverage audit[/bold]")
    rprint(
        f"  open-compound coverage: [green]{m['open_compound_coverage_pct']}%[/green] "
        f"({m['open_compound_secured_n']}/{m['open_compound_universe_n']}) "
        f"≥99%={m['meets_99pct_open_compound_target']}"
    )
    rprint(
        f"  studies curated={m.get('curated_public_studies')} "
        f"expanded={m.get('expanded_public_studies')} "
        f"(+{m.get('mw_discovered_additional')})"
    )
    rprint(f"  integrity files hashed: {report['security'].get('n_hashed_files')}")
    rprint(f"  report: {out_dir / 'COVERAGE_AUDIT.md'}")


@app.command("eval-lit-compare")
def eval_lit_compare_cmd(
    out_dir: Path = typer.Option(Path("runs/lit_compare")),
    demo_max_patients: Optional[int] = typer.Option(
        200, help="Patients for demographics PCA (default 200 for speed)"
    ),
    n_diseases: int = typer.Option(50, help="Diseases in connection suite"),
):
    """Literature concordance + demographics PCA + disease VOC connections."""
    from .eval.lit_compare import run_lit_demo_disease100

    report = run_lit_demo_disease100(
        out_dir=out_dir,
        demo_max_patients=demo_max_patients,
        n_diseases=n_diseases,
    )
    lit = report["literature"]
    demo = report["demographics_pca"]
    d100 = report["disease100"]
    rprint("[bold]Literature / demographics / disease suite[/bold]")
    rprint(
        f"  lit concordance: [green]{lit.get('mean_concordance_pct')}%[/green] "
        f"({lit.get('n_diseases')} diseases)"
    )
    sil = demo.get("silhouette") or {}
    rprint(
        f"  demo PCA sil smoking={sil.get('demo_by_smoking')} "
        f"category={sil.get('demo_by_category')}"
    )
    rprint(f"  disease suite: {d100.get('n_diseases')} diseases")
    rprint(f"  report: {out_dir / 'LITERATURE_COMPARE.md'}")


@app.command("eval-malaria-diagnostic")
def eval_malaria_diagnostic_cmd(
    out_dir: Path = typer.Option(Path("runs/malaria_diagnostic")),
    feature_map: str = typer.Option(
        "malaria_lit",
        "--feature-map",
        help="atlas | malaria_lit (terpene/thioether/alkane remaps)",
    ),
    max_features: int = typer.Option(8, "--max-features"),
    seed: int = typer.Option(42, "--seed"),
    n_splits: int = typer.Option(5, "--n-splits"),
    external_matrix: Optional[Path] = typer.Option(
        None, "--external-matrix", help="Optional external patient×VOC CSV"
    ),
    external_labels: Optional[Path] = typer.Option(
        None,
        "--external-labels",
        help="Optional external labels CSV (0/1); use bundled csiro_chmi_labels.csv",
    ),
    loso: bool = typer.Option(
        False,
        "--loso",
        help="Leave-one-study-out vs ST000883 when external intensity matrix is provided",
    ),
    pool_external: bool = typer.Option(
        False,
        "--pool-external",
        help="Also pool external with ST000883 for nested metrics (default: pool only when not --loso)",
    ),
):
    """
    Malaria upgrade pack: literature peak remaps + nested sparse vs transferable
    signature + fixed-sens curves + learning curve / CSIRO labels / LOSO.
    """
    from .eval.malaria_diagnostic import evaluate_malaria_diagnostic

    if feature_map not in {"atlas", "malaria_lit"}:
        raise typer.BadParameter("feature_map must be atlas|malaria_lit")
    report = evaluate_malaria_diagnostic(
        out_dir=out_dir,
        feature_map=feature_map,
        max_features=max_features,
        seed=seed,
        n_splits=n_splits,
        external_matrix=external_matrix,
        external_labels=external_labels,
        loso=loso,
        pool_external=pool_external,
    )
    t = report.get("transferable_signature") or {}
    sk = report.get("fit_on_cohort_sparse_kbest") or {}
    cmp_ = report.get("comparison_to_baseline_atlas_map") or {}
    rprint("[bold]Malaria diagnostic upgrade[/bold]")
    rprint(
        f"  features: atlas={cmp_.get('n_voc_features_atlas')} → "
        f"malaria_lit={cmp_.get('n_voc_features_malaria_lit')} "
        f"(+{len(cmp_.get('newly_mapped_vocs') or [])})"
    )
    rprint(
        f"  transferable nested AUROC={t.get('nested_auroc')} "
        f"gap={t.get('optimism_gap')}"
    )
    rprint(
        f"  sparse nested AUROC={sk.get('mean_test_auroc')} "
        f"gap={sk.get('optimism_gap')} "
        f"feats={sk.get('consensus_features')}"
    )
    for p in (t.get("fixed_sensitivity") or {}).get("points") or []:
        rprint(
            f"  fixed-sens≥{p.get('target_sensitivity')}: "
            f"spec={p.get('specificity')} CI={p.get('specificity_ci95')}"
        )
    csiro = report.get("csiro_chmi") or {}
    loso_r = report.get("leave_one_study_out") or {}
    rprint(
        f"  CSIRO labels bundled={csiro.get('labels_bundled')} "
        f"intensity={csiro.get('intensity_status')}"
    )
    rprint(
        f"  LOSO status={loso_r.get('status')} "
        f"sig={loso_r.get('mean_auroc_signature')} "
        f"sparse={loso_r.get('mean_auroc_sparse')}"
    )
    rprint(f"  report → {out_dir / 'MALARIA_DIAGNOSTIC_UPGRADE.md'}")
    rprint(f"  stop-chasing → {out_dir / 'STOP_CHASING_ST000883.md'}")


@app.command("eval-confounder-ptr")
def eval_confounder_ptr_cmd(
    out_dir: Path = typer.Option(Path("runs/confounder_ptr")),
    skip_scidata: bool = typer.Option(
        False, "--skip-scidata", help="Skip Sci Data age/sex strata block"
    ),
):
    """
    Metadata-rich confounder pack: ST003200 smoking/sex/age + Sci Data GC-MS strata.

    Healthy PTR (n=504) measures VOC→confounder signal — not disease AUROC.
    """
    from .eval.confounder_ptr import evaluate_confounder_ptr

    report = evaluate_confounder_ptr(
        out_dir=out_dir, include_scidata=not skip_scidata
    )
    s = (report.get("st003200") or {}).get("confounder_proxy_auroc") or {}
    rprint("[bold]Confounder / demographics pack[/bold]")
    rprint(f"  smoking current vs never AUROC={s.get('smoking_current_vs_never')}")
    rprint(f"  sex male vs female AUROC={s.get('sex_male_vs_female')}")
    rprint(f"  report → {out_dir / 'CONFOUNDER_PTR.md'}")


@app.command("demo-close")
def demo_close_cmd(
    disease: str = typer.Option(
        "malaria", "--disease", "-d", help="malaria | asthma (fixed closing demos)"
    ),
    note: Optional[str] = typer.Option(
        None, "--note", "-n", help="Patient vignette (default: bundled demo note)"
    ),
    out_dir: Path = typer.Option(Path("runs/closing_demo")),
    no_loso: bool = typer.Option(False, "--no-loso", help="Skip partner LOSO block"),
):
    """
    One demo that closes: paste note → ledger-cited VOCs + optimism gap + paper zip.

    Lead message: cut the cost of wrong VOC panels — not clinical AUROC theater.
    """
    from .eval.closing_demo import run_closing_demo

    if disease not in {"malaria", "asthma"}:
        raise typer.BadParameter("disease must be malaria|asthma")
    report = run_closing_demo(
        disease=disease,  # type: ignore[arg-type]
        note=note,
        out_dir=out_dir,
        include_partner_loso=not no_loso,
    )
    gap = report.get("optimism_gap") or {}
    rprint("[bold]Closing demo — cut the cost of wrong VOC panels[/bold]")
    rprint(f"  disease={report.get('disease')}")
    rprint(
        f"  nested AUROC={gap.get('nested_auroc')}  "
        f"optimism_gap={gap.get('optimism_gap')}"
    )
    rprint(f"  citations={report.get('citation_summary')}")
    pack = report.get("paper_pack") or {}
    rprint(f"  paper zip → {pack.get('zip_path')}")
    rprint(f"  summary → {out_dir / 'CLOSING_DEMO.md'}")


@app.command("diligence-loso")
def diligence_loso_cmd(
    site_a: Optional[Path] = typer.Option(
        None, "--site-a", help="OMNI-style CSV (default: bundled fixture)"
    ),
    site_b: Optional[Path] = typer.Option(None, "--site-b"),
    disease_id: str = typer.Option("malaria", "--disease-id"),
    out: Path = typer.Option(Path("runs/partner_loso/PARTNER_LOSO.json"), "--out"),
    rebuild_fixture: bool = typer.Option(
        False, "--rebuild-fixture", help="Regenerate ST000883 split OMNI fixture"
    ),
):
    """Partner OMNI-style LOSO diligence slide (fixture or your two site CSVs)."""
    from .gcms.omni_partner import (
        build_partner_loso_fixture,
        run_partner_loso_diligence,
    )

    if rebuild_fixture:
        build_partner_loso_fixture()
    report = run_partner_loso_diligence(
        site_a=site_a, site_b=site_b, disease_id=disease_id
    )
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, default=str))
    slide = report.get("buyer_slide") or {}
    rprint("[bold]Partner LOSO diligence[/bold]")
    rprint(f"  {slide.get('headline')}")
    rprint(
        f"  signature LOSO={slide.get('mean_auroc_signature')}  "
        f"sparse LOSO={slide.get('mean_auroc_sparse')}"
    )
    rprint(f"  report → {out}")


@app.command("import-breathvoc")
def import_breathvoc_cmd(
    path: Path = typer.Argument(..., help="BreathVOC JSON or OMNI-style CSV"),
    fmt: str = typer.Option("auto", "--format", help="auto|breathvoc|omni"),
    disease_id: str = typer.Option("malaria", "--disease-id"),
    out_dir: Path = typer.Option(Path("runs/imported_matrix"), "--out-dir"),
):
    """One-click import of partner feature table (BreathVOC JSON or OMNI CSV)."""
    from .gcms.interchange import export_breathvoc
    from .oem.kit import import_feature_table

    if fmt not in {"auto", "breathvoc", "omni"}:
        raise typer.BadParameter("format must be auto|breathvoc|omni")
    packed = import_feature_table(
        path, fmt=fmt, disease_id=disease_id  # type: ignore[arg-type]
    )
    matrix = packed["matrix"]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    export_breathvoc(matrix, out_dir / f"{matrix.study_id}.breathvoc.json")
    matrix.matrix.to_csv(out_dir / f"{matrix.study_id}_matrix.csv")
    matrix.labels.to_frame("label").to_csv(out_dir / f"{matrix.study_id}_labels.csv")
    meta = {k: v for k, v in packed.items() if k != "matrix"}
    (out_dir / "IMPORT_META.json").write_text(json.dumps(meta, indent=2, default=str))
    rprint("[bold]Imported partner matrix[/bold]")
    rprint(f"  format={meta.get('format')} n={meta.get('n_subjects')}×{meta.get('n_vocs')}")
    rprint(f"  → {out_dir}")


@app.command("serve-api")
def serve_api_cmd(
    host: str = typer.Option("127.0.0.1", "--host"),
    port: int = typer.Option(8787, "--port"),
):
    """Procurement/OEM JSON API (score-sample, verify-split, import, demo-close)."""
    from .oem.api_server import serve_api

    serve_api(host=host, port=port)


@app.command("eval-industry-pack")
def eval_industry_pack_cmd(
    out_dir: Path = typer.Option(Path("runs/industry_pack")),
    quick: bool = typer.Option(
        False, "--quick", help="Skip heavier lit-compare / multi-cohort steps"
    ),
):
    """
    Comprehensive industry-standard smoke + benchmark pack with commercial framing.

    Runs nested AUROC / optimism gap / stratified metrics / stack holdout / coverage /
    Sci Data per-sample / paper-pack completeness — then writes BUYER_BRIEF.md with
    honest TAM comps (research enablement, not clinical SOTA).
    """
    from .eval.industry_benchmark import run_industry_pack

    report = run_industry_pack(out_dir=out_dir, quick=quick)
    score = report.get("industry_scorecard") or {}
    rprint("[bold]Industry benchmark pack[/bold]")
    for row in score.get("rows") or []:
        rprint(f"  {row['metric']}: {row['value']}  ({row['standard']})")
    rprint(f"  buyer brief → {out_dir / 'BUYER_BRIEF.md'}")
    rprint(f"  full report → {out_dir / 'INDUSTRY_PACK.md'}")


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


@app.command("eval-vision")
def eval_vision_cmd(
    out_dir: Path = typer.Option(Path("runs/vision_eval")),
    top_k: int = typer.Option(15, help="Top-k for hallmark VOC recall"),
    mode: str = typer.Option("hybrid"),
    limit: Optional[int] = typer.Option(
        None, help="Optional profile limit (debug); default runs all 50"
    ),
    profiles: Optional[Path] = typer.Option(
        None, help="Override path to vision_50_profiles.json"
    ),
):
    """
    Full vision suite: 50 disease patient profiles vs measured / high-accuracy GT.

    Prints headline % scores for closeness to the original ExhalePath vision
    (disease + location + patient → ranked VOCs + pathway/cell/site explainability).
    """
    from .eval.vision import evaluate_vision

    report = evaluate_vision(
        top_k=top_k,
        mode=mode,
        out_dir=out_dir,
        profiles_path=profiles,
        limit=limit,
    )
    o = report["overall"]
    rprint("[bold]Vision evaluation (50 diseases)[/bold]")
    rprint(f"  profiles: {o['n_profiles']}")
    rprint(f"  [green]vision fidelity:[/green] {o['vision_fidelity_pct']}%")
    rprint(f"  evidence-backed (A+B): {o['evidence_backed_pct']}%")
    rprint(f"  held-out style (A+B+C): {o['held_out_style_pct']}%")
    rprint(f"  Grade D prior consistency: {o['grade_D_prior_consistency_pct']}%")
    mm = o.get("metric_means") or {}
    rprint(
        "  means — "
        f"dir={_fmt_pct(mm.get('direction'))} "
        f"fold={_fmt_pct(mm.get('fold'))} "
        f"topk={_fmt_pct(mm.get('topk'))} "
        f"pathway={_fmt_pct(mm.get('pathway'))} "
        f"cell={_fmt_pct(mm.get('cell'))} "
        f"site={_fmt_pct(mm.get('site'))}"
    )
    for g, row in (o.get("by_grade") or {}).items():
        rprint(f"  Grade {g} (n={row['n']}): {row['mean_composite_pct']}%")
    rprint(f"  report: {out_dir / 'VISION_EVAL.md'}")
    rprint(f"  json: {out_dir / 'vision_eval.json'}")


@app.command("eval-stress-hard")
def eval_stress_hard_cmd(
    out_dir: Path = typer.Option(Path("runs/stress_hard")),
    profiles: Optional[Path] = typer.Option(
        None, help="Override path to stress_hard_50.json"
    ),
):
    """Adversarial 50-case stress suite — empty inputs, gene traps, alias collisions, SZ suppress."""
    from .eval.stress_hard import evaluate_stress_hard

    report = evaluate_stress_hard(out_dir=out_dir, profiles_path=profiles)
    o = report["overall"]
    rprint("[bold]Stress-hard evaluation (50 adversarial cases)[/bold]")
    rprint(f"  pass rate: [green]{o['pass_pct']}%[/green] ({o['n_passed']}/{o['n_cases']})")
    for cat, row in (o.get("by_category") or {}).items():
        rprint(f"  {cat}: {row['passed']}/{row['n']} ({100*row['pass_rate']:.0f}%)")
    if o.get("failed_ids"):
        rprint(f"  failed: {', '.join(o['failed_ids'][:12])}{'…' if len(o['failed_ids'])>12 else ''}")
    rprint(f"  report: {out_dir / 'STRESS_HARD.md'}")


@app.command("eval-patient-cohort")
def eval_patient_cohort_cmd(
    out_dir: Path = typer.Option(Path("runs/patient_cohort_1000")),
    cohort: Optional[Path] = typer.Option(
        None, help="Override path to patient_cohort_1000.json"
    ),
    max_patients: Optional[int] = typer.Option(
        None, help="Optional limit for smoke tests"
    ),
):
    """
    1000-patient diversity cohort: dense VOC CSVs, breakdown flags,
    PCA/UMAP of VOC profiles and genetic-shift vectors (COPD / bronchitis / lung cancer).
    """
    from .eval.patient_cohort import run_cohort

    report = run_cohort(out_dir=out_dir, cohort_path=cohort, max_patients=max_patients)
    o = report["overall"]
    rprint("[bold]Patient cohort evaluation[/bold]")
    rprint(f"  ok: [green]{o['ok_pct']}%[/green] ({o['n_ok']}/{o['n_patients']})")
    rprint(f"  unresolved: {o.get('n_unresolved')}")
    flags = o.get("breakdown_flag_counts") or {}
    if flags:
        top = ", ".join(f"{k}={v}" for k, v in list(flags.items())[:6])
        rprint(f"  breakdown: {top}")
    dist = o.get("focus_centroid_distances") or {}
    if dist:
        rprint("  focus centroid L2:")
        for k, v in sorted(dist.items()):
            rprint(f"    {k}: {v:.3f}")
    rprint(f"  report: {out_dir / 'COHORT_1000.md'}")
    rprint(f"  dense csv: {out_dir / 'patients_dense.csv'}")
    rprint(f"  figures: {out_dir / 'figures'}")


def _fmt_pct(v) -> str:
    if v is None:
        return "—"
    return f"{100.0 * float(v):.0f}%"


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


def _parse_kv_floats(spec: Optional[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    if not spec:
        return out
    for part in spec.split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        try:
            out[k.strip()] = float(v)
        except ValueError:
            continue
    return out


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
    description: Optional[str] = typer.Option(
        None,
        "--description",
        "-d",
        help="Phenotype/mechanism text for zero-shot (e.g. 'ketotic mitochondrial stress')",
    ),
    pathway_overrides: Optional[str] = typer.Option(
        None,
        help="pathway=score pairs, e.g. ketone_body_metabolism=1.8,lipid_peroxidation=1.4",
    ),
    cell_fractions: Optional[str] = typer.Option(
        None,
        help="cell_state=fraction pairs, e.g. hepatocyte_ketogenic=0.5,oxidative_stress_cell=0.3",
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
    copy_voc_priors: bool = typer.Option(
        False,
        help="Allow ontology NN to copy VOC priors (default: mechanism-only transfer)",
    ),
    out_dir: Optional[Path] = typer.Option(
        None,
        help="Write full visual report here (default: auto runs/voc_<disease>_…)",
    ),
    no_save: bool = typer.Option(
        False, help="Do not write HTML/dashboard/CSV artifacts"
    ),
    no_opentargets: bool = typer.Option(False, help="Disable Open Targets enrichment"),
    no_explain: bool = typer.Option(
        False, help="Skip VOC mechanism explanations (pathways/cells/genes)"
    ),
):
    """
    Whole-body exhaled biomarker prediction: disease + any cell location → top VOCs (ppb).

    One-liner friendly (also the default when you run `voc "disease" …`):
      voc "depression" -l brain -c obesity --age 24 --sex male

    Auto-writes a comprehensive visual pack (HTML + dashboard PNG + CSV/JSON)
    unless --no-save is set.
    """
    from .biomarker import ExhaleBiomarkerEngine
    from .viz.dashboard import (
        default_report_dir,
        print_comprehensive_console,
        save_visual_dashboard,
    )

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
        pathway_overrides=_parse_kv_floats(pathway_overrides) or None,
        cell_state_fractions=_parse_kv_floats(cell_fractions) or None,
        description=description,
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
        copy_voc_priors_from_neighbor=copy_voc_priors,
    )
    print_comprehensive_console(report, top_display=min(15, top))
    if not no_save:
        dest = Path(out_dir) if out_dir else default_report_dir(disease, location)
        paths = save_visual_dashboard(report, dest)
        rprint(f"[green]Full visual report:[/green] {paths.get('report_html')}")
        rprint(f"[green]Dashboard:[/green] {paths.get('dashboard')}")
        rprint(f"[dim]Open REPORT.html in a browser · artifacts → {dest}[/dim]")


def _run_stack_and_report(
    *,
    disease: str,
    location: Optional[str],
    comorbidities: Optional[list[str]],
    age: Optional[float],
    sex: Optional[str],
    genes: Optional[list[str]],
    description: Optional[str],
    smoking: Optional[str],
    top: int,
    out_dir: Optional[Path],
    no_save: bool,
    json_out: bool,
    patient_template: Optional[dict] = None,
):
    from .great_stack import run_great_stack
    from .great_stack.report import (
        default_stack_dir,
        print_stack_console,
        result_to_dict,
        save_stack_report,
    )

    result = run_great_stack(
        disease,
        location=location,
        comorbidities=list(comorbidities or []),
        age=age,
        sex=sex,
        genes=genes,
        description=description,
        smoking=smoking,
        top_n=top,
    )
    payload = result_to_dict(result)
    if patient_template is not None:
        payload["patient_template"] = patient_template
    if json_out:
        typer.echo(json.dumps(payload, indent=2, default=str))
        return result
    print_stack_console(result, top_display=min(15, top))
    if not no_save:
        dest = Path(out_dir) if out_dir else default_stack_dir(disease, location)
        paths = save_stack_report(result, dest)
        if patient_template is not None:
            (dest / "PATIENT_TEMPLATE.json").write_text(
                json.dumps(patient_template, indent=2, default=str)
            )
            rprint(f"[cyan]Patient template:[/cyan] {dest / 'PATIENT_TEMPLATE.json'}")
        rprint(f"[green]Stack HTML report:[/green] {paths.get('html')}")
        rprint(f"[green]Dashboard:[/green] {paths.get('dashboard')}")
        rprint(f"[dim]Full pack → {dest}[/dim]")
    return result


@app.command("stack")
def great_stack_cmd(
    disease: Optional[str] = typer.Argument(
        None, help="Disease name (atlas or novel). Optional if --nl is set."
    ),
    location: Optional[str] = typer.Option(None, "--location", "-l"),
    comorbidity: Optional[list[str]] = typer.Option(None, "--comorbidity", "-c"),
    age: Optional[float] = typer.Option(None, "--age"),
    sex: Optional[str] = typer.Option(None, "--sex"),
    genes: Optional[str] = typer.Option(None, help="Comma-separated genes"),
    description: Optional[str] = typer.Option(None, "--description", "-d"),
    smoking: Optional[str] = typer.Option(None, help="never|former|current"),
    nl: Optional[str] = typer.Option(
        None,
        "--nl",
        help="Naturalistic patient vignette/note → PatientTemplate → stack",
    ),
    llm: str = typer.Option(
        "rules",
        "--llm",
        help="NL backend for --nl: rules | auto | ollama | openai",
    ),
    top: int = typer.Option(20, "--top", "-k"),
    out_dir: Optional[Path] = typer.Option(
        None, "--out-dir", help="Write STACK_REPORT.html pack (default: runs/stack_…)"
    ),
    no_save: bool = typer.Option(False, "--no-save"),
    json_out: bool = typer.Option(False, "--json"),
    show_template: bool = typer.Option(
        False, "--show-template", help="Print parsed PatientTemplate JSON and exit"
    ),
):
    """
    Great Disease Prediction Stack — fuse 20 models across VOC + disease biology.

    Structured flags **or** naturalistic input:
      voc stack "schizophrenia" -l brain --age 35 --sex male
      voc stack --nl "35M with schizophrenia, smokes, on olanzapine, BMI 32"
    """
    template_dump = None
    if nl:
        from .nl import parse_patient_template

        tpl = parse_patient_template(nl, disease_catalog=_disease_catalog(), llm=llm)
        template_dump = tpl.model_dump()
        rprint("[cyan]Parsed patient template[/cyan]")
        for line in tpl.summary_lines():
            rprint(f"  {line}")
        if show_template:
            typer.echo(json.dumps(template_dump, indent=2, default=str))
            raise typer.Exit(code=0)
        if tpl.missing_required():
            raise typer.BadParameter(
                "Could not resolve primary_disease from --nl; try a clearer vignette "
                "or pass disease as the first argument"
            )
        kw = tpl.to_stack_kwargs()
        # CLI flags override template when explicitly provided
        _run_stack_and_report(
            disease=disease or kw["disease"],
            location=location or kw.get("location"),
            comorbidities=list(comorbidity or kw.get("comorbidities") or []),
            age=age if age is not None else kw.get("age"),
            sex=sex or kw.get("sex"),
            genes=(
                [g.strip().upper() for g in genes.split(",") if g.strip()]
                if genes
                else kw.get("genes")
            ),
            description=description or kw.get("description"),
            smoking=smoking or kw.get("smoking"),
            top=top,
            out_dir=out_dir,
            no_save=no_save,
            json_out=json_out,
            patient_template=template_dump,
        )
        return

    if not disease:
        raise typer.BadParameter("Provide a disease argument or --nl vignette")
    if show_template:
        raise typer.BadParameter("--show-template requires --nl")
    gene_list = [g.strip().upper() for g in genes.split(",") if g.strip()] if genes else None
    _run_stack_and_report(
        disease=disease,
        location=location,
        comorbidities=list(comorbidity or []),
        age=age,
        sex=sex,
        genes=gene_list,
        description=description,
        smoking=smoking,
        top=top,
        out_dir=out_dir,
        no_save=no_save,
        json_out=json_out,
    )


@app.command("patient")
def patient_cmd(
    text: Optional[str] = typer.Argument(
        None,
        help="Naturalistic patient vignette, clinic note, or JSON/key=value block",
    ),
    file: Optional[Path] = typer.Option(
        None, "--file", "-f", help="Read patient text from a file"
    ),
    llm: str = typer.Option(
        "auto",
        "--llm",
        help="rules | auto | ollama | openai — enrich gaps when available",
    ),
    engine: str = typer.Option(
        "stack",
        "--engine",
        help="stack (20-model fusion) | biomarker (single ExhalePath engine)",
    ),
    show_template: bool = typer.Option(
        False, "--show-template", help="Print PatientTemplate JSON only"
    ),
    yes: bool = typer.Option(False, "--yes", "-y", help="Skip confirmation"),
    out_dir: Optional[Path] = typer.Option(None, "--out-dir"),
    no_save: bool = typer.Option(False, "--no-save"),
    json_out: bool = typer.Option(False, "--json"),
    top: int = typer.Option(20, "--top", "-k"),
):
    """
    Naturalistic patient input → PatientTemplate → prediction.

    Accepts vignettes, clinic-note sections (CC/HPI/PMH/Meds), key=value lines,
    or JSON. Example:

      voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations"
      voc patient -f note.txt --engine stack
      voc patient --show-template "24yo obese female with depression and insomnia"
    """
    from .nl import confirm_slots, parse_patient_template

    if file is not None:
        text = Path(file).read_text()
    if text is None or not str(text).strip():
        # read stdin if piped
        import sys

        if not sys.stdin.isatty():
            text = sys.stdin.read()
    if not text or not str(text).strip():
        raise typer.BadParameter("Provide vignette text, --file, or stdin")

    tpl = parse_patient_template(
        text, disease_catalog=_disease_catalog(), llm=llm
    )
    rprint("[cyan]Patient template[/cyan]")
    for line in tpl.summary_lines():
        rprint(f"  {line}")

    if show_template:
        typer.echo(tpl.model_dump_json(indent=2))
        raise typer.Exit(code=0)

    if tpl.missing_required():
        raise typer.BadParameter(
            "primary_disease unresolved — include a condition name or diagnosis"
        )

    if not yes:
        slots = tpl.to_query_slots()
        if not confirm_slots(slots):
            rprint("[yellow]Cancelled[/yellow]")
            raise typer.Exit(code=0)

    tpl.top = top
    if engine == "biomarker":
        slots = tpl.to_query_slots()
        # stash description via engine predict kwargs
        from .biomarker import ExhaleBiomarkerEngine
        from .viz.dashboard import (
            default_report_dir,
            print_comprehensive_console,
            save_visual_dashboard,
        )

        kw = tpl.to_biomarker_kwargs()
        eng = ExhaleBiomarkerEngine(use_opentargets=True, reload_knowledge=True)
        report = eng.predict(**kw, explain=True)
        if json_out:
            typer.echo(
                json.dumps(
                    {
                        "patient_template": tpl.model_dump(),
                        "top_vocs": [p.model_dump() for p in report.top_vocs],
                    },
                    indent=2,
                    default=str,
                )
            )
            return
        print_comprehensive_console(report, top_display=min(15, top))
        if not no_save:
            dest = Path(out_dir) if out_dir else default_report_dir(
                report.disease_name, slots.location
            )
            paths = save_visual_dashboard(report, dest)
            (dest / "PATIENT_TEMPLATE.json").write_text(tpl.model_dump_json(indent=2))
            rprint(f"[green]Report:[/green] {paths.get('report_html')}")
        return

    # default: full stack
    kw = tpl.to_stack_kwargs()
    _run_stack_and_report(
        disease=kw["disease"],
        location=kw.get("location"),
        comorbidities=kw.get("comorbidities"),
        age=kw.get("age"),
        sex=kw.get("sex"),
        genes=kw.get("genes"),
        description=kw.get("description"),
        smoking=kw.get("smoking"),
        top=top,
        out_dir=out_dir,
        no_save=no_save,
        json_out=json_out,
        patient_template=tpl.model_dump(),
    )


@app.command("predict-novel")
def predict_novel_cmd(
    disease: str = typer.Argument(..., help="Novel / unseen disease name"),
    location: Optional[str] = typer.Option(
        None, "--location", "-l", help="Affected tissue / organ"
    ),
    description: Optional[str] = typer.Option(
        None, "--description", "-d", help="Phenotype / mechanism description"
    ),
    genes: Optional[str] = typer.Option(
        None, help="Comma-separated mechanism genes, e.g. HMGCS2,CPT1A"
    ),
    pathway_overrides: Optional[str] = typer.Option(
        None, help="pathway=score pairs for first-class mechanism input"
    ),
    cell_fractions: Optional[str] = typer.Option(
        None, help="cell_state=fraction pairs"
    ),
    affected_fraction: Optional[float] = typer.Option(
        0.35, help="Default affected-cell density at location for novel diseases"
    ),
    top: int = typer.Option(15, help="Top-N VOCs to display"),
    mode: str = typer.Option("hybrid"),
    out_dir: Optional[Path] = typer.Option(
        None, help="Write visual report here (default: auto under runs/)"
    ),
    no_save: bool = typer.Option(
        False, help="Do not write HTML/dashboard/CSV artifacts"
    ),
    copy_voc_priors: bool = typer.Option(
        False, help="Allow VOC prior copy from ontology neighbor (off by default)"
    ),
):
    """
    Zero-shot novel-disease path: mechanism-first prediction with full visual report.

    Prefer genes / pathway overrides / tissue / description over bare names.
    VOC priors are not copied from atlas neighbors unless --copy-voc-priors.
    Same output pack as `voc \"disease\" …` (REPORT.html + dashboard.png).
    """
    from .biomarker import ExhaleBiomarkerEngine
    from .viz.dashboard import (
        default_report_dir,
        print_comprehensive_console,
        save_visual_dashboard,
    )

    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    mode_norm = mode if mode in {"physiology", "hybrid", "legacy"} else "hybrid"
    report = engine.predict(
        disease,
        location=location,
        top_n=top,
        genes=[g.strip().upper() for g in genes.split(",")] if genes else None,
        pathway_overrides=_parse_kv_floats(pathway_overrides) or None,
        cell_state_fractions=_parse_kv_floats(cell_fractions) or None,
        description=description,
        affected_fraction=affected_fraction,
        mode=mode_norm,
        explain=True,
        copy_voc_priors_from_neighbor=copy_voc_priors,
    )
    print_comprehensive_console(report, top_display=min(15, top))
    if not no_save:
        dest = Path(out_dir) if out_dir else default_report_dir(disease, location)
        paths = save_visual_dashboard(report, dest)
        rprint(f"[green]Full visual report:[/green] {paths.get('report_html')}")
        rprint(f"[green]Dashboard:[/green] {paths.get('dashboard')}")
        rprint(f"[dim]Open REPORT.html in a browser · artifacts → {dest}[/dim]")


@app.command("eval-zero-shot-reliability")
def eval_zero_shot_reliability_cmd(
    out_dir: Path = typer.Option(Path("runs/zero_shot_reliability")),
    mode: str = typer.Option("hybrid"),
    profiles: Optional[Path] = typer.Option(
        None, help="Override zero_shot_holdout_profiles.json"
    ),
    max_ldo: Optional[int] = typer.Option(
        None, help="Optional cap on leave-disease-out atlas diseases"
    ),
):
    """Held-out rare-disease zero-shot eval + leave-disease-out prior-direction audit."""
    from .eval.zero_shot_reliability import evaluate_zero_shot_reliability

    report = evaluate_zero_shot_reliability(
        out_dir=out_dir,
        profiles_path=profiles,
        mode=mode,
        max_ldo_diseases=max_ldo,
    )
    h = report["holdout"]
    ldo = report["leave_disease_out"]
    rprint("[bold]Zero-shot reliability[/bold]")
    rprint(
        f"  holdout: [green]{100 * h['pass_rate']:.0f}%[/green] "
        f"({h['n_passed']}/{h['n_profiles']})"
    )
    if ldo.get("skipped"):
        rprint(f"  leave-disease-out: skipped ({ldo.get('reason')})")
    else:
        rprint(
            f"  leave-disease-out: [green]{100 * ldo['pass_rate']:.0f}%[/green] "
            f"mean_dir={ldo['mean_directional_accuracy']:.3f} "
            f"({ldo['n_passed']}/{ldo['n_diseases']})"
        )
    rprint(f"  overall pass: {report['pass']}")
    rprint(f"  report: {out_dir / 'ZERO_SHOT_RELIABILITY.md'}")
    if not report["pass"]:
        raise typer.Exit(code=1)


@app.command("eval-implementation-readiness")
def eval_implementation_readiness_cmd(
    out_dir: Path = typer.Option(Path("runs/implementation_readiness")),
    max_cohort: Optional[int] = typer.Option(
        1000, help="Patient cohort size (default full 1000)"
    ),
    skip_pytest: bool = typer.Option(False, help="Skip pytest layer"),
):
    """Comprehensive multi-layer readiness gate across all eval suites."""
    from .eval.implementation_readiness import evaluate_implementation_readiness

    report = evaluate_implementation_readiness(
        out_dir=out_dir,
        max_cohort_patients=max_cohort,
        skip_pytest=skip_pytest,
    )
    rprint("[bold]Implementation readiness[/bold]")
    ready = report["ready_for_implementation"]
    color = "green" if ready else "red"
    rprint(
        f"  ready: [{color}]{ready}[/{color}]  "
        f"score={report['n_gates_passed']}/{report['n_gates']}"
    )
    for k, v in report["gates"].items():
        rprint(f"  {'✓' if v else '✗'} {k}")
    rprint(f"  recommendation: {report['recommendation']}")
    rprint(f"  report: {out_dir / 'IMPLEMENTATION_READINESS.md'}")
    if not ready:
        raise typer.Exit(code=1)


@app.command("eval-zero-shot-hard")
def eval_zero_shot_hard_cmd(
    out_dir: Path = typer.Option(Path("runs/zero_shot_hard_break")),
    mode: str = typer.Option("hybrid"),
    profiles: Optional[Path] = typer.Option(
        None, help="Override zero_shot_hard_break.json"
    ),
):
    """Adversarial hard-break suite for zero-shot (ambiguous, eponym, conflict, traps)."""
    from .eval.zero_shot_hard import evaluate_zero_shot_hard

    report = evaluate_zero_shot_hard(
        out_dir=out_dir, profiles_path=profiles, mode=mode
    )
    rprint("[bold]Zero-shot hard-break[/bold]")
    rprint(
        f"  pass rate: {100 * report['pass_rate']:.0f}% "
        f"({report['n_passed']}/{report['n_profiles']})"
    )
    rprint(f"  failed: {', '.join(report.get('failed_ids') or []) or '—'}")
    for t, row in (report.get('by_trap') or {}).items():
        rprint(
            f"  {t}: {row['n'] - row['failed']}/{row['n']} ({100 * row['pass_rate']:.0f}%)"
        )
    rprint(f"  report: {out_dir / 'ZERO_SHOT_HARD_BREAK.md'}")
    if not report["pass"]:
        raise typer.Exit(code=1)


@app.command("cjd-profile")
def cjd_profile_cmd(
    age: float = typer.Option(62, help="Patient age"),
    sex: str = typer.Option("female", help="female|male|other"),
    genes: str = typer.Option("PRNP", help="Comma-separated genes"),
    out_dir: Path = typer.Option(
        Path("runs/cjd_clinical_profile"), help="Write JSON + Markdown profile"
    ),
):
    """
    Whole clinical VOC tempo for Creutzfeldt–Jakob disease (incubating → terminal).

    Research hypothesis — no validated exhaled VOC signature for CJD exists.
    """
    import sys
    from pathlib import Path as P

    # Prefer installed package script path; fall back to repo scripts/
    repo_scripts = P(__file__).resolve().parents[2] / "scripts"
    if str(repo_scripts) not in sys.path:
        sys.path.insert(0, str(repo_scripts.parent))
    from scripts.cjd_clinical_profile import run_profile

    sex_n = sex if sex in {"female", "male", "other"} else "female"
    gene_list = [g.strip().upper() for g in genes.split(",") if g.strip()]
    profile = run_profile(
        age=age, sex=sex_n, genes=gene_list or ["PRNP"], out_dir=out_dir
    )
    s = profile["summary"]
    rprint("[bold]CJD clinical VOC profile[/bold]")
    rprint(s["interpretation"])
    rprint(
        f"  first hint: {s['first_meaningful_hint_phase']} · "
        f"first obvious: {s['first_obvious_phase']}"
    )
    for ph in profile["phases"]:
        ox = ph["oxidative_panel"]
        rprint(
            f"  • {ph['label']}: [cyan]{ox['verdict']}[/cyan] "
            f"mean|log2FC|={ox['mean_abs_log2fc_oxidative']}"
        )
    rprint(f"[green]Wrote[/green] {out_dir / 'cjd_clinical_profile.md'}")


@app.command("diabetes-profile")
def diabetes_profile_cmd(
    age: float = typer.Option(55, help="Patient age"),
    sex: str = typer.Option("male", help="female|male|other"),
    genes: str = typer.Option("TCF7L2,PPARG", help="Comma-separated genes"),
    out_dir: Path = typer.Option(
        Path("runs/diabetes_clinical_profile"),
        help="Write JSON + Markdown profile",
    ),
):
    """
    Clinical VOC tempo for type 2 diabetes (prediabetes → ketotic) + literature eval.

    Gates on published breath acetone / ketone expectations (PMID:21903721).
    """
    import sys
    from pathlib import Path as P

    repo_scripts = P(__file__).resolve().parents[2] / "scripts"
    if str(repo_scripts) not in sys.path:
        sys.path.insert(0, str(repo_scripts.parent))
    from scripts.diabetes_clinical_profile import run_profile

    sex_n = sex if sex in {"female", "male", "other"} else "male"
    gene_list = [g.strip().upper() for g in genes.split(",") if g.strip()]
    profile = run_profile(
        age=age,
        sex=sex_n,
        genes=gene_list or ["TCF7L2", "PPARG"],
        out_dir=out_dir,
    )
    s = profile["summary"]
    mark = "PASSED" if s["overall_passed"] else "FAILED"
    rprint(f"[bold]Diabetes clinical VOC profile — {mark}[/bold]")
    rprint(s["interpretation"])
    for pid, t in s["tempo"].items():
        af = t["acetone_fold"]
        rprint(
            f"  • {pid}: [cyan]{t['verdict']}[/cyan] "
            f"acetone={af:.2f}× lit={'yes' if t['literature_passed'] else 'no'}"
        )
    lit = s["primary_literature_eval"]
    rprint(
        f"  literature (poorly controlled): acetone_pass={lit['acetone_pass']} "
        f"dir_acc={lit['directional_accuracy']}"
    )
    rprint(f"[green]Wrote[/green] {out_dir / 'diabetes_clinical_profile.md'}")
    if not s["overall_passed"]:
        raise typer.Exit(code=1)


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


def main(argv: Optional[list[str]] = None):
    """
    Entry point for ``voc`` / ``exhalepath``.

    If the first argument is not a known subcommand, treat the invocation as
    ``voc biomarker …`` so users can run::

        voc "depression" -l brain -c obesity
    """
    import sys

    args = list(sys.argv[1:] if argv is None else argv)
    # Known top-level commands (keep in sync with @app.command names)
    commands = {
        "build-corpus",
        "build-real-corpus",
        "train",
        "eval-completion",
        "eval-priority10",
        "integrate-datasources",
        "harvest-public-breath",
        "eval-public-breath",
        "eval-stack-holdout",
        "eval-patient-diagnostic",
        "eval-scidata-samples",
        "export-metabolights",
        "export-paper-pack",
        "eval-coverage",
        "eval-lit-compare",
        "eval-malaria-diagnostic",
        "eval-confounder-ptr",
        "demo-close",
        "diligence-loso",
        "import-breathvoc",
        "serve-api",
        "eval-industry-pack",
        "lock-split",
        "score-sample",
        "harvest-clinical-comorbidity",
        "eval-comorbidity-clinical",
        "eval-vision",
        "eval-stress-hard",
        "eval-patient-cohort",
        "build-mechanism-packs",
        "explain",
        "harvest-chembl",
        "train-chembl",
        "predict",
        "biomarker",
        "predict-novel",
        "stack",
        "patient",
        "eval-zero-shot-reliability",
        "eval-implementation-readiness",
        "eval-zero-shot-hard",
        "ask",
        "nl",
        "cjd-profile",
        "diabetes-profile",
        "list-diseases",
        "list-locations",
        "list-vocs",
        "harvest-census",
        "audit",
        "eval-multisite",
        "help",
    }
    if args and not args[0].startswith("-") and args[0] not in commands:
        # Default to biomarker for disease-first UX
        args = ["biomarker", *args]
        sys.argv = [sys.argv[0], *args]
    app()


if __name__ == "__main__":
    main()
