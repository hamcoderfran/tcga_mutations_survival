# ExhalePath results index

Single entry point for past evaluations and reports on `main`.
Trained calibrators (`data/models/*.joblib`) are **not** modified by eval merges.

## One-line product usage (v1.4)

```bash
pip install -e exhalepath_atlas   # or: bash exhalepath_atlas/install_voc.sh
voc "depression" -l brain -c obesity --age 24 --sex male
# → rich console + runs/voc_*/REPORT.html + dashboard.png + REPORT.md + CSVs
voc eval-implementation-readiness   # research readiness gates
```

See [QUICKSTART.md](../../QUICKSTART.md).

## Headline reports (canonical)

| Result | Path | How to re-run |
|---|---|---|
| Visual review (HTML) | [`EXHALEPATH_REVIEW_REPORT.html`](EXHALEPATH_REVIEW_REPORT.html) | `python scripts/generate_review_report_figures.py` |
| Real-data pipeline | [`REAL_PIPELINE_REPORT.md`](REAL_PIPELINE_REPORT.md) | `python scripts/run_real_pipeline_eval.py` |
| Vision suite (50 diseases) | [`VISION_EVAL_REPORT.md`](VISION_EVAL_REPORT.md) · [`vision_eval_latest.json`](vision_eval_latest.json) | `voc eval-vision` |
| Stress-hard (50 adversarial) | [`STRESS_HARD_REPORT.md`](STRESS_HARD_REPORT.md) | `voc eval-stress-hard` |
| Stress baseline / after | [`STRESS_HARD_BASELINE.md`](STRESS_HARD_BASELINE.md) · [`STRESS_HARD_AFTER.md`](STRESS_HARD_AFTER.md) | — |
| 1000-patient cohort | [`cohort_1000/COHORT_1000.md`](cohort_1000/COHORT_1000.md) | `voc eval-patient-cohort` |
| Cohort embeddings / patients | [`cohort_1000/`](cohort_1000/) | PCA/UMAP PNGs + CSVs |

## Clinical profiles

| Profile | Example write-up | Script |
|---|---|---|
| Type 2 diabetes | [`../examples/diabetes_clinical_profile.md`](../examples/diabetes_clinical_profile.md) | `voc diabetes-profile` |
| Creutzfeldt–Jakob (CJD) | [`../examples/cjd_clinical_profile.md`](../examples/cjd_clinical_profile.md) | `voc cjd-profile` |

## Local run artifacts (gitignored `runs/` when present)

Reproducible under `exhalepath_atlas/runs/` after CLI evals:

- `runs/vision_eval/` · `runs/real_pipeline/` · `runs/review_report/`
- `runs/stress_hard/` · `runs/patient_cohort_1000/`
- `runs/diabetes_clinical_profile/` · `runs/cjd_clinical_profile/`
- `runs/comorbidity_clinical_eval/` · `runs/multisite_eval*` · `runs/completion/`

## Quick regression commands

```bash
cd exhalepath_atlas
pip install -e ".[dev]"
pytest -q
voc eval-stress-hard --out-dir runs/stress_hard
voc eval-patient-cohort --out-dir runs/patient_cohort_1000 --max-patients 12   # smoke
voc eval-vision --out-dir runs/vision_eval
```

## Branch hygiene

Active development lands on `main`. Historical feature branches that were
already squash/merged (or superseded by later PRs) are pruned from `origin`
so `main` remains the single source of truth for code **and** published results.
