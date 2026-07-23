# ExhalePath results index

Single entry point for past evaluations and reports on `main`.
Trained calibrators (`data/models/*.joblib`) are **not** modified by eval merges.

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
| VOC coverage audit (open corpus ≥99%) | [`COVERAGE_AUDIT.md`](COVERAGE_AUDIT.md) | `voc eval-coverage` |
| Extended VOC catalog | [`voc_extended_catalog.json`](voc_extended_catalog.json) | (from `eval-coverage`) |
| Integrity / anti-poisoning | [`../datasources/INTEGRITY_MANIFEST.json`](../datasources/INTEGRITY_MANIFEST.json) | SHA-256 of secured artifacts |

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
- `runs/coverage_audit/` · `runs/lit_compare/`

## Quick regression commands

```bash
cd exhalepath_atlas
pip install -e ".[dev]"
pytest -q
voc eval-coverage --offline
voc eval-stress-hard --out-dir runs/stress_hard
voc eval-patient-cohort --out-dir runs/patient_cohort_1000 --max-patients 12   # smoke
voc eval-vision --out-dir runs/vision_eval
```

## Branch hygiene

Active development lands on `main`. Historical feature branches that were
already squash/merged (or superseded by later PRs) are pruned from `origin`
so `main` remains the single source of truth for code **and** published results.
