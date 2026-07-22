# ExhalePath real-data pipeline evaluation

Generated: `2026-07-22T14:57:49.244528+00:00`

## Pipeline status

**Overall: PASSED**

### Corpus (real labels only)

- Type: `real_breath_only` · synthetic=False
- VOC training targets: **189** across **30** cases / **13** diseases
- Diseases: ards, asthma, cancer_prostate, cancer_stomach, copd, covid19, cystic_fibrosis, head_neck_cancer, heart_failure, malaria, pneumonia_bacterial, sleep_apnea, tuberculosis
- Sources: metabolomics_workbench:ST000587,ST000883,ST003200, scientific_data_figshare:23522490, literature_panels:priority10_voc_panels.json
- Label origins: `{"literature:mixed": 64, "literature:directional_only": 53, "measured_cohort_relative": 35, "measured_cohort": 19, "literature:quantified": 18}`

### Calibrator

- VOC models: **16** · mean MAE(log2fc)=0.1919

## Real-data evaluation scores

| Suite | Headline score | Detail |
|---|---:|---|
| Priority-10 literature panels | **100.0%** | 13 diseases · directional vs measured_log2fc |
| Public breath (all 7 cases) | **98.2%** | elev recall@15=74.8% |
| Public breath — literature | **100.0%** | elev recall@15=100.0% |
| Public breath — Sci Data 2024 | elev dir **88.9%** | elev recall@15=41.3% |
| Comorbidity clinical | **83.5%** | elev dir=100.0% · elev@k=95.0% · Magdeburg SZ + ST003181 |
| Multisite literature (20×5) | **98.8%** | dir=100.0% · site=93.8% · A — strong literature agreement for directional VOC modeling |
| Literature audit | **100.0%** | directional=100.0% · zero-shot=100.0% |
| Vision suite (50 diseases) | **98.02%** | evidence-backed=97.28% · held-out-style=97.37% |
| Completion CI gate | **PASSED** | all 8 gates green |
| Diabetes clinical profile | **PASSED** | acetone ≥1.8× gate (PMID:21903721) |
| CJD clinical profile | ran | incubating subtle → clinical obvious (oxidative tempo) |

## Completion gates

```json
{
  "pytest": true,
  "audit": true,
  "multisite_composite_ge_0.90": true,
  "lit_public_directional_ge_0.90": true,
  "lit_public_recall_ge_0.85": true,
  "sci_data_elev_directional_ge_0.60": true,
  "diseases_ge_100": true,
  "datasources_ge_12": true
}
```

## Honesty notes

- Priority-10 directional accuracy can look near-perfect when panels also inform atlas priors / training labels (partial circularity).
- Sci Data 2024 recall@k is lower because source tables lack healthy controls (cross-cohort differentials only).
- Vision Grade D is prior-consistency; prefer evidence-backed % for independent truthfulness.
- Research tool — not a medical device.

## Reproduce

```bash
cd exhalepath_atlas && pip install -e ".[dev]"
voc build-real-corpus && voc train
voc eval-priority10
voc eval-public-breath
voc eval-comorbidity-clinical
voc eval-multisite && voc audit
voc eval-vision
voc diabetes-profile --age 55 --sex male
voc cjd-profile --age 62 --sex female
voc eval-completion
# or: python scripts/run_real_pipeline_eval.py
```
