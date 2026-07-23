# Model improve + external-validated 100-disease suite

In-sample external concordance uses open corpora that also inform prior fusion (panels / public breath / HBDB proxies) and can look near-perfect. Held-out literature_benchmarks concordance is the less circular check. HBDB SQL / Owlstone still blocked until user dumps arrive. Calibrator *.joblib untouched.

**Held-out literature_benchmarks concordance: 100.0%** (12 diseases)
**In-sample external concordance: 99.67%** (25 diseases with external GT)
**Literature panel concordance: 99.75%** (19 diseases)
**Disease connection suite: 100 diseases**
**Audit directional accuracy: 1.0** (min-fold pass 1.0)

## Fuse summary

- Diseases updated: 20
- VOC prior updates: 77
- Evidence edges: 195
- Europe PMC DOIs (provenance): 97
- Expanded MW studies cataloged: 25
- VOC identity enrichments: 0

## External validation (vs open data)

- Held-out lit-bench mean concordance: 100.0%
- In-sample mean concordance: 99.67%
- Diseases with external GT (in-sample): 25
- VOC checks: 175 (hits 174)

Top disease concordances (in-sample):

- `alzheimer_disease`: 100.0% (4/4)
- `ards`: 100.0% (9/9)
- `asthma`: 100.0% (19/19)
- `breast_invasive_carcinoma`: 100.0% (2/2)
- `cancer_prostate`: 100.0% (6/6)
- `cancer_stomach`: 100.0% (9/9)
- `chronic_kidney_disease`: 100.0% (4/4)
- `chronic_liver_disease`: 100.0% (5/5)
- `colon_adenocarcinoma`: 100.0% (3/3)
- `copd`: 100.0% (20/20)
- `covid19`: 100.0% (8/8)
- `epilepsy`: 100.0% (3/3)
- `glioblastoma`: 100.0% (3/3)
- `head_neck_cancer`: 100.0% (8/8)
- `heart_failure`: 100.0% (5/5)

## Disease-100 connections (sample)

- copd ↔ lung_adenocarcinoma (cos=0.792)
- chronic_bronchitis ↔ copd (cos=0.933)
- gut_dysbiosis ↔ inflammatory_bowel_disease (cos=0.522)
- major_depressive_disorder ↔ schizophrenia (cos=0.915)
- malaria ↔ sepsis (cos=0.166)
- alzheimer_disease ↔ parkinson_disease (cos=0.994)
- asthma ↔ copd (cos=0.873)
- inflammatory_bowel_disease ↔ sibo (cos=0.524)

## Artifacts

- `data/knowledge/EXTERNAL_EVIDENCE_FUSE.json`
- `data/knowledge/MODEL_EXTERNAL_100DISEASE.json`
- `data/knowledge/MODEL_EXTERNAL_100DISEASE.md`
- `data/knowledge/lit_compare/` (connection matrices / figures)

Generated: 2026-07-23T06:16:40.884031+00:00
