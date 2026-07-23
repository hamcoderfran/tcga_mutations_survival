# Model improve + external-validated 100-disease suite

Held-out literature_benchmarks concordance is the primary non-circular metric. In-sample external concordance overlaps fuse sources (panels/public/HBDB) and can look near-perfect by construction — treat as coverage of open GT, not independent accuracy. Fuse always restarts from disease_voc_priors.pristine.json (idempotent). HBDB SQL / Owlstone still blocked until user dumps arrive. Calibrator *.joblib untouched.

**Held-out literature_benchmarks concordance: 100.0%** (13 diseases) ← primary non-circular metric
**In-sample external concordance: 99.69%** (27 diseases; circular=True)
**Literature panel concordance: 100.0%** (19 diseases)
**Disease connection suite: 100 diseases**
**Audit directional accuracy: 1.0** (min-fold pass 1.0)

## Fuse summary

- Diseases updated: 20
- VOC prior updates: 76
- Histology guard edits: 7
- Evidence edges: 195
- Europe PMC DOIs (provenance): 97
- Expanded MW studies cataloged: 25
- Blend from pristine (idempotent): True

## External validation (vs open data)

- Held-out lit-bench mean concordance: 100.0%
- In-sample mean concordance: 99.69% (includes fused sources — interpret cautiously)
- Diseases with external GT (in-sample): 27
- VOC checks: 189 (hits 188)

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

- copd ↔ lung_adenocarcinoma (cos=0.791)
- chronic_bronchitis ↔ copd (cos=0.931)
- gut_dysbiosis ↔ inflammatory_bowel_disease (cos=0.528)
- major_depressive_disorder ↔ schizophrenia (cos=0.915)
- malaria ↔ sepsis (cos=0.165)
- alzheimer_disease ↔ parkinson_disease (cos=0.991)
- asthma ↔ copd (cos=0.874)
- chronic_liver_disease ↔ hepatocellular_carcinoma (cos=0.757)

## Artifacts

- `data/knowledge/disease_voc_priors.pristine.json`
- `data/knowledge/EXTERNAL_EVIDENCE_FUSE.json`
- `data/knowledge/MODEL_EXTERNAL_100DISEASE.json`
- `data/knowledge/MODEL_EXTERNAL_100DISEASE.md`
- `data/knowledge/lit_compare/`

Generated: 2026-07-23T17:05:58.361567+00:00
