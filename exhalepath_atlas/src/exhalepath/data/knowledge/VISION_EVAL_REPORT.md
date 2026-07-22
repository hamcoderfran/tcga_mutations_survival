# ExhalePath vision evaluation (50 diseases)

**Vision fidelity (all 50): 88.25%**
**Evidence-backed (Grade A+B): 89.66%**
**Held-out style (Grade A+B+C): 90.0%**
**Grade D prior consistency: 85.14%**

## Metric means (all profiles)

- direction: 96.5%
- fold: 100.0%
- topk: 98.0%
- pathway: 80.3%
- cell: 59.0%
- site: 88.3%

## By evidence grade

- **Grade A** (n=6): composite 92.58% (dir=98%, fold=100%, topk=100%, pathway=96%)
- **Grade B** (n=20): composite 88.78% (dir=98%, fold=100%, topk=98%, pathway=72%)
- **Grade C** (n=6): composite 91.49% (dir=100%, fold=100%, topk=100%, pathway=81%)
- **Grade D** (n=18): composite 85.14% (dir=93%, fold=—, topk=96%, pathway=84%)

## Per-disease composite

| ID | Grade | Composite % | Direction | Fold | Top-k | Pathway | Cell | Site |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| asthma | A | 92.75 | 100% | 100% | 100% | 100% | 33% | 100% |
| copd | A | 89.13 | 100% | 100% | 100% | 100% | 0% | 100% |
| heart_failure | A | 92.3 | 100% | 100% | 100% | 75% | 67% | 100% |
| lung_adenocarcinoma | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| malaria | A | 81.27 | 86% | 100% | 100% | 100% | 0% | 50% |
| type_2_diabetes | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| alzheimer_disease | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| ards | B | 80.9 | 100% | — | 100% | 75% | 0% | 100% |
| breast_invasive_carcinoma | B | 93.48 | 100% | 100% | 67% | 100% | 100% | 100% |
| cancer_prostate | B | 76.85 | 100% | — | 100% | 33% | 33% | 100% |
| cancer_stomach | B | 79.75 | 100% | — | 100% | 25% | 67% | 100% |
| colon_adenocarcinoma | B | 89.13 | 100% | 100% | 100% | 67% | 50% | 100% |
| covid19 | B | 80.32 | 100% | — | 100% | 50% | 33% | 100% |
| cystic_fibrosis | B | 92.75 | 100% | 100% | 100% | 100% | 33% | 100% |
| glioblastoma | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| head_neck_cancer | B | 75.63 | 100% | 100% | 100% | 25% | 33% | 50% |
| hepatocellular_carcinoma | B | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| major_depressive_disorder | B | 91.54 | 95% | — | 100% | 67% | 100% | 100% |
| ovarian_cancer | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| pancreatic_adenocarcinoma | B | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| parkinson_disease | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| pneumonia_bacterial | B | 80.98 | 100% | 100% | 100% | 50% | 0% | 100% |
| schizophrenia | B | 85.0 | 71% | — | 100% | 67% | 100% | 100% |
| sleep_apnea | B | 74.07 | 100% | — | 100% | 50% | 33% | 50% |
| tuberculosis | B | 86.11 | 100% | — | 100% | 100% | 0% | 100% |
| type1_diabetes | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| chronic_kidney_disease | C | 86.41 | 100% | 100% | 100% | 50% | 50% | 100% |
| chronic_liver_disease | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| creutzfeldt_jakob | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| gut_dysbiosis | C | 89.13 | 100% | 100% | 100% | 67% | 50% | 100% |
| inflammatory_bowel_disease | C | 78.8 | 100% | 100% | 100% | 67% | 0% | 50% |
| sibo | C | 94.57 | 100% | 100% | 100% | 100% | 50% | 100% |
| autism_spectrum_disorder | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| endometriosis | D | 76.39 | 100% | — | 100% | 50% | 50% | 50% |
| epilepsy | D | 93.06 | 100% | — | 100% | 67% | 100% | 100% |
| esophageal_cancer | D | 84.72 | 100% | — | 67% | 100% | 50% | 100% |
| heart_disease | D | 79.17 | 100% | — | 100% | 67% | 0% | 100% |
| helicobacter_pylori_infection | D | 93.06 | 100% | — | 100% | 100% | 50% | 100% |
| hiv | D | 75.69 | 85% | — | 100% | 100% | 0% | 50% |
| hypertension | D | 76.39 | 100% | — | 100% | 50% | 50% | 50% |
| influenza | D | 81.94 | 85% | — | 100% | 100% | 0% | 100% |
| multiple_sclerosis | D | 93.06 | 100% | — | 100% | 100% | 50% | 100% |
| nafld | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| obesity | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| pcos | D | 95.83 | 100% | — | 100% | 100% | 100% | 67% |
| rheumatoid_arthritis | D | 62.27 | 75% | — | 67% | 67% | 33% | 50% |
| sepsis | D | 80.92 | 87% | — | 100% | 100% | 33% | 50% |
| sickle_cell | D | 76.39 | 100% | — | 100% | 50% | 50% | 50% |
| thyroid | D | 63.66 | 42% | — | 100% | 67% | 50% | 50% |
| traumatic_brain_injury | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |

## Notes

- Grade **A**: quantified measured fold changes.
- Grade **B**: strong directional literature / priority panels.
- Grade **C**: held-out / zero-shot literature expectations.
- Grade **D**: atlas-prior–aligned approximate GT (consistency, not independent accuracy).

vision_fidelity_pct is the headline % closeness to the original vision (disease+location+patient → VOC ranking + pathway/cell/site explainability) across all 50 profiles. Prefer evidence_backed_pct for independent truthfulness; Grade D inflates fidelity via prior circularity.
