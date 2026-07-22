# ExhalePath vision evaluation (50 diseases)

**Vision fidelity (all 50): 98.01%**
**Evidence-backed (Grade A+B): 97.26%**
**Held-out style (Grade A+B+C): 97.35%**
**Grade D prior consistency: 99.18%**

## Metric means (all profiles)

- direction: 98.9%
- fold: 100.0%
- topk: 99.3%
- pathway: 91.8%
- cell: 100.0%
- site: 99.3%

## By evidence grade

- **Grade A** (n=6): composite 99.8% (dir=99%, fold=100%, topk=100%, pathway=100%)
- **Grade B** (n=20): composite 96.5% (dir=98%, fold=100%, topk=98%, pathway=85%)
- **Grade C** (n=6): composite 97.74% (dir=100%, fold=100%, topk=100%, pathway=86%)
- **Grade D** (n=18): composite 99.18% (dir=99%, fold=—, topk=100%, pathway=98%)

## Per-disease composite

| ID | Grade | Composite % | Direction | Fold | Top-k | Pathway | Cell | Site |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| asthma | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| copd | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| heart_failure | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| lung_adenocarcinoma | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| malaria | A | 98.81 | 95% | 100% | 100% | 100% | 100% | 100% |
| type_2_diabetes | A | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| alzheimer_disease | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| ards | B | 94.79 | 100% | — | 100% | 75% | 100% | 100% |
| breast_invasive_carcinoma | B | 93.48 | 100% | 100% | 67% | 100% | 100% | 100% |
| cancer_prostate | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| cancer_stomach | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| colon_adenocarcinoma | B | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| covid19 | B | 89.58 | 100% | — | 100% | 50% | 100% | 100% |
| cystic_fibrosis | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| glioblastoma | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| head_neck_cancer | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| hepatocellular_carcinoma | B | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| major_depressive_disorder | B | 91.54 | 95% | — | 100% | 67% | 100% | 100% |
| ovarian_cancer | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| pancreatic_adenocarcinoma | B | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| parkinson_disease | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| pneumonia_bacterial | B | 91.85 | 100% | 100% | 100% | 50% | 100% | 100% |
| schizophrenia | B | 85.0 | 71% | — | 100% | 67% | 100% | 100% |
| sleep_apnea | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| tuberculosis | B | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| type1_diabetes | B | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| chronic_kidney_disease | C | 91.85 | 100% | 100% | 100% | 50% | 100% | 100% |
| chronic_liver_disease | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| creutzfeldt_jakob | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| gut_dysbiosis | C | 94.57 | 100% | 100% | 100% | 67% | 100% | 100% |
| inflammatory_bowel_disease | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| sibo | C | 100.0 | 100% | 100% | 100% | 100% | 100% | 100% |
| autism_spectrum_disorder | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| endometriosis | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| epilepsy | D | 93.06 | 100% | — | 100% | 67% | 100% | 100% |
| esophageal_cancer | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| heart_disease | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| helicobacter_pylori_infection | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| hiv | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| hypertension | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| influenza | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| multiple_sclerosis | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| nafld | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| obesity | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| pcos | D | 95.83 | 100% | — | 100% | 100% | 100% | 67% |
| rheumatoid_arthritis | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| sepsis | D | 96.43 | 87% | — | 100% | 100% | 100% | 100% |
| sickle_cell | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| thyroid | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |
| traumatic_brain_injury | D | 100.0 | 100% | — | 100% | 100% | 100% | 100% |

## Notes

- Grade **A**: quantified measured fold changes.
- Grade **B**: strong directional literature / priority panels.
- Grade **C**: held-out / zero-shot literature expectations.
- Grade **D**: atlas-prior–aligned approximate GT (consistency, not independent accuracy).

vision_fidelity_pct is the headline % closeness to the original vision (disease+location+patient → VOC ranking + pathway/cell/site explainability) across all 50 profiles. Prefer evidence_backed_pct for independent truthfulness; Grade D inflates fidelity via prior circularity.
