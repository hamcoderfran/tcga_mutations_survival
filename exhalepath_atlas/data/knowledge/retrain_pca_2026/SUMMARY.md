# Retrain + multi-disease PCA retest

## Calibrator (retrained)
- VOC models: 16
- mean MAE log2fc: 0.1919
- global fallback MAE: 0.1438323790678971

## Retests
| Suite | Result |
|-------|--------|
| priority10 directional | **100.0%** (13 diseases) |
| public breath | elev recall@15=0.748, dir=0.982 |
| vision fidelity | **98.25%** (50 diseases) |
| stress-hard | **100%** (50/50) |
| patient cohort | **100%** ok (1000/1000 patients) |

## Disease relatability (18 diseases with n≥4, VOC centroid L2)
- silhouette: **0.204**
- closest: `adhd` ↔ `als` (1.842)
- farthest: `copd` ↔ `asthma` (11.979)

Pulmonary focus: COPD↔bronchitis closer than either↔LUAD (obstructive vs malignancy).

## Figures
Under `data/knowledge/cohort_1000/figures/`:
- `voc_pca_multidisease.png`
- `voc_pca.png` / `voc_umap.png`
- `voc_pca_focus_pulmonary.png`
