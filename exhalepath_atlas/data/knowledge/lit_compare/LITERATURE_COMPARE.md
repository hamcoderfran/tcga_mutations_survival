# Literature compare · demographics PCA · 100-disease connections

## Does the 1000-profile atlas match literature?

**Mean directional concordance: 100.0%** across 19 diseases with curated elevate/suppress panels (priority10 + review-backed expectations for COPD / bronchitis / lung cancer / asthma / T2D / schizophrenia / HF / IBD / malaria).

### What matches well

- **Smoking → benzene/toluene (and ethylbenzene)**: classic exogenous BTEX breath signal ([Owlstone / smoking VOC reviews](https://www.owlstonemedical.com/about/blog/2023/jul/07/origins-of-vocs/); Filipiak et al. smoking VOC lists).
- **COPD → hexanal / ethane / pentane**: oxidative lipid-peroxidation alkanes & aldehydes (PMC8405872; PMC7796324).
- **Lung cancer → hexanal / heptanal / nonanal / 2-butanone**: aldehyde–ketone pattern (JTO breath VOC reviews).
- **T2D → acetone**: ketone-body breath literature; atlas min-fold gates.
- **Schizophrenia → ↓ acetone / isoprene, ↑ pentane/ethane**: Magdeburg PTR-MS panel.
- **COPD ↔ chronic bronchitis closer than either ↔ LUAD**: expected obstructive continuum.

### Where it is only partly aligned / cautious

- Chronic bronchitis has **few dedicated breath GC-MS panels**; we treat it as an inflammatory airway neighbor of COPD (feasible, but not gold-standard validated).
- Individual VOC markers are **non-specific** across airway diseases — literature emphasizes multi-VOC patterns, which is why cosine neighborhoods matter more than single-marker claims (PMC7796324).
- Hybrid physiology previously **attenuated** smoking BTEX (bug; now fixed).

### Per-disease concordance

- **ards**: 100.0% (9/9)
- **asthma**: 100.0% (19/19)
- **cancer_prostate**: 100.0% (7/7)
- **cancer_stomach**: 100.0% (9/9)
- **chronic_bronchitis**: 100.0% (5/5)
- **copd**: 100.0% (21/21)
- **covid19**: 100.0% (8/8)
- **cystic_fibrosis**: 100.0% (11/11)
- **head_neck_cancer**: 100.0% (8/8)
- **heart_failure**: 100.0% (7/7)
- **inflammatory_bowel_disease**: 100.0% (2/2)
- **lung_adenocarcinoma**: 100.0% (7/7)
- **lung_squamous_cell_carcinoma**: 100.0% (6/6)
- **malaria**: 100.0% (7/7)
- **pneumonia_bacterial**: 100.0% (11/11)
- **schizophrenia**: 100.0% (5/5)
- **sleep_apnea**: 100.0% (9/9)
- **tuberculosis**: 100.0% (9/9)
- **type_2_diabetes**: 100.0% (3/3)

## Demographics / identity PCA — do we see clusters?

- n=100 · demo PC var=[0.20497991485569836, 0.17494992897026504] · VOC PC var=[0.29029514839029064, 0.21153170237741162]
- Silhouette demo←smoking: **0.27728892995946436**
- Silhouette demo←disease category: **None** (should be weak if demographics are not disease-leaking)
- Silhouette VOC←smoking: **0.039513092777246656**
- Silhouette VOC←category: **None**

Demographics PCA should cluster by smoking/sex (identity axes), not by disease. VOC PCA should show both disease-category structure and a smoking axis (benzene/toluene literature).

Figures: `figures/demo_pca_by_smoking.png`, `demo_pca_by_sex.png`, `demo_pca_by_disease_category.png`, `voc_pca_by_smoking.png`, `voc_pca_by_age.png`.

## Bugs unearthed & fixed

- **ppb_ceiling_erases_exogenous**: ppb upper clip now includes headroom for smoking/age multipliers so saturated priors cannot nullify BTEX/age effects
- **smoking_hybrid_attenuation**: Smoking BTEX boost applied after physio/legacy blend (benzene/toluene/pentane/ethylbenzene)
- **demographics_ignored_mechanistically**: Mild age (oxidative VOCs) and sex modulators applied post-blend with clip headroom
- **lusc_copd_collapse**: LUSC priors shifted toward cancer aldehydes/ketones; ethane demoted to separate from COPD
- **cohort_nogene_error_marks_fail**: Gene-shift secondary failure no longer marks primary patient prediction as not-ok

### Smoking Δlog2fc after fix (current − never)

- copd: benzene=+0.450, toluene=+0.450, pentane=+0.450, ethylbenzene=+0.450, hexanal=+0.000
- lung_adenocarcinoma: benzene=+0.450, toluene=+0.450, pentane=+0.450, ethylbenzene=+0.450, hexanal=+0.000

### Focus centroid L2 after LUSC fix

- `chronic_bronchitis__vs__lung_adenocarcinoma`: 3.969
- `chronic_bronchitis__vs__lung_squamous_cell_carcinoma`: 2.476
- `copd__vs__chronic_bronchitis`: 2.909
- `copd__vs__lung_adenocarcinoma`: 4.105
- `copd__vs__lung_squamous_cell_carcinoma`: 3.242
- `lung_adenocarcinoma__vs__lung_squamous_cell_carcinoma`: 2.646

## 100-disease suite — interesting connections

Ran **100** distinct atlas diseases on a shared patient template (55y female, never-smoker, default site, driver genes when available).

### Top cross-category VOC neighbors

- `multiple_sclerosis` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.995
- `hypertension` (cardiovascular) ↔ `sickle_cell` (unspecified): cosine=0.992
- `hiv` (infectious) ↔ `hypertension` (cardiovascular): cosine=0.991
- `hypertension` (cardiovascular) ↔ `influenza` (infectious): cosine=0.991
- `ards` (pulmonary) ↔ `creutzfeldt_jakob` (neurological): cosine=0.991
- `migraine` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.991
- `anxiety` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.991
- `adhd` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.991
- `rheumatoid_arthritis` (inflammatory) ↔ `stroke` (neurological): cosine=0.991
- `ptsd` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.991
- `huntington` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.991
- `bipolar` (neurological) ↔ `rheumatoid_arthritis` (inflammatory): cosine=0.990

### Highlighted biological bridges

- **copd ↔ lung_adenocarcinoma** (cos=0.791, shared=5): Oxidative aldehydes overlap; cancer should still separate on ketone/aldehyde mix
  - VOCs: `heptanal|hexanal|nonanal|pentane|propionaldehyde`
- **chronic_bronchitis ↔ copd** (cos=0.931): Obstructive airway inflammation continuum (literature-expected)
- **gut_dysbiosis ↔ inflammatory_bowel_disease** (cos=0.528, shared=13): Microbial sulfur / putrefaction VOCs
  - VOCs: `acetaldehyde|ammonia|carbon_disulfide|dimethyl_amine|dimethyl_disulfide|ethanol|ethyl_acetate|hydrogen_sulfide|indole|methyl_mercaptan|phenol|propanol`
- **major_depressive_disorder ↔ schizophrenia** (cos=0.915): Neuro-oxidative + microbiome-adjacent breath features
- **malaria ↔ sepsis** (cos=0.165): Systemic oxidative / infectious breath stress
- **alzheimer_disease ↔ parkinson_disease** (cos=0.991): Neurodegeneration oxidative alkanes/aldehydes
- **asthma ↔ copd** (cos=0.874): Airway oxidative alkane overlap with distinct ketone/ester accents
- **chronic_liver_disease ↔ hepatocellular_carcinoma** (cos=0.757, shared=5): Hepatic sulfur / ammonia axis toward malignancy
  - VOCs: `2_butanone|ammonia|dms|hydrogen_sulfide|limonene`
- **inflammatory_bowel_disease ↔ sibo** (cos=0.530, shared=11): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|ammonia|dimethyl_amine|ethanol|ethyl_acetate|hydrogen_sulfide|indole|methyl_mercaptan|phenol|propanol|trimethylamine`
- **clostridioides_difficile_infection ↔ inflammatory_bowel_disease** (cos=0.511, shared=11): High shared elevated-VOC bridge across categories
  - VOCs: `ammonia|dimethyl_amine|dimethyl_disulfide|ethanol|ethyl_acetate|hydrogen_sulfide|indole|methyl_mercaptan|phenol|propanol|trimethylamine`
- **ards ↔ lung_squamous_cell_carcinoma** (cos=0.948, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **covid19 ↔ lung_squamous_cell_carcinoma** (cos=0.945, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ colon_adenocarcinoma** (cos=0.729, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **colon_adenocarcinoma ↔ covid19** (cos=0.728, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **covid19 ↔ pancreatic_adenocarcinoma** (cos=0.689, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ pancreatic_adenocarcinoma** (cos=0.686, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ glioblastoma** (cos=0.642, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`

## Outputs

- `LITERATURE_COMPARE.md` / `lit_demo_disease100.json`
- `demographics_pca_coords.csv`
- `disease100_voc_matrix.csv` / `disease100_cross_category_neighbors.csv` / `disease100_voc_bridges.csv`
- `figures/*`
