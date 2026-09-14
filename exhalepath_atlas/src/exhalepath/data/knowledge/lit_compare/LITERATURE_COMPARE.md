# Literature compare · demographics PCA · 100-disease connections

## Does the 1000-profile atlas match literature?

**Mean directional concordance: 100.0%** across 21 diseases with curated elevate/suppress panels (priority10 + MH Magdeburg/Gbaoui/OralChroma panels + review-backed expectations for COPD / bronchitis / lung cancer / asthma / T2D / schizophrenia / MDD / bipolar / HF / IBD / malaria).

### What matches well

- **Smoking → benzene/toluene (and ethylbenzene)**: classic exogenous BTEX breath signal ([Owlstone / smoking VOC reviews](https://www.owlstonemedical.com/about/blog/2023/jul/07/origins-of-vocs/); Filipiak et al. smoking VOC lists).
- **COPD → hexanal / ethane / pentane**: oxidative lipid-peroxidation alkanes & aldehydes (PMC8405872; PMC7796324).
- **Lung cancer → hexanal / heptanal / nonanal / 2-butanone**: aldehyde–ketone pattern (JTO breath VOC reviews).
- **T2D → acetone**: ketone-body breath literature; atlas min-fold gates.
- **Schizophrenia → ↓ acetone / isoprene / trimethylamine, ↑ pentane/ethane/CS2**: Magdeburg PTR-MS (doi:10.1080/15622975.2022.2040052; doi:10.1503/jpn.220139) + Phillips pentane/CS2.
- **MDD → ↑ ethanol/acetaldehyde, ↓ SCFAs / isoprene / TMA**: Magdeburg + Gbaoui breathomics (doi:10.3389/fpsyt.2022.1061326).
- **Bipolar → ↑ methyl_mercaptan (CH3SH)**: OralChroma VSC (doi:10.3390/jcm14062025); H2S/DMS assayed but unreported — not invented.
- **COPD ↔ chronic bronchitis closer than either ↔ LUAD**: expected obstructive continuum.

### Where it is only partly aligned / cautious

- Chronic bronchitis has **few dedicated breath GC-MS panels**; we treat it as an inflammatory airway neighbor of COPD (feasible, but not gold-standard validated).
- Individual VOC markers are **non-specific** across airway diseases — literature emphasizes multi-VOC patterns, which is why cosine neighborhoods matter more than single-marker claims (PMC7796324).
- Hybrid physiology previously **attenuated** smoking BTEX (bug; now fixed).
- Magdeburg psych figshare `19181742` is a **DOCX supplement**, not a patient×VOC intensity matrix — use `voc eval-mental-health` (panel-masked / mechanism-backed) rather than AUROC claims.

### Per-disease concordance

- **ards**: 100.0% (9/9)
- **asthma**: 100.0% (19/19)
- **bipolar**: 100.0% (2/2)
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
- **major_depressive_disorder**: 100.0% (8/8)
- **malaria**: 100.0% (7/7)
- **pneumonia_bacterial**: 100.0% (11/11)
- **schizophrenia**: 100.0% (10/10)
- **sleep_apnea**: 100.0% (9/9)
- **tuberculosis**: 100.0% (9/9)
- **type_2_diabetes**: 100.0% (3/3)

## Demographics / identity PCA — do we see clusters?

- n=80 · demo PC var=[0.2117228446350147, 0.17526251509713534] · VOC PC var=[0.25710262993130706, 0.19655930136815994]
- Silhouette demo←smoking: **0.278448309108771**
- Silhouette demo←disease category: **None** (should be weak if demographics are not disease-leaking)
- Silhouette VOC←smoking: **0.06015205286682022**
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

- `chronic_bronchitis__vs__lung_adenocarcinoma`: 4.125
- `chronic_bronchitis__vs__lung_squamous_cell_carcinoma`: 2.356
- `copd__vs__chronic_bronchitis`: 2.990
- `copd__vs__lung_adenocarcinoma`: 4.315
- `copd__vs__lung_squamous_cell_carcinoma`: 2.856
- `lung_adenocarcinoma__vs__lung_squamous_cell_carcinoma`: 2.817

## 100-disease suite — interesting connections

Ran **40** distinct atlas diseases on a shared patient template (55y female, never-smoker, default site, driver genes when available).

### Top cross-category VOC neighbors

- `ards` (pulmonary) ↔ `creutzfeldt_jakob` (neurological): cosine=0.990
- `alzheimer_disease` (neurological) ↔ `arrhythmia` (cardiovascular): cosine=0.989
- `anxiety` (neurological) ↔ `arrhythmia` (cardiovascular): cosine=0.988
- `als` (neurological) ↔ `arrhythmia` (cardiovascular): cosine=0.988
- `ards` (pulmonary) ↔ `covid19` (infectious): cosine=0.987
- `covid19` (infectious) ↔ `creutzfeldt_jakob` (neurological): cosine=0.984
- `arrhythmia` (cardiovascular) ↔ `endometriosis` (unspecified): cosine=0.982
- `als` (neurological) ↔ `atopic_dermatitis` (inflammatory): cosine=0.980
- `als` (neurological) ↔ `celiac` (inflammatory): cosine=0.980
- `arrhythmia` (cardiovascular) ↔ `atopic_dermatitis` (inflammatory): cosine=0.980
- `arrhythmia` (cardiovascular) ↔ `celiac` (inflammatory): cosine=0.980
- `adhd` (neurological) ↔ `arrhythmia` (cardiovascular): cosine=0.978

### Highlighted biological bridges

- **chronic_bronchitis ↔ copd** (cos=0.927): Obstructive airway inflammation continuum (literature-expected)
- **asthma ↔ copd** (cos=0.869): Airway oxidative alkane overlap with distinct ketone/ester accents
- **ards ↔ colon_adenocarcinoma** (cos=0.736, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **colon_adenocarcinoma ↔ covid19** (cos=0.724, shared=8): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ covid19** (cos=0.987, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `acetone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ cancer_kidney** (cos=0.951, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **cancer_kidney ↔ covid19** (cos=0.944, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ cervical_cancer** (cos=0.940, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ cancer_bladder** (cos=0.933, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **cervical_cancer ↔ covid19** (cos=0.932, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **ards ↔ endometrial_cancer** (cos=0.928, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `acetaldehyde|ethanol|heptanal|hexanal|nonanal|octanal|pentane`
- **cancer_bladder ↔ covid19** (cos=0.925, shared=7): High shared elevated-VOC bridge across categories
  - VOCs: `2_butanone|ethanol|heptanal|hexanal|nonanal|octanal|pentane`

## Outputs

- `LITERATURE_COMPARE.md` / `lit_demo_disease100.json`
- `demographics_pca_coords.csv`
- `disease100_voc_matrix.csv` / `disease100_cross_category_neighbors.csv` / `disease100_voc_bridges.csv`
- `figures/*`
