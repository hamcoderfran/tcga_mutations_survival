# ExhalePath adversarial stress suite (50 hard cases)

**Pass rate: 88.0%** (44/50)

## By category

- **robustness**: 4/7 (57%)
- **gene_trap**: 3/4 (75%)
- **alias**: 6/6 (100%)
- **location**: 5/5 (100%)
- **clinical_hard**: 7/9 (78%)
- **zero_shot**: 4/4 (100%)
- **conflict**: 9/9 (100%)
- **metastatic**: 6/6 (100%)

## Failure modes

- `wrong_or_missing_suppress`: 3
- `no_exception`: 2
- `location_unmatched`: 1
- `unresolved_or_not_breast`: 1
- `near_healthy`: 1

## Failed cases

- **symbol_location** (robustness): location_unmatched
- **comorbidity_weight_overflow** (robustness): no_exception
  - error: `ValidationError: 1 validation error for DiseaseQuery
comorbidity_weight
  Input should be less than or equal to 1.5 [type=less_than_equal, input_value=3.5, input_type=float]
    For further information visit https://errors.pydantic.dev/2.13/v/less_than_equal`
- **comorbidity_weight_negative** (robustness): no_exception
  - error: `ValidationError: 1 validation error for DiseaseQuery
comorbidity_weight
  Input should be greater than or equal to 0 [type=greater_than_equal, input_value=-1.0, input_type=float]
    For further information visit https://errors.pydantic.dev/2.13/v/greater_than_equal`
- **gene_BRCA1_as_disease** (gene_trap): unresolved_or_not_breast, near_healthy
- **sz_suppress_acetone** (clinical_hard): suppress:acetone
- **sz_heart_suppress_tma** (clinical_hard): suppress:acetone, suppress:trimethylamine

## All cases

| ID | Cat | Pass | Disease | max\|log2\| |
|---|---|---|---|---:|
| empty_disease | robustness | ✓ | custom::unknown | 0.167 |
| whitespace_disease | robustness | ✓ | custom::unknown | 0.167 |
| nonsense_disease | robustness | ✓ | custom::zxqy_not_a_real_disease_999 | 0.167 |
| symbol_location | robustness | ✗ | asthma | 4.637 |
| comorbidity_weight_overflow | robustness | ✗ |  | 0.000 |
| comorbidity_weight_negative | robustness | ✗ |  | 0.000 |
| twenty_comorbidities | robustness | ✓ | major_depressive_disorder | 5.068 |
| gene_BRCA1_as_disease | gene_trap | ✗ | breast_invasive_carcinoma | 3.355 |
| gene_KRAS_as_disease | gene_trap | ✓ | custom::kras | 0.167 |
| gene_TP53_as_disease | gene_trap | ✓ | custom::tp53 | 0.167 |
| luad_with_driver_genes | gene_trap | ✓ | lung_adenocarcinoma | 4.622 |
| bare_pneumonia | alias | ✓ | pneumonia_bacterial | 4.819 |
| bare_diabetes | alias | ✓ | type_2_diabetes | 3.819 |
| heart_failure_not_cad | alias | ✓ | heart_failure | 4.331 |
| nafld_not_cirrhosis | alias | ✓ | nafld | 2.370 |
| autism_not_heart_failure | alias | ✓ | autism_spectrum_disorder | 2.620 |
| bare_cancer_unresolved | alias | ✓ | custom::cancer | 0.167 |
| head_neck_exact | location | ✓ | head_neck_cancer | 4.123 |
| pharynx_osa | location | ✓ | sleep_apnea | 4.311 |
| joint_ra | location | ✓ | rheumatoid_arthritis | 4.594 |
| left_breast_free_text | location | ✓ | breast_invasive_carcinoma | 3.355 |
| retroperitoneal_mass | location | ✓ | ovarian_cancer | 2.834 |
| sz_suppress_acetone | clinical_hard | ✗ | schizophrenia | 3.451 |
| sz_heart_suppress_tma | clinical_hard | ✗ | schizophrenia | 4.198 |
| mdd_obesity_ketones | clinical_hard | ✓ | major_depressive_disorder | 3.061 |
| t2d_acetone_gate | clinical_hard | ✓ | type_2_diabetes | 3.819 |
| copd_ethane_hexanal | clinical_hard | ✓ | copd | 4.837 |
| malaria_benzene_acetone | clinical_hard | ✓ | malaria | 3.311 |
| cirrhosis_sulfur | clinical_hard | ✓ | chronic_liver_disease | 3.490 |
| ibd_h2s | clinical_hard | ✓ | inflammatory_bowel_disease | 3.834 |
| sibo_fermentation | clinical_hard | ✓ | sibo | 5.082 |
| zero_shot_freckling | zero_shot | ✓ | custom::benign_essential_freckling_syndrome_zx_0 | 0.167 |
| zero_shot_madeup_syndrome | zero_shot | ✓ | custom::quigley_hart_metabolic_vapor_syndrome | 0.167 |
| cjd_oxidative | zero_shot | ✓ | creutzfeldt_jakob | 4.837 |
| epilepsy_brain_energy | zero_shot | ✓ | epilepsy | 2.655 |
| obesity_alone_ketones | conflict | ✓ | obesity | 2.370 |
| paad_kras_tp53 | conflict | ✓ | pancreatic_adenocarcinoma | 5.059 |
| hcc_sulfur_not_only_ketone | conflict | ✓ | hepatocellular_carcinoma | 3.067 |
| covid_aldehydes | conflict | ✓ | covid19 | 4.339 |
| tb_oxidative | conflict | ✓ | tuberculosis | 4.375 |
| ards_ventilator | conflict | ✓ | ards | 4.837 |
| thyroid_isoprene | conflict | ✓ | thyroid | 2.530 |
| sickle_oxidative | conflict | ✓ | sickle_cell | 3.947 |
| parkinson_oxidative | conflict | ✓ | parkinson_disease | 4.543 |
| luad_brain_met | metastatic | ✓ | lung_adenocarcinoma | 4.500 |
| breast_male | metastatic | ✓ | breast_invasive_carcinoma | 3.355 |
| prostate_stage_high | metastatic | ✓ | cancer_prostate | 3.803 |
| glioblastoma_warburg | metastatic | ✓ | glioblastoma | 5.059 |
| sepsis_icu | metastatic | ✓ | sepsis | 3.782 |
| influenza_not_bacterial | metastatic | ✓ | influenza | 4.132 |
