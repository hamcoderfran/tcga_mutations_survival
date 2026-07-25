# Novel diseases (no VOC background)

**Passed: True**

- novel queries: 10
- unresolved (custom::): 10/10
- near-healthy (peak |log2fc| < 0.45): 10/10
- mean peak |log2fc|: 0.167
- mean confidence: 0.608
- acetone Δppb (novel mean): 13.42
- acetone Δppb (T2D): 6374.0965545601575
- separation (T2D − novel): 6360.67 ppb

## Gates

- ✓ `all_novel_unresolved`
- ✓ `most_near_healthy`
- ✓ `mean_peak_abs_log2fc_lt_0_6`
- ✓ `acetone_sep_vs_t2d_gt_50ppb`
- ✓ `known_not_unresolved`

## Novel cases (top VOCs)

- **completely novel ailment QZZ-42** → `custom::completely_novel_ailment_qzz_42` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **Xylophage syndrome type Zeta** → `custom::xylophage_syndrome_type_zeta` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **hyperblue mitochondrial spark disease** → `custom::hyperblue_mitochondrial_spark_disease` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **quantum itch disorder** → `custom::quantum_itch_disorder` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **fibroquartz encephalopathy** → `custom::fibroquartz_encephalopathy` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **neon teal cholangioflux** → `custom::neon_teal_cholangioflux` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **sporadic purple glomerulopathy XYZ** → `custom::sporadic_purple_glomerulopathy_xyz` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **astral cartilage liquefaction** → `custom::astral_cartilage_liquefaction` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **cryptic umbra pancreatitis variant 9** → `custom::cryptic_umbra_pancreatitis_variant_9` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)
- **nonexistent pathogen Omega-7 breath plague** → `custom::nonexistent_pathogen_omega_7_breath_plague` unresolved=True peak|log2fc|=0.167 near_healthy=True
  - top: acetone(+0.04), ammonia(+0.04), isoprene(-0.03), methanol(-0.01)

## Known contrast

- **type 2 diabetes** → `type_2_diabetes` peak|log2fc|=3.844 acetone_Δppb=6374.0965545601575
- **asthma** → `asthma` peak|log2fc|=4.837 acetone_Δppb=134.1931579595705
- **COPD** → `copd` peak|log2fc|=4.837 acetone_Δppb=156.03005185354914
- **lung adenocarcinoma** → `lung_adenocarcinoma` peak|log2fc|=4.141 acetone_Δppb=94.66280886427307
