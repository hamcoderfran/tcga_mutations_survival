# Type 2 diabetes — clinical VOC profile

**Patient:** 55yo male, genes=TCF7L2,PPARG, location=systemic.

> Research hypothesis only — not a diagnostic or glycemic monitor. Breath acetone tracks ketosis / metabolic stress, not HbA1c directly. Absolute ppb are model estimates.

## Bottom line

Prediabetes acetone fold=3.05× (obvious). First meaningful ketone VOC hint: **prediabetes**. Becomes obvious: **prediabetes**. Poorly controlled T2D literature gate: PASS (acetone 33.54× vs ≥1.8×; dir_acc=1.0). With obesity comorbidity acetone=33.54×; ketotic exacerbation acetone=33.54×. Breath acetone reflects ketosis / FAO stress — complementary to, not a replacement for, glucose/HbA1c.

- Overall literature gate: **PASSED**
- First meaningful hint: `prediabetes`
- First obvious: `prediabetes`

## Phase tempo

| Phase | Verdict | Acetone fold | Lit pass |
|---|---|---:|:---:|
| prediabetes | obvious | 3.05× | no |
| controlled_t2d | obvious | 11.48× | yes |
| poorly_controlled | obvious | 33.54× | yes |
| t2d_obesity | obvious | 33.54× | yes |
| ketotic | obvious | 33.54× | yes |
| type1_reference | obvious | 33.54× | yes |

## Phases

### Prediabetes / insulin resistance

*HbA1c ~5.7–6.4%; fasting glucose elevated; no ketosis*

Insulin resistance with mild fatty-acid oxidation stress; breath acetone may be subtly up vs healthy.

- Burden: stage I, frac=0.12, activity=1.0, prior_scale=0.35
- Comorbidities: —
- Ketone panel: **obvious** (acetone 3.05×, mean |log2FC| ketone=0.8196)
- Literature eval: FAIL (dir_acc=0.6666666666666666)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 1454.62 | +977.62 | 3.05x |
| 2 | Ammonia | 250.00 | 469.29 | +219.29 | 1.88x |
| 3 | Isoprene | 105.00 | 92.44 | -12.56 | 0.88x |
| 4 | n-Pentane | 8.50 | 20.03 | +11.53 | 2.36x |
| 5 | Methanol | 160.00 | 152.34 | -7.66 | 0.95x |
| 6 | Ethanol | 120.00 | 126.77 | +6.77 | 1.06x |
| 7 | 2-Propanol | 25.00 | 27.58 | +2.58 | 1.10x |
| 8 | 2-Butanone | 1.20 | 3.26 | +2.06 | 2.71x |
| 9 | Acetaldehyde | 15.00 | 14.31 | -0.69 | 0.95x |
| 10 | Methyl acetate | 2.00 | 1.73 | -0.27 | 0.87x |

Mechanisms:
- Acetone is predicted elevated in Type 2 diabetes mellitus; primarily via Ketone body metabolism (score=0.07; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Metabolically stressed neuron cells in brain (density×activity=0.121, Census fraction≈0.393); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 2 diabetes mellitus; primarily via Urea cycle / ammonia detoxification (score=0.00; genes CPS1, OTC, ASS1, ASL); sourced from Metabolically stressed neuron cells in brain (density×activity=0.121, Census fraction≈0.393); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Isoprene is predicted reduced in Type 2 diabetes mellitus; primarily via Mevalonate / cholesterol biosynthesis (score=0.00; genes HMGCR, MVK, MVD, FDPS); sourced from Metabolically stressed neuron cells in brain (density×activity=0.121, Census fraction≈0.393); biosynthetic chain(s): isoprene_mevalonate; key genetic/pathway nodes: HMGCR, MVK, MVD, PIK3CA, PTEN, AKT1; Census enriched population example: oligodendrocyte in brain; KEGG C16521 (Terpenoid backbone biosynthesis, Biosynthesis of terpenoids and steroids); EC 4.2.3.27; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

### Controlled T2D

*HbA1c ~6.5–7.5% on diet/metformin; euglycemic most days*

Established type 2 diabetes with reasonable control; ketone-body VOC axis mildly activated.

- Burden: stage II, frac=0.25, activity=1.25, prior_scale=0.65
- Comorbidities: —
- Ketone panel: **obvious** (acetone 11.48×, mean |log2FC| ketone=1.9316)
- Literature eval: PASS (dir_acc=1.0)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 5477.69 | +5000.69 | 11.48x |
| 2 | Ammonia | 250.00 | 1166.90 | +916.90 | 4.67x |
| 3 | Ethanol | 120.00 | 138.72 | +18.72 | 1.16x |
| 4 | 2-Butanone | 1.20 | 15.58 | +14.38 | 12.98x |
| 5 | n-Pentane | 8.50 | 20.21 | +11.71 | 2.38x |
| 6 | Isoprene | 105.00 | 95.00 | -10.00 | 0.90x |
| 7 | 2-Propanol | 25.00 | 31.04 | +6.04 | 1.24x |
| 8 | Methanol | 160.00 | 161.51 | +1.51 | 1.01x |
| 9 | Hexanal | 4.20 | 5.13 | +0.93 | 1.22x |
| 10 | Phenol | 3.50 | 3.84 | +0.34 | 1.10x |

Mechanisms:
- Acetone is predicted elevated in Type 2 diabetes mellitus; primarily via Ketone body metabolism (score=0.19; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.312); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 2 diabetes mellitus; primarily via Urea cycle / ammonia detoxification (score=0.00; genes CPS1, OTC, ASS1, ASL); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.312); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ethanol is predicted elevated in Type 2 diabetes mellitus; primarily via Glycolysis / Warburg metabolism (score=0.00; genes HK2, PKM, LDHA, PFKP); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.312); biosynthetic chain(s): ethanol_fermentation; key genetic/pathway nodes: HK2, PKM, LDHA, KRAS, NRAS, HRAS; Census enriched population example: oligodendrocyte in brain; KEGG C00469 (Glycolysis / Gluconeogenesis, Pyruvate metabolism); EC 1.1.1.1, 1.1.1.2, 1.1.1.71; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

### Poorly controlled T2D

*HbA1c ≥8.5%; hyperglycemia; possible mild ketonemia*

Chronic hyperglycemia with stronger ketone-body and FAO pathway bias; acetone / 2-butanone / isopropanol expected clearly elevated.

- Burden: stage III, frac=0.45, activity=1.7, prior_scale=1.0
- Comorbidities: —
- Ketone panel: **obvious** (acetone 33.54×, mean |log2FC| ketone=2.5326)
- Literature eval: PASS (dir_acc=1.0)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 16000.00 | +15523.00 | 33.54x |
| 2 | Ammonia | 250.00 | 2936.51 | +2686.51 | 11.75x |
| 3 | Ethanol | 120.00 | 151.51 | +31.51 | 1.26x |
| 4 | 2-Butanone | 1.20 | 20.00 | +18.80 | 16.67x |
| 5 | Methanol | 160.00 | 172.81 | +12.81 | 1.08x |
| 6 | n-Pentane | 8.50 | 20.45 | +11.95 | 2.41x |
| 7 | 2-Propanol | 25.00 | 36.03 | +11.03 | 1.44x |
| 8 | Isoprene | 105.00 | 95.34 | -9.66 | 0.91x |
| 9 | Hexanal | 4.20 | 5.18 | +0.98 | 1.23x |
| 10 | Phenol | 3.50 | 4.11 | +0.61 | 1.17x |

Mechanisms:
- Acetone is predicted elevated in Type 2 diabetes mellitus; primarily via Ketone body metabolism (score=0.45; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.765); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 2 diabetes mellitus; primarily via Urea cycle / ammonia detoxification (score=0.00; genes CPS1, OTC, ASS1, ASL); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.765); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ethanol is predicted elevated in Type 2 diabetes mellitus; primarily via Glycolysis / Warburg metabolism (score=0.00; genes HK2, PKM, LDHA, PFKP); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.765); biosynthetic chain(s): ethanol_fermentation; key genetic/pathway nodes: HK2, PKM, LDHA, KRAS, NRAS, HRAS; Census enriched population example: oligodendrocyte in brain; KEGG C00469 (Glycolysis / Gluconeogenesis, Pyruvate metabolism); EC 1.1.1.1, 1.1.1.2, 1.1.1.71; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

### T2D + obesity comorbidity

*Poorly controlled T2D with BMI ≥30; adipose + hepatic burden*

Comorbid obesity fuses adipose lipolysis / FAO priors onto the T2D ketone axis (common clinical phenotype).

- Burden: stage III, frac=0.5, activity=1.8, prior_scale=1.0
- Comorbidities: obesity
- Ketone panel: **obvious** (acetone 33.54×, mean |log2FC| ketone=2.6596)
- Literature eval: PASS (dir_acc=1.0)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 16000.00 | +15523.00 | 33.54x |
| 2 | Ammonia | 250.00 | 5020.36 | +4770.36 | 20.08x |
| 3 | Ethanol | 120.00 | 151.83 | +31.83 | 1.27x |
| 4 | 2-Propanol | 25.00 | 44.40 | +19.40 | 1.78x |
| 5 | 2-Butanone | 1.20 | 20.00 | +18.80 | 16.67x |
| 6 | Methanol | 160.00 | 173.93 | +13.93 | 1.09x |
| 7 | Isoprene | 105.00 | 95.97 | -9.03 | 0.91x |
| 8 | n-Pentane | 8.50 | 17.34 | +8.84 | 2.04x |
| 9 | 2-Pentanone | 1.40 | 2.25 | +0.85 | 1.61x |
| 10 | Phenol | 3.50 | 4.14 | +0.64 | 1.18x |

Mechanisms:
- Acetone is predicted elevated in Type 2 diabetes mellitus; primarily via Ketone body metabolism (score=0.46; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.900); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 2 diabetes mellitus; primarily via Urea cycle / ammonia detoxification (score=0.07; genes CPS1, OTC, ASS1, ASL); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.900); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ethanol is predicted elevated in Type 2 diabetes mellitus; primarily via Glycolysis / Warburg metabolism (score=0.00; genes HK2, PKM, LDHA, PFKP); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.900); biosynthetic chain(s): ethanol_fermentation; key genetic/pathway nodes: HK2, PKM, LDHA, KRAS, NRAS, HRAS; Census enriched population example: oligodendrocyte in brain; KEGG C00469 (Glycolysis / Gluconeogenesis, Pyruvate metabolism); EC 1.1.1.1, 1.1.1.2, 1.1.1.71; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

### Ketotic exacerbation / DKA-adjacent

*Marked ketosis (illness, SGLT2i, or T1D-like DKA)*

Strong ketone-body metabolism activation; breath acetone is the classic marker (literature folds often ≫2× healthy).

- Burden: stage IV, frac=0.7, activity=2.2, prior_scale=1.25
- Comorbidities: —
- Ketone panel: **obvious** (acetone 33.54×, mean |log2FC| ketone=2.6297)
- Literature eval: PASS (dir_acc=1.0)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 16000.00 | +15523.00 | 33.54x |
| 2 | Ammonia | 250.00 | 6123.12 | +5873.12 | 24.49x |
| 3 | Ethanol | 120.00 | 163.33 | +43.33 | 1.36x |
| 4 | Methanol | 160.00 | 185.86 | +25.86 | 1.16x |
| 5 | 2-Butanone | 1.20 | 20.00 | +18.80 | 16.67x |
| 6 | 2-Propanol | 25.00 | 41.76 | +16.76 | 1.67x |
| 7 | n-Pentane | 8.50 | 20.75 | +12.25 | 2.44x |
| 8 | Isoprene | 105.00 | 100.04 | -4.96 | 0.95x |
| 9 | Hexanal | 4.20 | 5.28 | +1.08 | 1.26x |
| 10 | Phenol | 3.50 | 4.42 | +0.92 | 1.26x |

Mechanisms:
- Acetone is predicted elevated in Type 2 diabetes mellitus; primarily via Ketone body metabolism (score=0.78; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.540); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 2 diabetes mellitus; primarily via Urea cycle / ammonia detoxification (score=0.00; genes CPS1, OTC, ASS1, ASL); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.540); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ethanol is predicted elevated in Type 2 diabetes mellitus; primarily via Glycolysis / Warburg metabolism (score=0.00; genes HK2, PKM, LDHA, PFKP); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.540); biosynthetic chain(s): ethanol_fermentation; key genetic/pathway nodes: HK2, PKM, LDHA, KRAS, NRAS, HRAS; Census enriched population example: oligodendrocyte in brain; KEGG C00469 (Glycolysis / Gluconeogenesis, Pyruvate metabolism); EC 1.1.1.1, 1.1.1.2, 1.1.1.71; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

### Type 1 diabetes (reference)

*Absolute insulin deficiency; high ketosis risk*

T1D atlas entry for comparison — acetone / ketone alcohols elevated; not the primary T2D vignette.

- Burden: stage III, frac=0.55, activity=1.9, prior_scale=1.0
- Comorbidities: —
- Ketone panel: **obvious** (acetone 33.54×, mean |log2FC| ketone=2.5465)
- Literature eval: PASS (dir_acc=1.0)

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 16000.00 | +15523.00 | 33.54x |
| 2 | Ammonia | 250.00 | 6464.76 | +6214.76 | 25.86x |
| 3 | 2-Butanone | 1.20 | 20.00 | +18.80 | 16.67x |
| 4 | Methanol | 160.00 | 175.07 | +15.07 | 1.09x |
| 5 | 2-Propanol | 25.00 | 37.56 | +12.56 | 1.50x |
| 6 | Ethanol | 120.00 | 123.64 | +3.64 | 1.03x |
| 7 | Isoprene | 105.00 | 107.18 | +2.18 | 1.02x |
| 8 | Phenol | 3.50 | 4.16 | +0.66 | 1.19x |
| 9 | n-Pentane | 8.50 | 9.06 | +0.56 | 1.07x |
| 10 | 2-Pentanone | 1.40 | 1.94 | +0.54 | 1.39x |

Mechanisms:
- Acetone is predicted elevated in Type 1 diabetes; primarily via Ketone body metabolism (score=0.36; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.045); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- Ammonia is predicted elevated in Type 1 diabetes; primarily via Urea cycle / ammonia detoxification (score=0.12; genes CPS1, OTC, ASS1, ASL); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.045); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: CPS1, OTC, ASS1, CASP3, CASP8, CASP9; Census enriched population example: oligodendrocyte in brain; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.
- 2-Butanone is predicted elevated in Type 1 diabetes; primarily via Ketone body metabolism (score=0.36; genes HMGCS2, HMGCL, BDH1, OXCT1); sourced from Ketogenic hepatocyte cells in liver (density×activity=1.045); biosynthetic chain(s): 2_butanone_fao; key genetic/pathway nodes: HMGCS2, HMGCL, BDH1, CPT1A, CPT2, ACADM; Census enriched population example: oligodendrocyte in brain; KEGG C02845; EC 4.1.2.46; VOLATILOME/HBDB breath compound + disease–VOC association in breathomics literature layer.

