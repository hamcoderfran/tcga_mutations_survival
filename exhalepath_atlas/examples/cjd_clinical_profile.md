# Creutzfeldt–Jakob disease — clinical VOC profile

**Patient:** 62.0yo female, genes=PRNP, location=brain.

> Research hypothesis only — not a diagnostic test. No peer-reviewed exhaled-VOC signature is established for CJD; gold-standard workup remains MRI DWI, EEG, CSF RT-QuIC / 14-3-3 / t-tau. Model absolute folds are upper-bound estimates; interpret direction and tempo.

## Bottom line

Incubating phase stays `subtle_or_absent` (mean |log2FC|=0.0667). First meaningful VOC hint: **Prodromal / very early clinical** (0–4 weeks from clinical onset), verdict `obvious`. Becomes obvious on the oxidative breath panel: **Prodromal / very early clinical** (0–4 weeks from clinical onset). Even when model VOCs move early after clinical onset, they are a nonspecific oxidative/necrosis signature — not prion-specific. MRI DWI + CSF RT-QuIC remain the appropriate early diagnostic path.

- First meaningful hint phase: `prodromal`
- First obvious phase: `prodromal`
- CJD vs AD (matched early burden) mean |log2FC|: 1.3772 vs 1.2109

## Phases

### Incubating / preclinical (months–years before onset (silent))

Prion propagation without frank dementia; MRI/EEG normal; CSF RT-QuIC may convert late in incubation. Breath VOCs expected near baseline.

- Model burden: stage I, affected_fraction=0.02, activity=0.7, prior_scale=0.15
- Oxidative panel verdict: **subtle_or_absent** (mean |log2FC|=0.0667)
- Hint VOCs: pentane
- Obvious VOCs: isoprene

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Isoprene | 105.00 | 62.37 | -42.63 | 0.59x |
| 2 | Acetone | 477.00 | 466.62 | -10.38 | 0.98x |
| 3 | Methanol | 160.00 | 150.04 | -9.96 | 0.94x |
| 4 | Ammonia | 250.00 | 244.41 | -5.59 | 0.98x |
| 5 | n-Pentane | 8.50 | 10.80 | +2.30 | 1.27x |
| 6 | Ethanol | 120.00 | 117.84 | -2.16 | 0.98x |
| 7 | 1-Propanol | 12.00 | 11.33 | -0.67 | 0.94x |
| 8 | Acetaldehyde | 15.00 | 14.43 | -0.57 | 0.96x |
| 9 | 2-Propanol | 25.00 | 24.45 | -0.55 | 0.98x |
| 10 | Hydrogen sulfide | 8.00 | 7.66 | -0.34 | 0.96x |

Mechanisms:
- Isoprene is predicted reduced in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.03; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.068); biosynthetic chain(s): isoprene_mevalonate; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, HMGCR, MVK, MVD; KEGG C16521 (Terpenoid backbone biosynthesis, Biosynthesis of terpenoids and steroids); EC 4.2.3.27; VOLATILOME/HBDB breath compound.
- Acetone is predicted reduced in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.03; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.068); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, CASP3, CASP8, CASP9; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound.
- Methanol is predicted reduced in Creutzfeldt-Jakob disease; primarily via One-carbon / folate metabolism (score=0.00; genes MTHFR, MTR, SHMT1, SHMT2); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.068); key genetic/pathway nodes: MTHFR, MTR, SHMT1, FFAR2, FFAR3, SLC5A8; KEGG C00132 (Butanoate metabolism, Methane metabolism); EC 1.1.1.244, 1.1.2.7, 1.1.2.10; VOLATILOME/HBDB breath compound.
- Ammonia is predicted reduced in Creutzfeldt-Jakob disease; primarily via Neuroinflammation / microglial activation (score=0.05; genes TNF, IL1B, IL6, NFKB1); sourced from Ketogenic hepatocyte cells in liver (density×activity=0.068); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: TNF, IL1B, IL6, MAOA, MAOB, COMT; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound.

### Prodromal / very early clinical (0–4 weeks from clinical onset)

Subtle cognitive slowing, anxiety/depression, sleep disturbance; MRI/EEG often still nondiagnostic; RT-QuIC may already be positive.

- Model burden: stage I, affected_fraction=0.1, activity=1.0, prior_scale=0.55
- Oxidative panel verdict: **obvious** (mean |log2FC|=0.751)
- Hint VOCs: nonanal, ethane
- Obvious VOCs: hexanal, pentane

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Acetone | 477.00 | 566.03 | +89.03 | 1.19x |
| 2 | n-Pentane | 8.50 | 75.26 | +66.76 | 8.85x |
| 3 | Ammonia | 250.00 | 287.78 | +37.78 | 1.15x |
| 4 | Hexanal | 4.20 | 21.99 | +17.79 | 5.24x |
| 5 | Isoprene | 105.00 | 87.58 | -17.42 | 0.83x |
| 6 | Methanol | 160.00 | 151.66 | -8.34 | 0.95x |
| 7 | Nonanal | 5.00 | 7.23 | +2.23 | 1.45x |
| 8 | Ethanol | 120.00 | 118.20 | -1.80 | 0.98x |
| 9 | Acetaldehyde | 15.00 | 16.46 | +1.46 | 1.10x |
| 10 | 1-Propanol | 12.00 | 12.62 | +0.62 | 1.05x |

Mechanisms:
- Acetone is predicted elevated in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.15; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.100); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, CASP3, CASP8, CASP9; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound.
- n-Pentane is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.18; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.100); biosynthetic chain(s): pentane_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C13388; VOLATILOME/HBDB breath compound.
- Ammonia is predicted elevated in Creutzfeldt-Jakob disease; primarily via Neuroinflammation / microglial activation (score=0.27; genes TNF, IL1B, IL6, NFKB1); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.100); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: TNF, IL1B, IL6, CASP3, CASP8, CASP9; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound.
- Hexanal is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.18; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.100); biosynthetic chain(s): hexanal_lipid_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C02233; VOLATILOME/HBDB breath compound.

### Early clinical (4–8 weeks)

Rapid cognitive decline, ataxia, emerging myoclonus; diffusion MRI cortical/basal-ganglia ribboning often appears here.

- Model burden: stage II, affected_fraction=0.25, activity=1.4, prior_scale=0.85
- Oxidative panel verdict: **obvious** (mean |log2FC|=1.3772)
- Hint VOCs: octanal, ethane, acetone, ammonia
- Obvious VOCs: hexanal, heptanal, nonanal, pentane

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | n-Pentane | 8.50 | 240.00 | +231.50 | 28.24x |
| 2 | Acetone | 477.00 | 644.21 | +167.21 | 1.35x |
| 3 | Hexanal | 4.20 | 120.00 | +115.80 | 28.57x |
| 4 | Ammonia | 250.00 | 324.20 | +74.20 | 1.30x |
| 5 | Isoprene | 105.00 | 91.50 | -13.50 | 0.87x |
| 6 | Nonanal | 5.00 | 9.93 | +4.93 | 1.99x |
| 7 | Methanol | 160.00 | 163.84 | +3.84 | 1.02x |
| 8 | Acetaldehyde | 15.00 | 17.90 | +2.90 | 1.19x |
| 9 | Heptanal | 3.10 | 5.36 | +2.26 | 1.73x |
| 10 | 1-Propanol | 12.00 | 13.42 | +1.42 | 1.12x |

Mechanisms:
- n-Pentane is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.39; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.350); biosynthetic chain(s): pentane_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C13388; VOLATILOME/HBDB breath compound.
- Acetone is predicted elevated in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.32; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.350); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, CASP3, CASP8, CASP9; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound.
- Hexanal is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.39; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.350); biosynthetic chain(s): hexanal_lipid_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C02233; VOLATILOME/HBDB breath compound.
- Ammonia is predicted elevated in Creutzfeldt-Jakob disease; primarily via Neuroinflammation / microglial activation (score=0.60; genes TNF, IL1B, IL6, NFKB1); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.350); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: TNF, IL1B, IL6, CASP3, CASP8, CASP9; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound.

### Fulminant (8–16 weeks)

Dementia, startle myoclonus, pyramidal/extrapyramidal signs; periodic sharp waves on EEG; clinical diagnosis usually clear.

- Model burden: stage III, affected_fraction=0.5, activity=1.95, prior_scale=1.0
- Oxidative panel verdict: **obvious** (mean |log2FC|=1.5431)
- Hint VOCs: octanal, acetone, ammonia, hydrogen_sulfide, acetaldehyde
- Obvious VOCs: hexanal, heptanal, nonanal, pentane, ethane

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | n-Pentane | 8.50 | 240.00 | +231.50 | 28.24x |
| 2 | Acetone | 477.00 | 708.36 | +231.36 | 1.49x |
| 3 | Hexanal | 4.20 | 120.00 | +115.80 | 28.57x |
| 4 | Ammonia | 250.00 | 355.27 | +105.27 | 1.42x |
| 5 | Methanol | 160.00 | 178.83 | +18.83 | 1.12x |
| 6 | Nonanal | 5.00 | 13.24 | +8.24 | 2.65x |
| 7 | Isoprene | 105.00 | 100.11 | -4.89 | 0.95x |
| 8 | Ethanol | 120.00 | 124.29 | +4.29 | 1.04x |
| 9 | Acetaldehyde | 15.00 | 18.95 | +3.95 | 1.26x |
| 10 | Heptanal | 3.10 | 6.63 | +3.53 | 2.14x |

Mechanisms:
- n-Pentane is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.65; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.975); biosynthetic chain(s): pentane_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C13388; VOLATILOME/HBDB breath compound.
- Acetone is predicted elevated in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.53; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.975); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, CASP3, CASP8, CASP9; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound.
- Hexanal is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.65; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.975); biosynthetic chain(s): hexanal_lipid_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C02233; VOLATILOME/HBDB breath compound.
- Ammonia is predicted elevated in Creutzfeldt-Jakob disease; primarily via Neuroinflammation / microglial activation (score=1.01; genes TNF, IL1B, IL6, NFKB1); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=0.975); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: TNF, IL1B, IL6, CASP3, CASP8, CASP9; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound.

### Terminal / akinetic mute (16–26 weeks (typical sCJD))

Akinetic mutism with severe neuronal loss and gliosis; median sCJD survival ~4–6 months from clinical onset.

- Model burden: stage IV, affected_fraction=0.75, activity=2.4, prior_scale=1.1
- Oxidative panel verdict: **obvious** (mean |log2FC|=1.7274)
- Hint VOCs: hydrogen_sulfide, acetaldehyde
- Obvious VOCs: hexanal, heptanal, nonanal, octanal, pentane, ethane, acetone, ammonia

| Rank | VOC | Healthy ppb | Pred ppb | Δppb | Fold |
|---:|---|---:|---:|---:|---:|
| 1 | Isoprene | 105.00 | 1400.93 | +1295.93 | 13.34x |
| 2 | Acetone | 477.00 | 774.49 | +297.49 | 1.62x |
| 3 | n-Pentane | 8.50 | 240.00 | +231.50 | 28.24x |
| 4 | Ammonia | 250.00 | 410.68 | +160.68 | 1.64x |
| 5 | Hexanal | 4.20 | 120.00 | +115.80 | 28.57x |
| 6 | Methanol | 160.00 | 195.31 | +35.31 | 1.22x |
| 7 | Nonanal | 5.00 | 17.44 | +12.44 | 3.49x |
| 8 | Ethanol | 120.00 | 127.93 | +7.93 | 1.07x |
| 9 | Acetaldehyde | 15.00 | 19.98 | +4.98 | 1.33x |
| 10 | Heptanal | 3.10 | 8.02 | +4.92 | 2.59x |

Mechanisms:
- Isoprene is predicted elevated in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.77; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=1.800); biosynthetic chain(s): isoprene_mevalonate; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, HMGCR, MVK, MVD; KEGG C16521 (Terpenoid backbone biosynthesis, Biosynthesis of terpenoids and steroids); EC 4.2.3.27; VOLATILOME/HBDB breath compound.
- Acetone is predicted elevated in Creutzfeldt-Jakob disease; primarily via Brain energy / mitochondrial metabolism (score=0.77; genes PPARGC1A, TFAM, SOD2, PARK7); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=1.800); biosynthetic chain(s): acetone_ketogenesis; key genetic/pathway nodes: PPARGC1A, TFAM, SOD2, CASP3, CASP8, CASP9; KEGG C00207 (Butanoate metabolism, Metabolic pathways); EC 1.1.1.80, 1.14.13.226, 1.14.14.141; VOLATILOME/HBDB breath compound.
- n-Pentane is predicted elevated in Creutzfeldt-Jakob disease; primarily via Lipid peroxidation / oxidative stress (score=0.95; genes GPX4, SOD2, CAT, ACSL4); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=1.800); biosynthetic chain(s): pentane_peroxidation; key genetic/pathway nodes: GPX4, SOD2, CAT, TNF, IL1B, IL6; KEGG C13388; VOLATILOME/HBDB breath compound.
- Ammonia is predicted elevated in Creutzfeldt-Jakob disease; primarily via Neuroinflammation / microglial activation (score=1.48; genes TNF, IL1B, IL6, NFKB1); sourced from Oxidative-stress / lipid-peroxidizing cell cells in multi (density×activity=1.800); biosynthetic chain(s): ammonia_urea; key genetic/pathway nodes: TNF, IL1B, IL6, CASP3, CASP8, CASP9; KEGG C00014 (Arginine biosynthesis, Purine metabolism); EC 1.3.7.8, 1.4.1.1, 1.4.1.2; VOLATILOME/HBDB breath compound.

