# Malaria breath diagnostic upgrade

Generated: 2026-08-09T19:13:30.214066+00:00
Study: `ST000883` · feature_map=`malaria_lit`

> Transferable signature vs fit-on-cohort sparse — research enablement only.

## Feature remapping

- Library metabolites mapped: {'n_library_metabolites': 147, 'n_library_metabolites_mapped': 70}
- Atlas-only VOC count: 17
- Malaria-lit VOC count: 25
- Newly mapped: alpha_pinene, cyclohexanone, delta_3_carene, methyl_undecane, tridecane, trimethyl_alkene, trimethyl_decane, xylene

## Transferable (hybrid + malaria literature signature)

- Nested AUROC: **58.3%**
- Non-nested AUROC: **58.8%**
- Optimism gap: **0.5%**
- AUPRC: **64.1%** CI95=(0.41965453725815444, 0.8373129258547501)

### Fixed-sensitivity operating points (signature)

- sens≥0.8: achieved sens=82.4%, spec=27.8% CI95=(0.02065217391304354, 0.7030882352941176)
- sens≥0.9: achieved sens=100.0%, spec=11.1% CI95=(0.0, 0.4362092391304347)
- sens≥0.95: achieved sens=100.0%, spec=11.1% CI95=(0.0, 0.375)

## Fit-on-cohort nested sparse (ceiling)

- SelectKBest nested AUROC: **73.3%** (gap=22.7%)
- L1 nested AUROC: **65.0%** (gap=32.4%)
- Consensus features (kbest): 2_pentanone, methyl_undecane, propanol, tridecane, cyclohexanone, isoprene, nonanal, delta_3_carene

## Learning curve (n-limitation)

| n | mean nested-ish AUROC | std |
|---:|---:|---:|
| 12 | 51.9% | 0.16691387915795244 |
| 16 | 48.7% | 0.16059944293165793 |
| 20 | 57.4% | 0.1473066378064415 |
| 24 | 57.6% | 0.15649633360100024 |
| 28 | 58.7% | 0.11368804784242811 |
| 32 | 67.8% | 0.10230304775639347 |
| 35 | 69.9% | 0.09030583609558028 |

## External cohorts

- Catalog: `runs/malaria_diagnostic/external/EXTERNAL_MALARIA_CATALOG.json`
- Pooled: ST000883 only (no open second intensity table)

## Caveats

- ST000883 n≈35 — AUROC CIs remain wide (JBR: curves flatten near n≈50).
- Schaber 2018: thioethers largely absent; terpenes (pinene/carene) are the remap unlock.
- CSIRO CHMI and JID 2024 Malawi lack open intensity tables — catalog recorded for when deposits appear.
- Fit-on-cohort sparse AUROC is a ceiling, not a transferable clinical claim.
- Research enablement only — not a diagnostic device claim.
