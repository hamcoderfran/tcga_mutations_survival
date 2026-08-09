# Confounder / demographics pack

Generated: 2026-08-09T19:27:02.558780+00:00

> ST003200 is healthy-only — VOC→smoking/sex AUROC measures confounder signal in the feature space, not disease discrimination. Sci Data OVR strata are pulmonary case-mix, not case–control vs healthy.

## ST003200 healthy PTR (n=504)

- Features: 8 VOCs — acetaldehyde, acetone, ethanol, isoprene, methanol, phenol, propanol, toluene
- Smoking current vs never — confounder-proxy AUROC: **72.5%** (n=490)
- Sex male vs female — confounder-proxy AUROC: **63.2%** (n=504)

## Sci Data GC-MS age/sex strata

- Cohorts reported: 3
- Smoking available: False
- `SCIDATA2024_asthma_ovr` n=121 overall AUROC=51.9%
- `SCIDATA2024_copd_ovr` n=121 overall AUROC=58.0%
- `SCIDATA2024_bronchiectasis_ovr` n=121 overall AUROC=71.1%

## Caveats

- ST003200: no disease labels — use as confounder effect-size reference only.
- Sci Data: one-vs-rest pulmonary cohorts; smoking metadata absent.
- Do not treat confounder-proxy AUROC as diagnostic performance.
