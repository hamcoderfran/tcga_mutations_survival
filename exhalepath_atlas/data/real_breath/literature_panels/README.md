# Priority breath VOC literature panels

This directory contains `priority10_voc_panels.json`, a conservative, literature-backed VOC fold-change panel for 13 diseases:

- asthma
- copd
- covid19
- pneumonia_bacterial (VAP/HAP/CAP)
- tuberculosis
- cystic_fibrosis
- sleep_apnea (OSA)
- cancer_stomach (gastric)
- head_neck_cancer
- cancer_prostate
- heart_failure
- malaria
- ards

## Curation rules

- Only real published breath, exhaled-gas, ventilator-exhalate, or exhaled breath condensate studies were used.
- Numerical log2 fold changes are calculated only from ratios explicitly reported in cited sources.
- If a cited source reports only a compound direction, the JSON uses conservative placeholders:
  - `+0.6` for increased
  - `-0.4` for decreased
- Direction-only placeholders are marked in `voc_evidence` with `"evidence": "directional_only"`.
- Mixed or null compound-specific evidence is marked as `mixed`, `quantified_null`, or `quantified_null_or_no_difference`.
- Healthy context is included as Metabolomics Workbench `ST003200` (504 healthy subjects, PTR-TOF-MS), but disease fold changes use the comparator published in each disease study when available.

## Evidence-grade counts

Strict disease-level `evidence_grade` counts:

- Quantified: 4 diseases
- Directional-only: 6 diseases
- Mixed: 3 diseases

Diseases with at least one explicit published numerical ratio or median contrast in `measured_log2fc`: asthma, copd, cystic_fibrosis, head_neck_cancer, cancer_prostate, heart_failure.

## Preferred public breath accessions included as context

- Healthy baseline: `ST003200`
- Heart failure: `ST000587`
- Malaria: `ST000883`
- Cystic fibrosis: `ST001164`
- Pneumonia/CAP comparator volatile dataset: `ST002449`
- Asthma/COPD clinical breathomics: figshare `23522490.v6`
