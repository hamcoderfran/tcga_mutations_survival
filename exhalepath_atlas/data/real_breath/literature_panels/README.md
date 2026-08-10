# Breath VOC literature panels

This directory holds curated, DOI-backed VOC direction/fold-change panels:

| File | Scope |
|------|--------|
| `priority10_voc_panels.json` | Airway/infection/cancer/HF/malaria/ARDS priority set (13 diseases) |
| `mental_health_voc_panels.json` | Schizophrenia, major depression, bipolar (+ honesty notes for thin MH conditions) |
| `magdeburg_ptrms_mz_map.json` | Magdeburg PTR-MS m/z → atlas VOC identity map |
| `mental_health_readiness.json` | Honest MH readiness / hole summary |

Loaders merge every `*_voc_panels.json` for claim ledger, literature overlays, and lit_compare.

## Curation rules

- Only real published **exhaled breath / alveolar gas / OralChroma breath VSC** studies.
- Numerical log2 fold changes only from ratios explicitly reported in cited sources.
- Direction-only placeholders: `+0.6` increased / `-0.4` decreased, marked `directional_only`.
- Fecal/urine VOC studies and unnamed stress-fingerprint VOCs are **not** promoted to measured panels.

## Mental health (honest bar)

**Measured panels (DOI-backed):**

- **Schizophrenia** — Magdeburg PTR-MS (Jiang 2022; JPN 2023; Molecules 2023) + Phillips pentane/CS2; butyric acid (m/z 90) in atlas catalog
- **Major depressive disorder** — Magdeburg MDD PTR-MS + Gbaoui; **butyric quantified** `log2(116/169)`
- **Bipolar** — OralChroma CH3SH quantified (`log2(18.62/9.45 ppb)`), plus directional pentane accent

**Thin / no measured exhaled panel (documented in `thin_evidence_conditions`):**

- Anxiety, PTSD — acute-stress volatilomics / trauma-cue odour work exists; no named disease fold-change panel suitable for `measured_log2fc`
- ADHD — no dedicated ADHD breath disease cohort found
- Autism — published VOC signatures are primarily fecal/urine microbiome volatilomes, not exhaled alveolar gas

Those four keep **de-cloned atlas priors** with explicit `atlas_source` labels instead of template clones.

**Open intensity data:** Magdeburg figshare `19181742` is a DOCX supplement only — no patient×VOC intensity matrix is bundled.

## Preferred public breath accessions (priority set context)

- Healthy baseline: `ST003200`
- Heart failure: `ST000587`
- Malaria: `ST000883`
- Cystic fibrosis: `ST001164`
- Pneumonia/CAP comparator: `ST002449`
- Asthma/COPD clinical breathomics: figshare `23522490.v6`
