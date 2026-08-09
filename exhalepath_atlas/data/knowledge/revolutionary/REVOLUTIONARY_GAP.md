# Revolutionary gap evaluation

Generated: 2026-08-09T19:32:40.715247+00:00

> Research / hypothesis-generation enablement. Not clinical diagnostic SOTA, not FDA claims, not a substitute for prospective multi-site trials.

## Verdict

ExhalePath already leads the *open research-reasoning* axis (mechanism fusion, NL, paper packs, claim ledger, closed loop). It does not lead instrument preprocessing or commercial on-breath identification (BreathXplorer / ptairMS / Owlstone). Revolutionary path = deepen the reasoning axis + bridge to peak tables, not chase fit-on-cohort AUROC.

**ExhalePath scorecard:** 29 / 36 (80.6%)

### Competitor totals (same scale)

- `breathxplorer`: 4
- `ptairms`: 8
- `owlstone_omni`: 17
- `scireports_xgboost_baseline`: 6

## Capability matrix

| Capability | Us | BreathXplorer | ptairMS | Owlstone | SciRep XGB |
|---|---:|---:|---:|---:|---:|
| Mechanism → VOC multi-model fusion | 3 | 0 | 0 | 1 | 0 |
| Patient-level nested CV + locked splits | 3 | 0 | 1 | 2 | 2 |
| Raw SESI/PTR/GC-MS peak picking | 0 | 3 | 3 | 3 | 0 |
| On-breath vs blank / MSI identification | 1 | 1 | 2 | 3 | 0 |
| Cross-disease differential diagnosis | 2 | 0 | 0 | 1 | 3 |
| Confounder residualization (age/sex/smoking) | 3 | 0 | 0 | 1 | 0 |
| Evidence-graded VOC↔disease claim ledger | 3 | 0 | 0 | 2 | 0 |
| Open BreathVOC interchange format | 3 | 0 | 1 | 1 | 0 |
| Mechanism ↔ cohort closed-loop eval | 3 | 0 | 0 | 0 | 0 |
| Multi-cohort leave-one-study-out | 2 | 0 | 1 | 2 | 0 |
| TRIPOD+AI / paper-pack automation | 3 | 0 | 0 | 1 | 1 |
| Natural-language clinical → VOC prediction | 3 | 0 | 0 | 0 | 0 |

## Live differentiators (this run)

- Claim ledger: **958** claims — {'atlas_prior': 823, 'directional_only': 57, 'directional_only_mixed': 1, 'mixed': 52, 'quantified': 23, 'quantified_null': 2}
- BreathVOC interchange: `BreathVOC-1.1` roundtrip_ok=True
- Sci Data DDx: fit-on-cohort macro AUROC **100.0%** vs mechanism hybrid **57.3%** (top1=28.1%)
- Age/sex residualization: mean VOC R²=0.048; residual mechanism macro=68.3%
- Mechanism↔cohort loop:
  - `ST000883` hybrid nested=53.3%
  - `ST000587` hybrid nested=60.0%

## Roadmap (what remains)

1. **Second open intensity cohort (malaria or pulmonary)** — LOSO / external validation is the #1 credibility gap vs closed multi-site e-nose trials (`labels_ready_intensity_pending`)
2. **Raw-spectrum bridge (mzML / Agilent .D → BreathVOC)** — Without peak picking we cannot ingest the majority of public deposits (`not_started`)
3. **Claim-ledger-grounded NL answers** — voc ask should cite quantified DOIs, not only atlas priors (`ledger_shipped_nl_ungrounded`)
4. **On-breath MSI enrichment vs Owlstone/NIST** — Identification confidence is what converts research → translational trust (`schema_ready`)
5. **Prospective multi-site protocol template** — Field failures (e-nose CRC external validation) are protocol failures (`partial_multisite_literature_only`)
6. **Do not: hill-climb ST000883 or Sci Data fit-on-cohort AUROC** — Sci Reports already owns the black-box ceiling (~0.998); our wedge is transferable mechanism DDx (`documented_stop_chasing`)

## What 'revolutionary' means here

- Be the open *reasoning layer* of breathomics — not another peak picker
- Make mechanism hypotheses falsifiable on locked patient matrices
- Make every VOC↔disease edge citeable with an evidence grade
- Make BreathVOC the interchange labs actually ship between instruments and papers
- Refuse AUROC theater on n≈35 / fit-on-cohort ceilings
