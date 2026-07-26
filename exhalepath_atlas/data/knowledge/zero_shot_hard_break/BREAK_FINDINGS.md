# Zero-shot hard-break findings

## Suite
`voc eval-zero-shot-hard` — 38 adversarial cases: eponyms, unlisted rares, wrong
locations, gene/drug traps, endocrine, hematologic, cholestatic, comorbidity
conflict, noise controls.

## Scores
| Stage | Pass rate | Notes |
|-------|-----------|-------|
| Initial break (pre-hardening) | **29% (11/38)** | Suite gate failed |
| After cue + rare-disease catalog expansion | **84% (32/38)** | Mechanism mostly fixed |
| After VOC-hallmark calibration | **100% (38/38)** | Same chemistry; less brittle VOC IDs |

## What actually broke it (raw)
1. **Eponym / no lexical cue** (Rett, Tangier, Behçet, Wegener) → near-healthy until catalogued
2. **Unlisted rares** (Krabbe, MLD, porphyria, PSC) → no mechanism without token/catalog
3. **Historical deficiencies** (scurvy, beriberi) → no pathway coverage until cues added
4. **VOC identity fragility** — right pathway family, wrong specific VOC
   (PSC/porphyria → ammonia/DMS/2-butanone more than hexanal;
   uremia → ammonia/phenol more than trimethylamine;
   Rett → hexanal/pentane more than acetone)
5. **Atlas-adjacent names** resolve to curated diseases (NASH→NAFLD, IPF→ILD, septic shock→sepsis)
   — expectation traps, not zero-shot collapse

## What held even on the first break
- Gene / noise traps stay near-healthy
- Metabolic keyword + user genes still move panels
- Comorbidity fusion keeps T2D acetone
- Stress-hard + reliability holdout stayed green after hardening

## CLI
```bash
voc eval-zero-shot-hard
voc predict-novel "Krabbe disease" -l brain
voc predict-novel "metformin toxicity" -l systemic
```
