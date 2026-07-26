# Zero-shot novel disease stack

## What shipped
1. Custom mechanism resolver (`knowledge/custom_resolver.py`) — token cues + category pathway templates (no VOC prior copy)
2. Offline disease→gene index (`zero_shot_disease_genes.json` + OT free-text)
3. Ontology nearest-neighbor transfer (pathway/site/cell-state donors only)
4. First-class CLI: `voc predict-novel` + `voc biomarker --description/--genes/--pathway-overrides`
5. Reliability eval: `voc eval-zero-shot-reliability` (held-out rare diseases + leave-disease-out)

## Results
- Holdout rare diseases: **12/12**
- Leave-disease-out prior direction: mean **0.854** (108/110)
- Stress-hard: **50/50**
