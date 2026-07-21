# ExhalePath

**Pathway-informed prediction of exhaled volatile organic compounds (VOCs) in parts-per-billion (ppb) for any disease — including tumor stage, site, histology, and driver mutations.**

ExhalePath connects large-scale disease molecular data (TCGA/GDC mutations + clinical covariates, Reactome pathways, Open Targets associations) to a mechanistic VOC emission model grounded in metabolic biochemistry and breath-biomarker literature. Optional gradient-boosting calibrators learn pathway→ppb mappings on a multi-cohort training corpus.

> Research / hypothesis-generation tool. Not a medical device. Predicted ppb values should be validated against breath GC-MS / PTR-MS cohorts before clinical use.

## Why this exists

Breath VOC panels are promising non-invasive biomarkers, but most work is empirical and disease-specific. ExhalePath asks a harder, more transferable question:

> Given a disease context (and optional tumor stage / site / histology / mutated genes), which metabolic pathways are perturbed, and what exhaled VOC concentration shifts (ppb) does that biochemistry imply?

## Architecture

```
Disease + tumor context + mutations
        │
        ▼
 Pathway activity scoring  ← Reactome gene sets + Open Targets
        │
        ▼
 Mechanistic VOC emission priors (pathway × VOC coefficients)
        │
        ▼
 Optional ML calibrators (per-VOC HistGradientBoosting)
        │
        ▼
 Exhaled VOC panel: healthy ppb → predicted ppb, Δppb, fold-change, CIs
```

### Data scale

`build-corpus` harvests clinical cases and mutation edges across dozens of GDC projects, expands pathway gene sets via Reactome, and materializes case×pathway×VOC training matrices that reach **hundreds of thousands to millions of rows** depending on project count.

VOC regression targets in the default corpus are **mechanistic-prior synthetic labels** (literature disease fold-changes + pathway scores + stage effects + noise) so the calibrator can be trained without a private breath biobank. Swap in measured breath cohorts when you have them — the feature schema is unchanged.

## Install

```bash
cd exhalepath
pip install -e ".[dev]"
```

## Quickstart

```bash
# 1) Mechanistic prediction (no download required)
python -m exhalepath predict "pancreatic adenocarcinoma" \
  --stage III --site pancreas --histology adenocarcinoma \
  --genes KRAS,TP53,CDKN2A \
  --out-dir ./runs/paad_demo

# 2a) Build a large offline training corpus (no network; hundreds of thousands of rows)
python -m exhalepath build-corpus --offline-demo --demo-cases 5000

# 2b) Or pull live multi-cohort TCGA/GDC + Reactome (slower, real molecular data)
python -m exhalepath build-corpus --max-projects 8

# 3) Train per-VOC calibrators
python -m exhalepath train

# 4) Predict again with the warm model
python -m exhalepath predict LUAD --stage IV --site lung --genes TP53,KRAS --smoking current
```

### Python API

```python
from exhalepath import ExhalePathPredictor, DiseaseQuery, TumorContext

predictor = ExhalePathPredictor()
result = predictor.predict(
    DiseaseQuery(
        disease="breast cancer",
        tumor=TumorContext(stage="IIA", primary_site="breast", histology="invasive ductal carcinoma"),
        mutated_genes=["PIK3CA", "TP53"],
        sex="female",
        age_years=58,
    )
)
print(result.top(10))
```

## CLI

| Command | Purpose |
|---|---|
| `predict` | VOC ppb panel for a disease / tumor context |
| `build-corpus` | Pull GDC + pathway data; write training matrices |
| `train` | Fit per-VOC calibrators |
| `list-diseases` | Curated disease atlas |
| `list-vocs` | VOC catalog with healthy breath baselines |
| `audit` | Literature accuracy, calibrator holdout, zero-shot + stress tests |

```bash
python -m exhalepath audit --out runs/audit/audit_report.json
```

## Knowledge bases

Bundled under `data/knowledge/`:

- `voc_catalog.json` — 20 endogenous breath VOCs with healthy ppb medians/ranges + literature anchors
- `pathway_voc_map.json` — metabolic pathways → VOC emission coefficients (lipid peroxidation, Warburg glycolysis, mevalonate/isoprene, methionine→DMS, urea cycle→ammonia, …)
- `disease_voc_priors.json` — disease-level pathway biases + VOC log2FC priors (cancers + metabolic / inflammatory diseases)

Any unrecognized disease still runs: Open Targets enrichments (when online) + generic pathway scoring.

### Dysbiosis & brain conditions

First-class atlas entries (not just zero-shot fallbacks) include:

- **Microbiome / dysbiosis:** gut dysbiosis, SIBO, *C. difficile*, *H. pylori*, IBD  
- **Neurological:** Alzheimer, Parkinson, depression, schizophrenia, epilepsy, MS, TBI, autism spectrum, glioblastoma  

These use dedicated pathways (`gut_microbiome_fermentation`, `microbial_proteolysis_putrefaction`, `neuroinflammation`, `neurotransmitter_metabolism`, `brain_energy_metabolism`) and microbiome VOCs (`indole`, `phenol`, `dimethyl_disulfide`). Predictions for these categories up-weight the mechanistic pathway model over the cancer-trained calibrator.

## Tumor covariates

`TumorContext` modulates pathway scores and VOC output via:

- **Stage** (I→IV burden multipliers)
- **Primary site / position** (liver boosts sulfur/ammonia pathways; lung boosts peroxidation / CYP, etc.)
- **Histology** (e.g. squamous vs adenocarcinoma)
- **Metastatic flag** and optional tumor-burden proxy

## Outputs

Prediction reports include:

- `voc_predictions.csv` — healthy ppb, predicted ppb, Δppb, fold-change, CI, pathway drivers
- `pathway_scores.csv` — ranked metabolic dysregulation
- `prediction_bundle.json` — full provenance
- `voc_fold_changes.png` — ranked VOC shift plot

## Scientific caveats

- Breath VOC concentrations vary with diet, microbiome, smoking, sampling protocol, and instrument.
- Gene-length artifacts in “top mutated gene” lists are handled by using curated metabolic pathway gene sets, not raw frequency ranks alone.
- Synthetic training labels enable end-to-end ML plumbing; for publication-grade calibration, replace `voc_training_targets.csv` with measured breath concentrations joined on disease/stage/molecular features.

## Related

This repository also contains `project1_tcga_mutations_survival/` (mutation frequency × overall survival for LUAD/BRCA/PAAD/COAD). ExhalePath complements that work by mapping the same molecular disease contexts onto exhaled VOC chemistry.
