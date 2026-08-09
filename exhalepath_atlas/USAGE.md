# ExhalePath Atlas / `voc-breath` — usage guide

Research / hypothesis-generation for exhaled VOC biology. **Not a medical device.**

Package version **1.6.0** on `main`: 20-model Great Disease Stack + GC-MS patient diagnostic research layer + analysis exports (literature overlay, Methods, GraphPad CSV).

## Install (ready to run — no training required)

```bash
pip install "voc-breath[stack] @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

Or from a clone:

```bash
git clone https://github.com/hamcoderfran/tcga_mutations_survival.git
cd tcga_mutations_survival/exhalepath_atlas
pip install -e ".[dev,stack]"
# or: bash install_voc.sh
```

You get the `voc` command (aliases: `exhalepath`, `voc-breath`).

**Bundled and ready:**
- Knowledge atlas + disease/VOC/pathway priors
- Trained `voc_calibrator.joblib` + `chembl_aux_model.joblib`
- Great-stack prior packs (PrimeKG-lite, OPERA panel, flux proxy, AGORA routes, calibrated fusion weights)
- Optional: `pip install "voc-breath[gsmm]"` for full COBRApy FBA (Human-GEM SBML not bundled)

## Everyday prediction

```bash
# One-line biomarker pack → runs/voc_*/REPORT.html
voc "depression" -l brain -c obesity --age 24 --sex male

# 20-model fused stack → runs/stack_*/STACK_REPORT.html
voc stack "depression" -l brain -c obesity --age 24 --sex male

# Naturalistic clinical text
voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations"
voc stack --nl "Maple syrup urine disease, genes BCKDHA BCKDHB"

# Novel / zero-shot
voc predict-novel "unknown phenotype" -l liver -d "ketotic mitochondrial stress" --genes HMGCS2,CPT1A
```

### Analysis pack contents (what you get)

| Artifact | Biomarker (`voc`) | Stack (`voc stack`) |
|---|---|---|
| HTML report | `REPORT.html` | `STACK_REPORT.html` |
| Markdown | `REPORT.md` | `STACK_REPORT.md` |
| Dashboard PNG | `dashboard.png` | `stack_dashboard.png` + vote heatmap |
| VOC tables | `top_voc_biomarkers.csv` | `fused_vocs.csv` (+ epistemic CI) |
| **Literature overlay** | `literature_overlay.csv` | same |
| **GraphPad/Prism long CSV** | `graphpad_voc_long.csv` | same |
| **Methods draft** | `METHODS.md` | same |
| **Next experiments** | `NEXT_EXPERIMENTS.md` | same |
| Machine JSON | `prediction_bundle.json` | `STACK_RESULT.json` |

## GC-MS diagnostic research

```bash
# Patient-level AUROC on bundled MW malaria breath GC-MS
voc eval-patient-diagnostic --study ST000883 --signature hybrid

# Compare hybrid / stack / literature signatures
voc eval-patient-diagnostic --benchmark --study ST000883

# All bundled studies
voc eval-patient-diagnostic --all --signature stack

# Preregistration-style locked splits (cite the sha256)
voc lock-split --study ST000883 --out data/knowledge/gcms_splits/ST000883_split_v1.json

# Score an observed VOC vector against a disease signature
voc score-sample malaria --signature hybrid \
  --vocs '{"acetone":0.4,"pentane":0.5,"isoprene":-0.3,"benzene":0.6}'
```

## Evaluation / readiness

```bash
voc eval-stack-holdout              # stack vs hybrid + PatientTemplate breaks
voc eval-implementation-readiness
voc eval-public-breath
voc eval-priority10
```

Canonical archived reports: [`data/knowledge/RESULTS_INDEX.md`](data/knowledge/RESULTS_INDEX.md).

## Optional retrain (not required)

```bash
voc build-real-corpus
voc train                 # refreshes voc_calibrator.joblib
voc train-chembl          # refreshes chembl aux model
```

## How this builds on what already exists

| Existing ecosystem | How ExhalePath / voc-breath extends it |
|---|---|
| Metabolomics Workbench / Sci Data peak tables | Ingests real cohorts; adds **patient-level AUROC + locked splits** (not only group medians) |
| Literature breath VOC panels (HBDB-style, curated DOIs) | Uses them as priors **and** as an explicit **literature overlay** in every report (with circularity warnings) |
| Physiology / Farhi / cell-state models | Hybrid engine + calibrator residual; fused with genetics/flux/ADME/microbiome heads |
| Knowledge graphs (PrimeKG / OmniPath / Open Targets) | Lightweight offline packs in the 20-model stack — no heavyweight graph DB required to start |
| Owlstone / Peppermint / TRIPOD+AI reporting norms | Ships TRIPOD stub + Methods.md + nested vs non-nested optimism gap for GC-MS studies |
| Generic ML on peak tables | Mechanism signatures as **transferable templates** + logistic baseline ceiling on the same features |

## Making analysis even more useful (roadmap hooks already started)

1. **Surface uncertainty everywhere** — epistemic CI columns in stack HTML/MD (done); biomarker HTML next  
2. **Literature overlay + DOIs in every pack** — done for stack + biomarker exports  
3. **Prism/GraphPad long CSV** — done (`graphpad_voc_long.csv`)  
4. **Next-experiment bullets** — done from disagreements / missing mappings / model gaps  
5. **Next ups:** Sci Data per-sample adapters, smoking-stratified AUCs when metadata exist, mzTab-M export, blank-ratio filters  

## Honesty bar

- Cutting-edge **open systems / research enablement** stack  
- **Not** clinical SOTA vs prospective multi-site GC-MS/PTR trials  
- See [`great_disease_stack/SOTA.md`](great_disease_stack/SOTA.md) and [`RESEARCH.md`](RESEARCH.md)
