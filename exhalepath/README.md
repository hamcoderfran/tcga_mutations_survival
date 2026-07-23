# ExhalePath

**Whole-body exhaled VOC biomarker prediction — think AlphaFold-style structure prediction, but for breath: disease + any location of affected cells → ranked top-50 VOC concentration shifts (ppb).**

ExhalePath models the full human body tissue map (~66 Census-informed sites), ~100 high-burden diseases, and a **50-VOC exhaled panel**. Pathway chains × affected cell states → tissue production → blood transfer → Farhi alveolar release → exhaled ppb. Optional ML calibrators refine well-studied VOCs.

> Research / hypothesis-generation tool. Not a medical device. Predicted ppb values should be validated against breath GC-MS / PTR-MS cohorts before clinical use.

## Why this exists

Breath VOC panels are promising non-invasive biomarkers, but most work is empirical and disease-specific. ExhalePath asks a harder, more transferable question:

> Given a disease context (and optional tumor stage / site / histology / mutated genes), which metabolic pathways are perturbed, and what exhaled VOC concentration shifts (ppb) does that biochemistry imply?

## Architecture

```
Disease + tumor context + mutations (+ optional scRNA cell fractions)
        │
        ▼
 Pathway activity scoring  ← Reactome gene sets + Open Targets
        │
        ├─► Multi-step VOC pathway chains (e.g. ketogenesis→acetone, Warburg→acetaldehyde)
        │
        ├─► Affected cell states: density × metabolic activity
        │         │
        │         ▼
        │   Tissue VOC production
        │         │
        │         ▼
        │   Blood transfer (perfusion, hepatic first-pass)
        │         │
        │         ▼
        │   Farhi alveolar release (λ blood:air, VA, Q) → ppb
        │
        ▼
 Optional ML calibrators (hybrid blend)
        │
        ▼
 Exhaled VOC panel: healthy ppb → predicted ppb, Δppb, fold-change, CIs + physio traces
```

### Physiology mode

```bash
# Full cell→blood→alveolar model (diabetes acetone / cancer acetaldehyde chains)
python -m exhalepath predict "type 2 diabetes" --mode physiology --out-dir ./runs/t2d_physio
python -m exhalepath predict LUAD --mode physiology --stage III --site lung \
  --genes KRAS,TP53,HK2,LDHA --out-dir ./runs/luad_physio

# Override cell-state fractions from single-cell analysis
python -m exhalepath predict PAAD --mode physiology \
  --cell-fractions tumor_epithelial_warburg=0.55,hepatocyte_ketogenic=0.4
```

### Public breath validation + mechanism “WHY” packs

```bash
# Download Scientific Data 2024 breathomics (asthma/COPD/bronchiectasis) + build cases
python -m exhalepath harvest-public-breath
python -m exhalepath eval-public-breath --top-k 15 --out-dir runs/public_breath_eval

# Explain VOC changes via pathways → Census cell populations → driver genes
python -m exhalepath explain "lung adenocarcinoma" --location lung --genes KRAS,TP53 --top 10

# Build mechanism packs for the full ~100-disease atlas
python -m exhalepath build-mechanism-packs --top-vocs 15
```

Each VOC mechanism states **why** it moves: dysregulated pathways, biosynthetic chains, affected cell states (with Census fractions / top cell types), and genetic nodes.

### ChEMBL chemogenomic training (millions of activity rows)

ChEMBL (∼2.9M compounds, ∼**24M** bioactivities) does **not** label exhaled ppb. It trains the **chemogenomic middle layer**: which chemicals potently modulate VOC-pathway enzymes, plus VOC physicochemical priors (AlogP → blood–air λ hints).

```bash
# Offline multi-million-row drill (no network)
python -m exhalepath harvest-chembl --offline-demo --demo-rows 1000000
python -m exhalepath train-chembl --max-rows 2000000

# Live API harvest for pathway seed genes (cap as needed; raise for fuller slice)
python -m exhalepath harvest-chembl --max-per-target 5000 --max-rows 500000
python -m exhalepath train-chembl
```

Outputs: `data/chembl/chembl_pathway_activities.csv`, distilled `data/knowledge/chembl_pathway_priors.json`, aux model `data/models/chembl_aux_model.joblib`. Predict-time features include `chembl_lig_*` / `chembl_voc_ligandability`.

For the full database dump, download ChEMBL SQLite from EBI FTP and filter to ExhalePath seed-gene targets — same activity schema feeds `train-chembl`.

### Single-cell Census harvest (top US diseases + body-wide healthy tissues)

ExhalePath pulls **CELLxGENE Census** primary human cells (~97M), keeps **composition summaries** (disease × tissue × cell type), leans into diseases/tissues with the most cells, and calibrates cell-state densities.

```bash
pip install cellxgene-census
python -m exhalepath harvest-census --top-n 100 --min-cells 5000
```

Outputs land in `data/census/` (`census_manifest.json`, healthy tissue tables, disease→ExhalePath state fractions). Predictions auto-use Census fractions when available.

### Data scale

`build-corpus` harvests clinical cases and mutation edges across dozens of GDC projects, expands pathway gene sets via Reactome, and materializes case×pathway×VOC training matrices that reach **hundreds of thousands to millions of rows** depending on project count.

VOC regression targets in the default corpus are **mechanistic-prior synthetic labels** (literature disease fold-changes + pathway scores + stage effects + noise) so the calibrator can be trained without a private breath biobank. Swap in measured breath cohorts when you have them — the feature schema is unchanged.

## Install

```bash
cd exhalepath
pip install -e ".[dev]"
```

## Quickstart — biomarker engine (recommended)

```bash
# Disease + any affected-cell location → top 50 VOCs with Δppb quantities
python -m exhalepath biomarker "lung adenocarcinoma" \
  --location lung --genes KRAS,TP53 --top 50 \
  --out-dir ./runs/biomarker_luad

python -m exhalepath biomarker "Alzheimer's" --location brain --top 50 \
  --out-dir ./runs/biomarker_ad

python -m exhalepath biomarker "type 2 diabetes" --location pancreas \
  --affected-fraction 0.4 --top 50

python -m exhalepath list-locations
python -m exhalepath list-diseases
```

### Python API (biomarker)

```python
from exhalepath import ExhaleBiomarkerEngine

engine = ExhaleBiomarkerEngine()
report = engine.predict(
    "breast invasive carcinoma",
    location="left breast upper outer",
    genes=["PIK3CA", "TP53"],
    top_n=50,
)
print(report.summary())
report.to_dataframe()  # rank, voc, healthy/predicted/Δ ppb, fold, physio traces
```

## Quickstart — legacy predict CLI

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
| `biomarker` | **Whole-body engine:** disease + location → top-N VOCs by \|Δppb\| |
| `explain` | WHY VOC changes: pathways, Census cells, driver genes |
| `build-mechanism-packs` | Mechanism packs for ~100 atlas diseases |
| `harvest-public-breath` | Scientific Data breathomics → public validation cases |
| `eval-public-breath` | Top-k / directional eval vs public breath sets |
| `predict` | VOC ppb panel for a disease / tumor context |
| `harvest-chembl` | ChEMBL activities for pathway genes + VOC physchem |
| `train-chembl` | Fit chemogenomic aux model on ChEMBL pChEMBL rows |
| `build-corpus` | Pull GDC + pathway data; write training matrices |
| `train` | Fit per-VOC calibrators |
| `list-diseases` | ~100-disease atlas |
| `list-locations` | Whole-body tissue map for affected-cell placement |
| `list-vocs` | 50-VOC catalog with healthy breath baselines |
| `harvest-census` | CELLxGENE Census compositions → cell-state calibration |
| `eval-multisite` | 20 diseases × 5 sites literature directional eval |
| `audit` | Literature accuracy, calibrator holdout, zero-shot + stress tests |

```bash
python -m exhalepath audit --out runs/audit/audit_report.json
```

## Knowledge bases

Bundled under `data/knowledge/`:

- `voc_catalog.json` — **50** endogenous breath VOCs with healthy ppb medians/ranges
- `whole_body_tissues.json` — full-body anatomic map (~66 tissues; Census healthy-cell counts)
- `pathway_voc_map.json` — metabolic pathways → VOC emission coefficients
- `voc_pathway_chains.json` — multi-step biosynthetic chains (ketogenesis→acetone, PUFA→aldehydes, …)
- `cell_state_atlas.json` — affected cell states (density × activity)
- `disease_voc_priors.json` — **~100** diseases with pathway biases + VOC log2FC priors
- `physio_constants.json` — λ blood:air, perfusion, hepatic first-pass, Farhi parameters

Rebuild atlas expansions with `python scripts/build_whole_body_atlas.py`.

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
