# Great Disease Prediction Stack

One comprehensive folder that fuses **20 models** across exhaled VOCs and the rest of disease biology — genetics, pathways, signaling, metabolic flux, ADME/PBPK transport, microbiome, cell states, pharmacology, comorbidity, literature priors, zero-shot evidence, phenotype/MONDO, null contrast, and a meta ensemble — on top of ExhalePath.

> Research / hypothesis-generation only. Not a medical device.

## One-line usage

```bash
pip install -e "exhalepath_atlas[dev,stack]"

# Structured
voc stack "depression" -l brain -c obesity --age 24 --sex male

# Naturalistic patient input → PatientTemplate → 20-model stack
voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations"
voc stack --nl "62F former smoker, stage II lung adenocarcinoma LLL, KRAS/TP53"
voc stack --nl "Maple syrup urine disease, genes BCKDHA BCKDHB"
```

`voc patient` accepts vignettes, clinic-note sections (`CC:`/`HPI:`/`PMH:`/`Meds:`),
`key: value` blocks, or JSON — then normalizes to a `PatientTemplate` the stack understands.
Use `--show-template` to inspect the parsed template without running models.

Writes `runs/stack_*/`:

| Artifact | Contents |
|---|---|
| `STACK_REPORT.html` | Full visual report |
| `stack_dashboard.png` | 4-panel consensus dashboard |
| `model_vote_heatmap.png` | Per-model VOC vote matrix |
| `fused_vocs.csv` | Consensus VOC table |
| `model_status.csv` | Which models fired |
| `STACK_RESULT.json` | Full machine-readable fusion |

## The 20 fused models

| # | Model ID | Aspect | Source family |
|---|---|---|---|
| 1 | `exhalepath_hybrid` | VOC quantity | ExhalePath hybrid physio+ML |
| 2 | `exhalepath_physiology` | VOC quantity | Farhi / cell-state physiology |
| 3 | `exhalepath_calibrator` | VOC quantity | Supervised calibrator residual |
| 4 | `zero_shot_mechanism` | Mechanism | Resolver + pathway→VOC projection |
| 5 | `zero_shot_evidence` | VOC prior | Literature themes / expected directions |
| 6 | `literature_voc_prior` | VOC prior | Curated disease→VOC atlas |
| 7 | `opentargets_genes` | Genetics | Open Targets associations |
| 8 | `humangem_flux` | Metabolic flux | Human-GEM-style flux proxy |
| 9 | `opera_physchem` | Transport / ADME | OPERA-style QSPR panel |
| 10 | `primekg_graph` | Knowledge graph | PrimeKG-lite multi-hop walks |
| 11 | `omnipath_signaling` | Signaling | OmniPath-style pathway diffusion |
| 12 | `agora_microbiome` | Microbiome | AGORA2/mVOC emission routes |
| 13 | `pbpk_transport` | Transport / ADME | Farhi/PBPK-lite alveolar delivery |
| 14 | `cell_census` | Cell state | CELLxGENE Census / atlas densities |
| 15 | `comorbidity_clinical` | Clinical | Comorbidity prior fusion |
| 16 | `chembl_pharm` | Pharmacology | ChEMBL pathway / physchem |
| 17 | `pathway_enrichment` | Pathways | KEGG/Reactome enrichment consensus |
| 18 | `phenotype_mondo` | Phenotype | Free-text / MONDO → VOC themes |
| 19 | `counterfactual_null` | Null | Deterministic null shrinkage |
| 20 | `meta_ensemble` | VOC quantity | Hybrid × literature meta head |

## Fusion method

1. Each model emits per-VOC `log2fc` + confidence (+ optional disease aspects).
2. **Mode-aware family multipliers** (atlas vs zero-shot regimes).
3. **Calibrated per-model weights** from `fusion_weights_calibrated.json`.
4. **Confidence-weighted average** of log2fc (effective weight × confidence).
5. **Anti-dilution anchor** toward hybrid/physio/calibrator when they agree.
6. **Reciprocal Rank Fusion (RRF)** across model rankings.
7. **Sign-agreement + epistemic std** → confidence and approximate 90% CI.
8. Aspects (genes, pathways, microbes, cell states, phenotypes) fused separately.

## Curated prior packs

Shipped under `exhalepath_atlas/data/great_stack/`:

- `model_registry.json` — model IDs + default fusion weights
- `fusion_weights_calibrated.json` — holdout-aware calibrated weights
- `primekg_lite_graph.json` — disease–gene–pathway–VOC graph
- `opera_adme_panel.json` — λ / exhalation efficiency / fu / Clint proxies
- `humangem_flux_proxy.json` — pathway capacity × VOC producers
- `agora_microbiome_routes.json` — microbial VOC routes by disease category
- `omnipath_signaling_lite.json` — pathway crosstalk (runtime-filled if empty)

## Optional heavier backends

```bash
pip install -e "exhalepath_atlas[stack,gsmm]"   # networkx + cobrapy
# Then place Human-GEM / Recon3D SBML under data/great_stack/gsmm/ for full FBA
# OPERA desktop app can refresh opera_adme_panel.json via scripts/refresh_opera.py
```

The stack **runs without** those installs — proxies use the bundled packs.

## Install / SOTA

- [INSTALL.md](INSTALL.md) — one-line pip + usage  
- [SOTA.md](SOTA.md) — honest cutting-edge vs clinical-SOTA assessment  
- Holdout: `voc eval-stack-holdout` → `STACK_HOLDOUT_REPORT.md`

## Layout

```
great_disease_stack/          ← you are here (docs + entry)
  README.md
  INSTALL.md
  SOTA.md
  docs/
  examples/
../src/exhalepath/great_stack/  ← engine, models, fusion, report
../data/great_stack/            ← curated multi-model priors
```

## Example

```bash
voc stack "Maple syrup urine disease" -l systemic --genes BCKDHA,BCKDHB,DBT
voc stack "lung adenocarcinoma" -l lung -c COPD --genes KRAS,TP53 --age 62 --sex male --smoking former
```
