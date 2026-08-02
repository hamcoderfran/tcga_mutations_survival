# Great Disease Prediction Stack

One comprehensive folder that fuses **16 models** across exhaled VOCs and the rest of disease biology — genetics, pathways, signaling, metabolic flux, ADME/PBPK transport, microbiome, cell states, pharmacology, comorbidity, and literature priors — on top of ExhalePath.

> Research / hypothesis-generation only. Not a medical device.

## One-line usage

```bash
pip install -e "exhalepath_atlas[dev,stack]"

# Structured
voc stack "depression" -l brain -c obesity --age 24 --sex male

# Naturalistic patient input → PatientTemplate → 16-model stack
voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations"
voc stack --nl "62F former smoker, stage II lung adenocarcinoma LLL, KRAS/TP53"
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

## The 16 fused models

| # | Model ID | Aspect | Source family |
|---|---|---|---|
| 1 | `exhalepath_hybrid` | VOC quantity | ExhalePath hybrid physio+ML |
| 2 | `exhalepath_physiology` | VOC quantity | Farhi / cell-state physiology |
| 3 | `exhalepath_calibrator` | VOC quantity | Supervised calibrator residual |
| 4 | `zero_shot_mechanism` | Mechanism | Custom resolver + ontology NN |
| 5 | `literature_voc_prior` | VOC prior | Curated disease→VOC atlas |
| 6 | `opentargets_genes` | Genetics | Open Targets associations |
| 7 | `humangem_flux` | Metabolic flux | Human-GEM-style flux proxy |
| 8 | `opera_physchem` | Transport / ADME | OPERA-style QSPR panel |
| 9 | `primekg_graph` | Knowledge graph | PrimeKG-lite multi-hop walks |
| 10 | `omnipath_signaling` | Signaling | OmniPath-style pathway diffusion |
| 11 | `agora_microbiome` | Microbiome | AGORA2/mVOC emission routes |
| 12 | `pbpk_transport` | Transport / ADME | Farhi/PBPK-lite alveolar delivery |
| 13 | `cell_census` | Cell state | CELLxGENE Census / atlas densities |
| 14 | `comorbidity_clinical` | Clinical | Comorbidity prior fusion |
| 15 | `chembl_pharm` | Pharmacology | ChEMBL pathway / physchem |
| 16 | `pathway_enrichment` | Pathways | KEGG/Reactome enrichment consensus |

## Fusion method

1. Each model emits per-VOC `log2fc` + confidence (+ optional disease aspects).
2. **Confidence-weighted average** of log2fc (weight × confidence).
3. **Reciprocal Rank Fusion (RRF)** across model rankings.
4. **Sign-agreement** shrinks confidence when models disagree on direction.
5. Aspects (genes, pathways, microbes, cell states, physiology) are fused separately into a multi-scale disease profile.

## Curated prior packs

Shipped under `exhalepath_atlas/data/great_stack/`:

- `model_registry.json` — model IDs + default fusion weights
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
