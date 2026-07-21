# ExhalePath Atlas (`voc`)

**Install once, predict exhaled VOC biomarkers:**

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
voc "depression" -l brain -c obesity --age 24 --sex male
```

See [QUICKSTART.md](QUICKSTART.md). Package name on pip: **`voc-breath`** · command: **`voc`**.

**Fully integrated exhaled VOC biomarker platform** — priority datasources **1–14**, comorbidities, whole-body physiology, and mechanism explainability across **100+ diseases**.

This folder is the GitHub deliverable that packages the complete stack (prediction + validation + all external data integrations).

> Research / hypothesis-generation tool. Not a medical device.

## Install

```bash
# From GitHub (recommended)
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"

# From this folder
cd exhalepath_atlas
pip install -e ".[dev]"

# One-liner script
bash install_voc.sh
```

After install you get the **`voc`** command (aliases: `exhalepath`, `voc-breath`).

## Real-data training (no synthetic VOC labels)

See [REAL_DATA.md](REAL_DATA.md). Quick path:

```bash
voc build-real-corpus    # MW + Sci Data + literature panels only
voc train                # exhalepath-calibrator-real-1.0
voc eval-priority10      # asthma…ARDS directional check
```

Priority supervised diseases: asthma, COPD, COVID-19, pneumonia/VAP/HAP, TB, CF, OSA, gastric/H&N/prostate cancer, heart failure, malaria, ARDS.

## What is integrated

| Pri | Datasource | Role |
|---|---|---|
| 1 | Metabolomics Workbench / MetaboLights / MassIVE | Quantified breath studies + healthy ppbv baselines |
| 2 | HMDB (+ PubChem bridge) | VOC identity / biofluid annotation for full panel |
| 3 | Blood:air λ partition table | Farhi alveolar physics |
| 4 | mVOC microbial emitters | Dysbiosis → breath VOC routes |
| 5 | PubChem physchem | MW / XLogP / TPSA for VOC panel |
| 6 | Reactome expansion | Pathway gene neighborhoods |
| 7 | GTEx enzyme×tissue priors | Where VOC enzymes are expressed |
| 8 | Open Targets / GWAS | Disease→gene priors for 100+ diseases |
| 9 | TCGA/GDC drivers | Cancer mutation panels |
| 10 | BindingDB pharmacology notes | Ligand-class VOC effects |
| 11 | Blood→breath proxy links | Systemic metabolite bridges (UKB-style, public) |
| 12 | NIST/GC-MS RI metadata | Peak→VOC identification aids |
| 13 | HBDB / EPA VOLATILOME | Full 777-compound breath catalog + disease links |
| 14 | KEGG | VOC → reaction → pathway → enzyme maps |
| + | ChEMBL / CELLxGENE Census / public breathomics | Bundled from ExhalePath core |

## Install

```bash
cd exhalepath_atlas
pip install -e ".[dev]"
```

## One-command datasource integration

```bash
# Live pulls where APIs allow + curated tables for the rest
python -m exhalepath integrate-datasources

# Offline / CI-safe curated fusion
python -m exhalepath integrate-datasources --offline
```

Artifacts land in:
- `data/datasources/<key>/` — raw harvests
- `data/knowledge/datasource_*.json` — fused knowledge fragments
- `data/datasources/INTEGRATION_MANIFEST.json` — full audit trail
- `data/knowledge/atlas_capability.json` — disease × datasource capability card

## Biomarker + WHY (100+ diseases)

```bash
# Short form (recommended) — disease is the first argument
voc "lung adenocarcinoma" -l lung --genes KRAS,TP53 --top 50 --out-dir runs/atlas_luad

# Comorbidities fuse pathway bias, VOC priors, and cell-state modulation
voc "depression" -l brain --age 24 --sex male -c obesity --top 20
voc "schizophrenia" -l brain --age 18 --sex male -c heart_disease --top 20

# Natural language OR interactive question fields
voc ask                                          # prompt: disease, comorbidities, location, …
voc nl "24yo obese male with depression"         # optional light LLM → slots → predict
voc nl "stage II LUAD left lower lobe KRAS TP53" --llm rules --yes
voc ask --nl "depression with obesity" --llm ollama   # tiny local model if Ollama is up

# Creutzfeldt–Jakob clinical VOC tempo (incubating → terminal; research hypothesis)
voc cjd-profile --age 62 --sex female
voc "Creutzfeldt-Jakob disease" -l brain --genes PRNP --age 62 --sex female

# Explicit subcommands still work
voc biomarker "type 2 diabetes" -l pancreas --top 20
voc explain "type 2 diabetes" --location pancreas --top 10
voc list-diseases
voc list-locations
voc list-vocs

voc harvest-public-breath
voc eval-public-breath --top-k 15
voc harvest-clinical-comorbidity
voc eval-comorbidity-clinical --top-k 15
```

### Natural language / ask mode

| Mode | Command | Notes |
|---|---|---|
| Questionnaire | `voc ask` | Impute disease, comorbidities, location, age, sex, stage, genes, smoking, mode |
| NL + optional LLM | `voc nl "…"` | `--llm auto\|rules\|ollama\|openai`; falls back to rules, then ask if disease missing |
| Seeded ask | `voc ask --nl "…"` | Parse NL first, then edit fields interactively |

Env: `VOC_LLM`, `VOC_OLLAMA_HOST`, `VOC_OLLAMA_MODEL` (default `qwen2.5:0.5b`), `OPENAI_API_KEY` / `VOC_OPENAI_API_KEY`.

## Architecture

```
Disease + comorbidities + location (+ genes)
   ├─ Fuse comorbidity pathway_bias / VOC priors / cell-state mults
   ├─ Open Targets / GDC drivers (pri 8–9)
   ├─ Census cell populations (CELLxGENE)
   ├─ Pathway scores + Reactome genes (pri 6)
   ├─ mVOC microbial routes (pri 4)
   ├─ GTEx enzyme localization (pri 7)
   ├─ Cell→blood→alveolar (λ from pri 3, PubChem pri 5)
   └─ Top-50 VOC Δppb + mechanism WHY text
```

## Tests & dataset completion gate

```bash
pytest -q

# Full dataset suite (pytest + audit + multisite + public breath) — must pass
python -m exhalepath eval-completion --out-dir runs/completion
```

Latest completion status is recorded in `data/knowledge/completion_status.json`.
