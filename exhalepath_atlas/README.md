# ExhalePath Atlas

**Fully integrated exhaled VOC biomarker platform** — priority datasources **1–14**, whole-body physiology, and mechanism explainability across **100+ diseases**.

This folder is the GitHub deliverable that packages the complete stack (prediction + validation + all external data integrations).

> Research / hypothesis-generation tool. Not a medical device.

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
python -m exhalepath biomarker "lung adenocarcinoma" \
  --location lung --genes KRAS,TP53 --top 50 \
  --out-dir runs/atlas_luad

python -m exhalepath explain "type 2 diabetes" --location pancreas --top 10
python -m exhalepath build-mechanism-packs --top-vocs 12

python -m exhalepath harvest-public-breath
python -m exhalepath eval-public-breath --top-k 15
```

## Architecture

```
Disease + location (+ genes)
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
