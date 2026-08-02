# Great Disease Stack — Holdout & Break Report

## SOTA / cutting-edge verdict

**Clinical breathomics SOTA?** `False`
**Cutting-edge open systems stack?** `True`

Cutting-edge as an open, multi-model *systems* breath-hypothesis stack (physiology + GSMM-style flux + KG + ADME + microbiome + zero-shot + fusion). Not SOTA as a validated clinical diagnostic vs GC-MS/PTR trial models (e.g. sensor-array or cohort-trained classifiers reporting ROC-AUC on held-out patients).

### Novel (open)

- 16-model fusion across VOC quantity + genetics + flux + ADME + microbiome + signaling
- Naturalistic PatientTemplate → stack path for diverse clinical text
- Mechanism explainability + uncertainty/agreement reporting
- Installable open package with offline curated priors

### Gaps to clinical SOTA

- No prospective multi-site GC-MS/PTR patient-level holdout with locked labels
- Public Sci Data panels are cross-cohort differentials, not healthy-controlled absolute ppb
- Literature/priority10 directional labels partially overlap atlas priors (circularity risk)
- Human-GEM/OPERA/PrimeKG are proxy packs, not full downloaded GEMs/OPERA binaries
- No head-to-head vs published breath ML baselines on identical splits

## Holdout directional scores (stack vs hybrid baseline)

| Benchmark | Hybrid dir | Stack dir | Stack recall@k |
|---|---:|---:|---:|
| Public breath | 98.2% | 100.0% | 95.2% |
| Literature | 100.0% | 100.0% | 97.9% |
| Priority-10 | 100.0% | 97.8% | — |

## PatientTemplate adversarial suite

- Hard cases: **13/13** (100.0%)
- Soft/expected-fragile: **5/5** (100.0%)

## Naturalistic → stack e2e

- Cases: 3 · mean models ok: 15.666666666666666

- `schizophrenia` · models 16/16 · top=pentane · 0.382s
- `depression` · models 16/16 · top=pentane · 0.292s
- `lung adenocarcinoma` · models 15/16 · top=pentane · 0.298s

## Break case details

- [PASS] `compact_ok` → disease='schizophrenia' conf=0.8500000000000001 fails=[]
- [PASS] `json_ok` → disease='asthma' conf=0.8 fails=[]
- [PASS] `kv_ok` → disease='type 2 diabetes' conf=0.75 fails=[]
- [PASS] `note_ok` → disease='copd' conf=0.9000000000000001 fails=[]
- [PASS] `empty` → disease=None conf=0.0 fails=[]
- [PASS] `whitespace` → disease=None conf=0.0 fails=[]
- [PASS] `no_disease_demographics_only` → disease=None conf=0.35 fails=[]
- [PASS] `gene_only_trap` → disease=None conf=0.2 fails=[]
- [PASS] `multi_disease_primary` → disease='depression' conf=0.65 fails=[]
- [PASS soft] `typo_disease` → disease='schizophrenia' conf=0.75 fails=[]
- [PASS] `emoji_noise` → disease='asthma' conf=0.8500000000000001 fails=[]
- [PASS soft] `spanish_mix` → disease='Type 2 diabetes mellitus' conf=0.65 fails=[]
- [PASS soft] `contradictory_sex` → disease='depression' conf=0.75 fails=[]
- [PASS] `huge_note` → disease='pancreatic adenocarcinoma' conf=0.8 fails=[]
- [PASS soft] `injectiony` → disease='copd' conf=0.6 fails=[]
- [PASS soft] `bare_cancer` → disease=None conf=0.15 fails=[]
- [PASS] `medication_as_disease_trap` → disease=None conf=0.2 fails=[]
- [PASS] `luad_lll_genes` → disease='lung adenocarcinoma' conf=0.65 fails=[]

> Research / hypothesis-generation. Not a medical device. Directional scores can be partly circular with atlas priors.

