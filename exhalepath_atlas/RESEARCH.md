# Research impact: GC-MS diagnostic enablement

**Goal:** Make ExhalePath Atlas / `voc-breath` the strongest *open research enablement* stack for exhaled-VOC / GC-MS disease studies — not by claiming clinical diagnostic SOTA, but by giving labs the scarce public goods that actually move papers and science forward.

## What creates the greatest scientific impact (2020–2026 evidence)

Ranked levers (see also `great_disease_stack/SOTA.md`):

1. **Patient-level nested / locked splits** — peak-level CV is the #1 overoptimism source in breath ML  
2. **Breath-native metadata** — device, fraction, blanks, batch, smoking/age/sex  
3. **On-breath vs background** — blank ratios + ID confidence (MSI levels)  
4. **Confounder-aware metrics** — disease AUC next to smoking/site AUC  
5. **Leakage-safe batch/QC correction** — fit on train/QC only  
6. **Public dataset adapters** — one schema for MW / Sci Data / MassIVE / OMNI-style tables  
7. **Identification evidence trail** — RI + library score + MSI level, not just a name  
8. **Paper-ready reporting** — TRIPOD+AI / STARD-oriented packs with split hashes  

This release implements the highest-leverage missing piece for *this* repo: **patient-level GC-MS matrices + AUROC diagnostic harness + locked splits + paper pack**, wired to the existing mechanism stack.

## What we shipped (`exhalepath.gcms` + CLI)

| Capability | Command / API | Why it matters |
|---|---|---|
| Patient×VOC matrix from MW | `load_mw_patient_matrix("ST000883")` | Restores per-patient intensities (not only group medians) |
| Locked patient splits | `voc lock-split --study ST000883` | Preregistration-style SHA256 manifests |
| Signature diagnostic scoring | `voc eval-patient-diagnostic` | Tests hybrid/stack/literature templates with AUROC/AUPRC/sens/spec |
| Observed vector scoring | `voc score-sample … --vocs '{…}'` | Lab peak tables → disease template match |
| Paper pack | `PATIENT_DIAGNOSTIC_REPORT.md` + ROC + TRIPOD checklist stub | Citation-grade methods scaffolding |

Bundled studies today:

- **ST000883** — malaria breath GC-MS (n≈35, Positive/Negative)  
- **ST000587** — heart-failure EBC (n≈23, HF/control)  

## How labs should use this (greatest impact path)

```bash
pip install -e "exhalepath_atlas[dev,stack]"

# 1) Lock a patient-level split (cite the sha256 in the paper)
voc lock-split --study ST000883 --seed 42 \
  --out data/knowledge/gcms_splits/ST000883_split_v1.json

# 2) Run patient diagnostic eval (mechanism signature + logistic baseline)
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-patient-diagnostic --all --signature stack

# 3) Score a new patient's mapped VOC vector
voc score-sample malaria --signature stack \
  --vocs '{"acetone":0.4,"pentane":0.5,"isoprene":-0.3,"benzene":0.6}'
```

**Recommended paper workflow**

1. Preregister / cite `content_sha256` from the split manifest  
2. Report **nested** AUROC (and optimism gap vs non-nested)  
3. Report Youden sens/spec + confusion + AUPRC  
4. State MSI / mapping limits (atlas VOC subset only)  
5. If smoking/age available, add stratified AUCs (schema fields ready; ST000883 factors lack them)  
6. Claim *research enablement / transferable signature test*, not clinical SOTA  

## Public datasets to prioritize next (adapters)

| Tier | Dataset | Why |
|---|---|---|
| A | Sci Data 2024 clinical breathomics (Figshare 23522490) | Asthma/COPD/bronchiectasis peak tables; ML baseline community |
| A | ST003200 healthy PTR (n=504) | Confounder effect-size reference (age/sex/smoking) |
| A | MSV000095340 pediatric asthma GC-qTOF (mzML) | Raw-spectrum pipeline demos |
| B | RADicA / ReCIVA blank-aware tables | Blank-first schema gold standard |
| B | Owlstone OMNI example | Industry feature-table shape reference |

## Current results on bundled patient GC-MS (regenerate anytime)

`voc eval-patient-diagnostic --benchmark --study ST000883`

| Signature | Method | AUROC | Nested AUROC | Nested logistic (same features) |
|---|---|---:|---:|---:|
| hybrid | cosine | ~57% | ~53% | ~73% |
| stack | cosine | ~58% | ~58% | ~73% |
| literature | cosine | ~66% | ~63% | ~73% |
| literature | dot | ~70% | ~65% | ~73% |

**How to read this for science**

- Nested logistic (~73%) is a **same-feature ceiling** when labels are fit on ST000883 itself.  
- Mechanism hybrid/stack signatures are **weaker but transferable** — they were not trained on these patient labels.  
- Literature panels score higher here partly because they encode published malaria breath directions (circularity risk).  
- Wide bootstrap CIs (n=35) mean these numbers are for **methods development**, not clinical claims.  
- The impact is the harness: locked splits, optimism gaps, paper packs other labs can reuse on larger cohorts.

## Honest boundary

- This stack **enables** GC-MS diagnostic *research* (design, eval, reporting).  
- It does **not** replace prospective multi-site clinical validation.  
- Directional literature panels remain partly circular with atlas priors; patient-level MW matrices are the cleaner stress test.

## Key references

- Kuo et al., *Sci Data* (2024): clinical breathomics dataset — doi:10.1038/s41597-024-03052-2  
- Nested vs non-nested breath ML optimism — JBR 2024 doi:10.1088/1752-7163/ad2b6e  
- Peppermint consortium sampling — doi:10.1088/1752-7163/aba130  
- TRIPOD+AI — BMJ 2024;385:e078378  
- ST003200 healthy breath PTR — doi:10.1007/s11306-024-02139-6  
- ST000883 malaria breath — Metabolomics Workbench; related doi:10.1093/infdis/jiy072  

See also: `great_disease_stack/SOTA.md`, `REAL_DATA.md`, `data/knowledge/gcms_diagnostic/`.
