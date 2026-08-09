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
| Sci Data **per-sample** peak tables | `voc eval-scidata-samples` / `load_scidata_ovr_matrix` | Asthma/COPD/bronchiectasis sample×VOC (not cohort means) |
| Locked patient splits | `voc lock-split --study ST000883` | Preregistration-style SHA256 manifests |
| Signature diagnostic scoring | `voc eval-patient-diagnostic` | Tests hybrid/stack/literature templates with AUROC/AUPRC/sens/spec |
| Age/sex/smoking stratified AUCs | auto in diagnostic reports when metadata exist | Confounder-aware reporting (smoking N/A on Sci Data 2024) |
| Blank-ratio / detection filters | `blank_ratio_filter` on peak tables | On-breath vs background when blanks present |
| MetaboLights scaffold | `voc export-metabolights` | mzTab-M-like + ISA-Tab investigation/study/assay |
| Observed vector scoring | `voc score-sample … --vocs '{…}'` | Lab peak tables → disease template match |
| One-zip paper pack | `voc export-paper-pack` | figures + Methods + overlay + split hash |

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

# 3) Sci Data per-sample (one-vs-rest; age/sex strata; paper zip)
voc eval-scidata-samples --cohort asthma
voc eval-scidata-samples --all

# 4) Deposit scaffolds + re-zip a run
voc export-metabolights --study ST000883
voc export-metabolights --study scidata:asthma
voc export-paper-pack runs/patient_diagnostic/ST000883

# 5) Score a new patient's mapped VOC vector
voc score-sample malaria --signature stack \
  --vocs '{"acetone":0.4,"pentane":0.5,"isoprene":-0.3,"benzene":0.6}'
```

**Recommended paper workflow**

1. Preregister / cite `content_sha256` from the split manifest  
2. Report **nested** AUROC (and optimism gap vs non-nested)  
3. Report Youden sens/spec + confusion + AUPRC  
4. State MSI / mapping limits (atlas VOC subset only)  
5. Report stratified AUCs when age/sex/smoking metadata exist (Sci Data: age/sex yes, smoking no)  
6. Attach the one-zip paper pack + MetaboLights scaffold if depositing  
7. Claim *research enablement / transferable signature test*, not clinical SOTA  

## Public datasets (adapters)

| Tier | Dataset | Status |
|---|---|---|
| A | Sci Data 2024 clinical breathomics (Figshare) | **Per-sample adapter shipped** (`eval-scidata-samples`); OVR only (no healthy arm) |
| A | ST000883 / ST000587 MW patient matrices | **Shipped** (`eval-patient-diagnostic`) |
| A | ST003200 healthy PTR (n=504) | Confounder effect-size reference (age/sex/smoking) — next |
| A | MSV000095340 pediatric asthma GC-qTOF (mzML) | Raw-spectrum pipeline demos — next |
| B | RADicA / ReCIVA blank-aware tables | Blank-first schema gold standard — next |
| B | Owlstone OMNI example | Industry feature-table shape reference — next |

## Draft branches — cherry-picked vs pruned

Do **not** merge wholesale draft PRs #13–#19 into `main`. This branch cherry-picked additive modules only:

| Cherry-picked | Skipped (regressive) |
|---|---|
| `secure_fetch` + HTTP allowlist (#14) | Wholesale `cli.py` / `predict.py` tip rewrites |
| `eval/coverage_audit.py` + COVERAGE_AUDIT artifacts (#14) | `disease_voc_priors.json` replacements |
| `eval/lit_compare.py` + lit_compare artifacts (#13) | Recalibrated `voc_calibrator.joblib` |
| Surgical post-blend smoking/age exo fix in `predict.py` | Full model-improve-100disease merge (#15) |
| `ds15_alt_breath_sources` + catalogs (#17) | Prior/CLI conflicts that strip zero-shot hardening |
| HBDB 60-disease JSON + extract script (#16) | |
| HMDB Wishart mirror constant + breath extracts (#18) | |

Industry diligence pack: `voc eval-industry-pack` → `data/knowledge/industry_pack/` (`INDUSTRY_PACK.md`, `BUYER_BRIEF.md`).

## Malaria metric upgrade (`voc eval-malaria-diagnostic`)

Literature-guided unlocks on ST000883 (Schaber 2018 / Berna 2015):

1. **Peak remaps** — α-pinene, 3-carene, cyclohexanone, tridecane, methyl-alkanes (+ thioether aliases)  
2. **Nested sparse** (inner FS) vs **transferable hybrid+lit signature**  
3. **Fixed-sensitivity** operating points with bootstrap spec CI95  
4. **External catalog + learning curve** — CSIRO CHMI / JID 2024 noted; open intensity table still only ST000883 (n=35)

Latest archived report: `data/knowledge/malaria_diagnostic/MALARIA_DIAGNOSTIC_UPGRADE.md`.

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
