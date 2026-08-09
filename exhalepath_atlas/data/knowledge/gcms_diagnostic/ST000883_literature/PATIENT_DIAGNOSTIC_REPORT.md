# Patient-level GC-MS diagnostic research report

**Generated:** 2026-08-02T02:19:57.086160+00:00
**Study:** `ST000883` · disease `malaria`
**Signature source:** literature · method `cosine`

> Research enablement report — not a clinical validation claim.

## Primary metrics (full-cohort signature match)

- n subjects: **35** (pos=17, neg=18)
- AUROC: **66.0%**  CI95=(0.46598639455782315, 0.8371376811594201)
- AUPRC: **70.5%**
- Sensitivity / Specificity (Youden): **82.4%** / **55.6%**
- PPV / NPV: **63.6%** / **76.9%**
- Confusion: `{'tn': 10, 'fp': 8, 'fn': 3, 'tp': 14}`
- Brier (rank-scaled): 0.28310717470580515

## Nested / locked-split evaluation

- Strategy: `stratified_kfold` · seed=42 · sha256=`e1b4271ace3b6d55…`
- Mean test AUROC across folds: **63.3%**
- Non-nested (full-fit) AUROC: **66.0%**
- Optimism gap (non-nested − nested): **2.7%**

### Per-fold

| Fold | n_test | AUROC |
|---|---:|---:|
| fold_0 | 7 | 58.3% |
| fold_1 | 7 | 100.0% |
| fold_2 | 7 | 33.3% |
| fold_3 | 7 | 50.0% |
| fold_4 | 7 | 75.0% |

## Why this matters for science

- Patient-level holdout (not peak-level) — the unit of clinical inference
- Locked split hash supports preregistration-style reporting
- Mechanism signature (ExhalePath/stack) tested as a transferable template
- Optimism gap surfaces leakage / overfit risk for papers

## Caveats

- Mapped atlas VOC subset only — many GC-MS peaks are unmapped and dropped.
- ST000883 lacks smoking/age in factors.json — confounder-stratified AUCs unavailable.
- Mechanism signature is independent of these patient labels, but VOC name mapping can still introduce circularity with literature priors.
- Small n (≈35) → wide bootstrap CIs; treat AUROC as research enablement evidence, not clinical validation.
- Not a medical device. No clinical diagnostic claim.
