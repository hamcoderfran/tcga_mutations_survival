# Patient-level GC-MS diagnostic research report

**Generated:** 2026-08-02T02:19:45.438659+00:00
**Study:** `ST000587` · disease `heart_failure`
**Signature source:** hybrid · method `cosine`

> Research enablement report — not a clinical validation claim.

## Primary metrics (full-cohort signature match)

- n subjects: **23** (pos=11, neg=12)
- AUROC: **55.3%**  CI95=(0.27852855477855476, 0.7950689935064935)
- AUPRC: **61.1%**
- Sensitivity / Specificity (Youden): **36.4%** / **91.7%**
- PPV / NPV: **80.0%** / **61.1%**
- Confusion: `{'tn': 11, 'fp': 1, 'fn': 7, 'tp': 4}`
- Brier (rank-scaled): 0.305119800337177

## Nested / locked-split evaluation

- Strategy: `stratified_kfold` · seed=42 · sha256=`c56043631bf59929…`
- Mean test AUROC across folds: **60.0%**
- Non-nested (full-fit) AUROC: **55.3%**
- Optimism gap (non-nested − nested): **-4.7%**

### Per-fold

| Fold | n_test | AUROC |
|---|---:|---:|
| fold_0 | 5 | 100.0% |
| fold_1 | 5 | 16.7% |
| fold_2 | 5 | 33.3% |
| fold_3 | 4 | 75.0% |
| fold_4 | 4 | 75.0% |

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
