# Patient-level GC-MS diagnostic research report

**Generated:** 2026-08-02T02:19:44.752455+00:00
**Study:** `ST000883` · disease `malaria`
**Signature source:** hybrid · method `cosine`

> Research enablement report — not a clinical validation claim.

## Primary metrics (full-cohort signature match)

- n subjects: **35** (pos=17, neg=18)
- AUROC: **56.9%**  CI95=(0.37494164332399627, 0.7612412587412585)
- AUPRC: **53.1%**
- Sensitivity / Specificity (Youden): **94.1%** / **33.3%**
- PPV / NPV: **57.1%** / **85.7%**
- Confusion: `{'tn': 6, 'fp': 12, 'fn': 1, 'tp': 16}`
- Brier (rank-scaled): 0.33271601583490906

## Nested / locked-split evaluation

- Strategy: `stratified_kfold` · seed=42 · sha256=`e1b4271ace3b6d55…`
- Mean test AUROC across folds: **53.3%**
- Non-nested (full-fit) AUROC: **56.9%**
- Optimism gap (non-nested − nested): **3.5%**

### Per-fold

| Fold | n_test | AUROC |
|---|---:|---:|
| fold_0 | 7 | 66.7% |
| fold_1 | 7 | 58.3% |
| fold_2 | 7 | 33.3% |
| fold_3 | 7 | 50.0% |
| fold_4 | 7 | 58.3% |

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
