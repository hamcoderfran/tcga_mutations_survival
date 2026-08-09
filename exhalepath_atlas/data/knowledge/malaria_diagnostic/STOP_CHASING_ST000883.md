# Stop chasing ST000883 AUROC

ST000883 (n≈35) is a **methods harness**, not a clinical benchmark to hill-climb.

## Report these two numbers

1. **Transferable nested AUROC** (hybrid+lit signature): 58.3%
2. **Fit-on-cohort sparse ceiling** (nested SelectKBest): 73.3%

The gap is expected: signatures were not trained on these labels; sparse logistic was.

## What to do instead

- Add a second **intensity** cohort and run `--loso` (CSIRO labels are bundled; export peaks).
- Widen n (JBR: learning curves flatten near n≈50).
- Report fixed-sensitivity operating points + AUPRC CI95.
- Use ST003200 / Sci Data for confounder (age/sex/smoking) diligence — not malaria AUROC.

## Do not

- Retune atlas priors or remaps solely to raise ST000883 nested AUROC.
- Market sparse ceiling (~70%+) as transferable diagnostic performance.
- Claim FDA/clinical readiness from this n.
