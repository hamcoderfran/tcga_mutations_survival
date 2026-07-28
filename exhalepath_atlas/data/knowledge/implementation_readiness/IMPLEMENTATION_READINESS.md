# Implementation readiness

**Ready:** True
**Score:** 12/12 gates (100%)
**Recommendation:** READY for research / hypothesis-generation deployment with uncertainty flags and mechanism audit trails.

## Gates
- ✓ `pytest`
- ✓ `audit`
- ✓ `calibrator_holdout`
- ✓ `priority10`
- ✓ `vision_fidelity_ge_90`
- ✓ `multisite_ge_0.90`
- ✓ `public_breath`
- ✓ `comorbidity`
- ✓ `stress_hard_ge_95`
- ✓ `zero_shot_hard_ge_85`
- ✓ `zero_shot_reliability`
- ✓ `cohort_ok_ge_99`

## Layer metrics
- **pytest**: passed=True, seconds=56.4
- **audit**: passed=True, seconds=10.37, calibrator_mae=0.16131710737762303 skipped=False
- **priority10**: passed=True, mean_directional_accuracy=1.0, seconds=0.85
- **vision**: passed=True, vision_fidelity_pct=98.25, n_profiles=50, seconds=6.17
- **multisite**: passed=True, composite_accuracy=0.9875, seconds=6.71
- **public_breath**: passed=True, mean_directional_accuracy=0.9821428571428571, seconds=0.92
- **comorbidity**: passed=True, mean_directional_accuracy=1.0, n_cases=4, seconds=0.29
- **stress_hard**: passed=True, pass_pct=100.0, n_cases=50, seconds=3.28
- **zero_shot_hard**: passed=True, pass_rate=1.0, n_profiles=58, failed_ids=[]
- **zero_shot_reliability**: passed=True, seconds=8.34
- **cohort**: passed=True, ok_pct=100.0, n_patients=1000, seconds=139.43

## Limitations (read before shipping)
- Zero-shot / hard-break suites score mechanism consistency and adversarial robustness — not external GC-MS clinical accuracy.
- Public scientific-data elevated recall@15 remains modest; directional agreement is the stronger public-breath signal today.
- Calibrator is fit on a small real corpus; treat MAE as small-n diagnostics, not a large held-out clinical trial.
- Owlstone and other gated breath atlases are not used for training.
