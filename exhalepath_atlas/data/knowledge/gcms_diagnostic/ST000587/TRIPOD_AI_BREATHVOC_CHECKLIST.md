# TRIPOD+AI / BreathVOC-1.0 lightweight checklist

Auto-filled fields from `voc eval-patient-diagnostic`. Complete remaining items before publication.

| Item | Status / value |
|---|---|
| Title identifies study aim | Patient-level GC-MS VOC diagnostic research enablement |
| Study / data source | ST000587 (heart_failure) |
| Eligibility / labels | Disease vs control from MW factors |
| Predictors | Mapped atlas VOC intensities; signature=hybrid |
| Outcome | Binary infection/disease status |
| Sample size | n=23 |
| Missing data handling | Column median fill for sparse mapped VOCs (documented) |
| Internal validation | Locked stratified_kfold · sha256=c56043631bf59929dfb037a3… |
| Nested vs non-nested | nested=0.6 non-nested=0.5530303030303031 gap=-0.046969696969696884 |
| Discrimination | AUROC + AUPRC reported |
| Calibration | Brier on rank-scaled scores (proxy) |
| Fairness / confounders | Smoking/age often unavailable in MW ST000883 — flagged as gap |
| Open science | Code path: `exhalepath.gcms` / `voc eval-patient-diagnostic` |
| Clinical use claim | **None** — research hypothesis / enablement only |

Generated: 2026-08-02T02:19:45.559146+00:00
