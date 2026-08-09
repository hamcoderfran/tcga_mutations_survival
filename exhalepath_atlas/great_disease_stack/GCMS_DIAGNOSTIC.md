# GC-MS patient diagnostic research layer

Companion to the Great Disease Stack. Turns mechanism VOC signatures into
**patient-level diagnostic research metrics** on public GC-MS cohorts.

```bash
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-patient-diagnostic --all --signature stack
voc eval-scidata-samples --cohort asthma          # per-sample Sci Data OVR
voc export-metabolights --study scidata:asthma    # mzTab-M + ISA-Tab scaffold
voc export-paper-pack runs/scidata_samples/asthma
voc lock-split --study ST000883
```

Outputs under `runs/patient_diagnostic/<STUDY>/` (and `runs/scidata_samples/<cohort>/`):

- `PATIENT_DIAGNOSTIC_REPORT.md` / `.json` (+ stratified AUROC when age/sex/smoking exist)
- `roc_curve.png`
- `METHODS.md` + `*_paper_pack.zip` (figures + Methods + overlay + split hash)
- `TRIPOD_AI_BREATHVOC_CHECKLIST.md`
- `*_split_manifest.json` (SHA256-locked folds)
- `filter_report.json` (blank-ratio / detection-fraction)
- `matrices/*_patient_voc_matrix.csv`

Sci Data caveat: labels are **one-vs-rest across pulmonary cohorts** (no healthy arm); smoking is absent in CBD metadata.

See [`RESEARCH.md`](../RESEARCH.md) for the scientific impact rationale.
