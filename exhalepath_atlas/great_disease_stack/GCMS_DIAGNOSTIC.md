# GC-MS patient diagnostic research layer

Companion to the Great Disease Stack. Turns mechanism VOC signatures into
**patient-level diagnostic research metrics** on public GC-MS cohorts.

```bash
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-patient-diagnostic --all --signature stack
voc lock-split --study ST000883
```

Outputs under `runs/patient_diagnostic/<STUDY>/`:

- `PATIENT_DIAGNOSTIC_REPORT.md` / `.json`
- `roc_curve.png`
- `TRIPOD_AI_BREATHVOC_CHECKLIST.md`
- `*_split_manifest.json` (SHA256-locked folds)
- `matrices/*_patient_voc_matrix.csv`

See `/exhalepath_atlas/RESEARCH.md` for the scientific impact rationale.
