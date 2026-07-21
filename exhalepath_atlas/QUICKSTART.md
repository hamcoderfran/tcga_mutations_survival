# Quick start — install & run in 30 seconds

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

Or from a downloaded wheel:

```bash
pip install voc_breath-*.whl
```

Then:

```bash
voc "depression" -l brain -c obesity --age 24 --sex male
voc "schizophrenia" -l brain -c heart_disease --age 18 --sex male
voc "lung adenocarcinoma" -l "left lower lobe" --stage II --genes KRAS,TP53
voc list-diseases
voc --help
```

No subcommand needed for predictions — `voc "<disease>" …` runs the biomarker engine.

Research / hypothesis tool only — not a medical device.
