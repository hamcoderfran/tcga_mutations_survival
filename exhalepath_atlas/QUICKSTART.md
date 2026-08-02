# Quick start — one line in, full visual out

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

## Great Disease Stack

```bash
voc stack "depression" -l brain -c obesity --age 24 --sex male
# → STACK_REPORT.html + multi-model consensus
```

## GC-MS patient diagnostic research

```bash
voc eval-patient-diagnostic --study ST000883          # AUROC + locked splits + paper pack
voc eval-patient-diagnostic --benchmark               # hybrid/stack/literature comparison
voc lock-split --study ST000883                       # preregistration-style SHA256 folds
```

See [RESEARCH.md](RESEARCH.md) and [great_disease_stack/README.md](great_disease_stack/README.md).

## One-line prediction (recommended)

```bash
voc "depression" -l brain -c obesity --age 24 --sex male
```

That single command:
1. Resolves the disease (atlas or zero-shot mechanism transfer)
2. Prints a rich terminal panel (top VOCs + pathways + why)
3. **Auto-writes a full visual pack** under `runs/voc_<disease>_…/`:
   - `REPORT.html` — open in a browser
   - `dashboard.png` — 4-panel overview
   - CSV/JSON of VOCs, pathways, cell states, mechanisms

```bash
# Novel / unseen disease
voc "Maple syrup urine disease" -l systemic
voc predict-novel "unknown phenotype" -l liver -d "ketotic mitochondrial stress" --genes HMGCS2,CPT1A

# Natural language
voc nl "24yo obese male with depression" --yes

# Skip saving artifacts
voc "COPD" -l lung --no-save
```

## Readiness / eval (optional)

```bash
voc eval-implementation-readiness
```

Research / hypothesis tool only — not a medical device.
