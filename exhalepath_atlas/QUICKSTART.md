# Quick start — one line in, full visual out

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

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
