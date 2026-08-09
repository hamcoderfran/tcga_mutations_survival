# Install & use ExhalePath Great Disease Stack

Research / hypothesis-generation tool. **Not a medical device.**

**Full usage guide:** [../USAGE.md](../USAGE.md)

## One-line install (from GitHub)

```bash
pip install "voc-breath[stack] @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

Optional extras:

```bash
pip install "voc-breath[dev,stack] @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

Or:

```bash
git clone https://github.com/hamcoderfran/tcga_mutations_survival.git
cd tcga_mutations_survival/exhalepath_atlas
pip install -e ".[dev,stack]"
# or
bash install_voc.sh
```

You get the `voc` command (aliases: `exhalepath`, `voc-breath`).

## Easy usage

```bash
# Naturalistic patient → template → 20-model stack
voc patient "35M with schizophrenia, smokes, on olanzapine, BMI 32, hallucinations"

# Inspect template only
voc patient --show-template "24yo obese male with depression"

# Structured stack
voc stack "depression" -l brain -c obesity --age 24 --sex male

# Holdout + break eval (public breath / literature / priority-10)
voc eval-stack-holdout
```

Artifacts land in `runs/stack_*/` and `runs/stack_holdout/`.

## What you get

- 20 fused models (VOC + genetics + flux + ADME + microbiome + …)
- PatientTemplate for vignettes, clinic notes, JSON, key=value
- HTML/PNG/CSV/JSON comprehensive reports

## Honesty bar

Cutting-edge as an **open systems / mechanism stack**.  
**Not** clinical SOTA vs prospective GC-MS/PTR diagnostic models. See `voc eval-stack-holdout` report.
