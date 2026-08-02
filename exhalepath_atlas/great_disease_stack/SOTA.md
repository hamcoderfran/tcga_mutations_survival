# Is this SOTA / cutting edge?

**Short answer:** Cutting-edge as an *open multi-model systems breath-hypothesis stack*. **Not** state-of-the-art as a clinically validated breath diagnostic.

## What “SOTA” would mean here

| Claim | Our status |
|---|---|
| Best published ROC-AUC on a locked multi-site GC-MS/PTR patient cohort | **No** — we do not claim this |
| Best open mechanism-aware disease→exhaled-VOC reasoning stack | **Strong / novel** among open tools |
| Beats hybrid ExhalePath alone on public directional panels | **Mixed** — see holdout (public breath ≥ hybrid; priority-10 slightly below) |

## Holdout evidence (regenerate with `voc eval-stack-holdout`)

See `data/knowledge/stack_holdout/STACK_HOLDOUT_REPORT.md` after running the eval.

Typical findings on current public panels:

- Public breath directional: stack ≈ hybrid (high; **partly circular** with priors)
- Literature directional: high for both
- Priority-10: hybrid can edge stack when fusion dilutes calibrated signals
- PatientTemplate adversarial hard pass: high-90%s after hardening

## Why not clinical SOTA (yet)

1. No prospective patient-level GC-MS/PTR holdout with locked labels  
2. Sci Data panels are cross-cohort differentials, not absolute healthy-controlled ppb  
3. Literature/priority panels overlap atlas priors → optimistic directional scores  
4. Human-GEM / OPERA / PrimeKG integrations are **proxy packs**, not full external binaries  
5. No head-to-head vs published sensor-array / cohort ML baselines on identical splits  

## Why it is still cutting-edge (open research)

- 16-model fusion across VOC quantity, genetics, flux, ADME/PBPK, microbiome, signaling, cell state, comorbidity, pharmacology  
- Naturalistic PatientTemplate → stack for diverse clinical text  
- Agreement / uncertainty and per-model vote audit trails  
- One-pip installable package with offline priors  

## Packaging for others

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
voc patient "…"
voc eval-stack-holdout
```

See [INSTALL.md](INSTALL.md).
