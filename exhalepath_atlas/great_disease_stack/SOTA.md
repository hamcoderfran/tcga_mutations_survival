# Is this SOTA / cutting edge?

**Short answer:** Cutting-edge as an *open multi-model systems breath-hypothesis + GC-MS research-enablement* stack. **Not** state-of-the-art as a clinically validated breath diagnostic.

## What “SOTA” would mean here

| Claim | Our status |
|---|---|
| Best published ROC-AUC on a locked multi-site GC-MS/PTR patient cohort | **No** — we do not claim this |
| Best open mechanism-aware disease→exhaled-VOC reasoning stack | **Strong / novel** among open tools |
| Patient-level GC-MS diagnostic *research harness* (AUROC, locked splits, paper pack) | **Yes — shipped** (`voc eval-patient-diagnostic`) |
| Beats hybrid ExhalePath alone on public directional panels | See `voc eval-stack-holdout` |

## GC-MS / diagnostic research enablement

The scarce public good is not another optimistic directional %. It is **patient-level**, leakage-aware evaluation that other labs can cite and extend.

```bash
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-patient-diagnostic --all --signature stack
voc lock-split --study ST000883
```

Reports land in `runs/patient_diagnostic/` and are copied under `data/knowledge/gcms_diagnostic/` after regeneration. See [RESEARCH.md](../RESEARCH.md) and [GCMS_DIAGNOSTIC.md](GCMS_DIAGNOSTIC.md).

## Holdout evidence (directional panels)

See `data/knowledge/stack_holdout/STACK_HOLDOUT_REPORT.md` after `voc eval-stack-holdout`.

Typical findings:

- Public breath / literature directional: high for hybrid and stack (**partly circular** with priors)
- PatientTemplate adversarial: hardened high pass rates
- Patient-level MW GC-MS (ST000883): AUROC + nested optimism gap now reported separately

## Why not clinical SOTA (yet)

1. No prospective multi-site GC-MS/PTR trial with locked external labels  
2. Bundled patient GC-MS n is small (ST000883 ≈35; ST000587 ≈23)  
3. Sci Data panels are cross-cohort differentials, not healthy-controlled absolute ppb  
4. Literature/priority panels can overlap atlas priors  
5. Human-GEM / OPERA / PrimeKG integrations remain **proxy packs**  
6. Confounders (smoking/age) often missing from public MW factors  

## Why it is cutting-edge (open research)

- Multi-model fusion across VOC quantity + systems biology  
- Naturalistic PatientTemplate → stack  
- **Patient-level GC-MS matrices + AUROC/AUPRC/sens/spec + locked SHA256 splits**  
- Paper-ready TRIPOD+AI / BreathVOC checklist stubs  
- One-pip installable package with offline priors  

## Packaging for others

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
voc patient "…"
voc eval-patient-diagnostic --study ST000883
voc eval-stack-holdout
```

See [INSTALL.md](INSTALL.md) and [RESEARCH.md](../RESEARCH.md).
