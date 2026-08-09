# Is this SOTA / cutting edge?

**Short answer:** Cutting-edge as an *open multi-model systems breath-hypothesis + GC-MS research-enablement* stack. **Not** clinical diagnostic SOTA (no locked multi-site GC-MS/PTR patient ROC claim).

## What “SOTA” would mean here

| Claim | Our status |
|---|---|
| Best published ROC-AUC on a locked multi-site GC-MS/PTR patient cohort | **No** — we do not claim this |
| Best open mechanism-aware disease→exhaled-VOC reasoning stack | **Strong / novel** — 20-head adaptive fusion |
| Patient-level GC-MS diagnostic *research harness* (AUROC, locked splits, paper pack) | **Yes — shipped** (`voc eval-patient-diagnostic`, `voc export-paper-pack`) |
| Sci Data per-sample adapters + stratified AUCs + MetaboLights scaffold | **Yes — shipped** (`voc eval-scidata-samples`, `voc export-metabolights`) |
| Beats hybrid ExhalePath alone on public directional panels | **Yes** on current public holdouts (see below) |
| Useful zero-shot VOC directions for unseen/rare diseases | **Yes** — pathway→VOC projection + literature theme evidence + phenotype/MONDO |

## Holdout evidence (`voc eval-stack-holdout`)

Latest regenerated report: `data/knowledge/stack_holdout/STACK_HOLDOUT_REPORT.md`

| Benchmark | Hybrid dir | Stack dir |
|---|---:|---:|
| Public breath | 98.2% | **100.0%** |
| Literature | 100.0% | **100.0%** |
| Priority-10 | 100.0% | **100.0%** |
| PatientTemplate hard | — | **13/13** |
| PatientTemplate soft | — | **5/5** |

Anti-dilution fusion closed the earlier priority-10 regression (malaria hexanal was flipped by proxy heads when physiology disagreed with hybrid/calibrator).

> Directional panels can still partially overlap atlas priors — treat as systems-level evidence, not clinical validation.

## GC-MS / diagnostic research enablement

The scarce public good is not another optimistic directional %. It is **patient-level**, leakage-aware evaluation that other labs can cite and extend.

```bash
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-patient-diagnostic --all --signature stack
voc eval-scidata-samples --all
voc export-metabolights --study scidata:asthma
voc export-paper-pack runs/scidata_samples/asthma
voc lock-split --study ST000883
```

Reports land in `runs/patient_diagnostic/` and under `data/knowledge/gcms_diagnostic/`. See [RESEARCH.md](../RESEARCH.md) and [GCMS_DIAGNOSTIC.md](GCMS_DIAGNOSTIC.md).

## Cutting-edge stack (v2 / voc-breath 1.6)

1. **Adaptive mode-aware fusion** — atlas vs zero-shot family multipliers  
2. **Anti-dilution anchor** — hybrid+calibrator (preferred) or full 3-anchor agreement  
3. **Epistemic UQ** — per-VOC weighted std + approximate 90% CI  
4. **Calibrated fusion weights** — `fusion_weights_calibrated.json`  
5. **`zero_shot_evidence`** — `expected_voc_direction` + mechanism themes  
6. Hardened **`zero_shot_mechanism`** — pathway_bias / gene-seeded `voc_effects` projection  
7. **`phenotype_mondo`**, **`counterfactual_null`**, **`meta_ensemble`**  
8. Expanded zero-shot prior pack (34 curated rare/novel disease entries)  
9. Patient-level GC-MS matrices + AUROC/AUPRC/sens/spec + locked SHA256 splits  
10. Paper-ready TRIPOD+AI / BreathVOC checklist stubs  

## Why not clinical SOTA (yet)

1. No prospective multi-site GC-MS/PTR trial with locked external labels  
2. Bundled patient GC-MS n is small (ST000883 ≈35; ST000587 ≈23)  
3. Sci Data panels are cross-cohort differentials, not absolute healthy-controlled ppb  
4. Literature/priority panels can overlap atlas priors → optimistic directional scores  
5. Human-GEM / OPERA / PrimeKG integrations remain **proxy packs**  
6. Confounders (smoking/age) often missing from public MW factors  

## Why it is cutting-edge (open research)

- 20-model fusion across VOC quantity + systems biology  
- Mode-aware anti-dilution fusion with epistemic uncertainty  
- Naturalistic PatientTemplate → stack  
- Dedicated zero-shot evidence channel  
- Patient-level GC-MS diagnostic research harness  
- One-pip installable package with offline priors + calibrated weights  

## Packaging for others

```bash
pip install "voc-breath[stack] @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
voc patient "…"
voc stack --nl "Maple syrup urine disease, genes BCKDHA BCKDHB"
voc eval-patient-diagnostic --study ST000883
voc eval-stack-holdout
```

See [USAGE.md](../USAGE.md), [INSTALL.md](INSTALL.md), and [RESEARCH.md](../RESEARCH.md).
