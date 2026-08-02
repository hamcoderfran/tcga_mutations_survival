# Is this SOTA / cutting edge?

**Short answer:** Among open mechanism-aware disease→exhaled-VOC systems, this stack is built to be **second to none for comprehensive + zero-shot hypothesis generation**. It is still **not** clinical diagnostic SOTA (no locked multi-site GC-MS/PTR patient ROC claim).

## What “SOTA” would mean here

| Claim | Our status |
|---|---|
| Best published ROC-AUC on a locked multi-site GC-MS/PTR patient cohort | **No** — we do not claim this |
| Best open mechanism-aware disease→exhaled-VOC reasoning stack | **Target / strong** — 20-head adaptive fusion |
| Beats hybrid ExhalePath alone on public directional panels | **Goal of anti-dilution fusion** — re-check with `voc eval-stack-holdout` |
| Useful zero-shot VOC directions for unseen/rare diseases | **Yes** — pathway→VOC projection + literature theme evidence + phenotype/MONDO |

## Cutting-edge stack (v2)

Beyond the original 16 heads, the stack now includes:

1. **Adaptive mode-aware fusion** — atlas vs zero-shot family multipliers  
2. **Anti-dilution anchor** — when hybrid / physiology / calibrator agree on sign, fused log2fc is pulled toward that anchor (stops weak proxies from washing out calibrated signals)  
3. **Epistemic UQ** — per-VOC weighted std + approximate 90% CI on fused effects  
4. **Calibrated fusion weights** — `fusion_weights_calibrated.json` (holdout-aware anchor boost)  
5. **`zero_shot_evidence`** — `expected_voc_direction` + mechanism themes → real VOC votes  
6. **Hardened `zero_shot_mechanism`** — pathway_bias / gene-seeded `voc_effects` projection (no empty ZS heads)  
7. **`phenotype_mondo`** — free-text phenotype / MONDO → VOC themes  
8. **`counterfactual_null`** — deterministic null shrinkage for overclaim control  
9. **`meta_ensemble`** — second-order hybrid × literature sign agreement  
10. **Expanded zero-shot prior pack** — 34 curated rare/novel disease gene→VOC entries  

## Holdout evidence (regenerate with `voc eval-stack-holdout`)

See `data/knowledge/stack_holdout/STACK_HOLDOUT_REPORT.md` after running the eval.

What to look for after this upgrade:

- Public breath / literature directional: stack **≥** hybrid (anti-dilution)  
- Priority-10: stack should **not** regress below hybrid due to dilution  
- Zero-shot diseases (MSUD, PKU, ALS, CF, …): non-empty fused VOC panels with mechanism themes  
- PatientTemplate adversarial: hard/soft parse + stack still runnable  

## Why not clinical SOTA (yet)

1. No prospective patient-level GC-MS/PTR holdout with locked labels  
2. Sci Data panels are cross-cohort differentials, not absolute healthy-controlled ppb  
3. Literature/priority panels can overlap atlas priors → optimistic directional scores  
4. Human-GEM / OPERA / PrimeKG integrations remain **proxy packs**, not full external binaries  
5. No head-to-head vs published sensor-array / cohort ML baselines on identical splits  

## Why it is cutting-edge (open research)

- 20-model fusion across VOC quantity, genetics, flux, ADME/PBPK, microbiome, signaling, cell state, comorbidity, pharmacology, phenotype, null, and meta  
- Mode-aware anti-dilution fusion with epistemic uncertainty  
- Naturalistic PatientTemplate → stack for diverse clinical text  
- Dedicated zero-shot evidence channel (not just atlas copy)  
- One-pip installable package with offline priors + calibrated weights  

## Packaging for others

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
voc patient "…"
voc stack --nl "Maple syrup urine disease, genes BCKDHA BCKDHB"
voc eval-stack-holdout
# optional recalibration
python scripts/calibrate_stack_fusion.py
```

See [INSTALL.md](INSTALL.md) and [README.md](README.md).
