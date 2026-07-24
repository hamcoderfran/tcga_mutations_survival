# Real breath VOC training (no synthetic labels)

## What changed

- **Removed** offline-demo synthetic VOC targets (mechanistic prior + Gaussian noise → fake ppb).
- **Training y** now comes only from:
  1. Metabolomics Workbench measured cohorts (`ST000587` heart failure, `ST000883` malaria, `ST003200` healthy baselines)
  2. Scientific Data clinical breathomics peak tables (asthma / COPD / bronchiectasis, Figshare `23522490`)
  3. Curated literature panels with DOIs (`data/real_breath/literature_panels/`)
- **Priority diseases** supervised: asthma, COPD, COVID-19, bacterial pneumonia (CAP/HAP/VAP), TB, CF, OSA, gastric / H&N / prostate cancer, heart failure, malaria, ARDS
- Pathway / gene / cell-state atlas still drives **zero-shot inference** for unseen diseases — but those priors are never written as training labels

## Commands

```bash
voc build-real-corpus          # measured + literature labels only
voc train                      # calibrator real-1.0
voc eval-priority10            # directional check vs literature panels
voc "asthma" -l lung --top 20
voc "COVID-19" -l lung
voc "malaria"
voc "ARDS" -l lung
```

`voc build-corpus --offline-demo` is **blocked** unless `VOC_ALLOW_SYNTHETIC=1`.

## Other bulk breath databases (access audit)

```bash
voc harvest-alt-breath-sources   # ACCESS_MATRIX + Zenodo adjuncts + PhysioNet catalog
```

| Source | Bulk status here |
|--------|------------------|
| Clinical Breathomics (Figshare) | **Open** — already in `data/public_breath/` |
| HBDB live HTML | Cloudflare 403; use Zenodo SQL (separate scrape PR) |
| Owlstone VOC Atlas | Registration + **AI train/validate license ban** — not ingested |
| HMDB downloads | Cloudflare 403; PubChem-bridged seed annotations kept |
| Breathomix BreathBase / Shirley Atlas | Gated / no public dump |
| PhysioNet | Open catalog (waveforms ≠ VOC chemistry) |
| Zenodo breathomics adjuncts | **Open** — common-78 features, Tedlar variability, HBDB eval index |

## Honesty note

`eval-priority10` directional accuracy can look very high because literature panel directions also inform `voc_log2fc_prior`. Treat it as consistency with published direction, not a fully held-out GC-MS challenge. External cohort holdouts are the next step for true generalization scores.
