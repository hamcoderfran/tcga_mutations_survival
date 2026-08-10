# Mental-health literature eval (de-circularized)

Generated: 2026-08-10T19:36:19.252783+00:00

> Raw directional accuracy can be curation-circular when panel VOCs are also in atlas priors. panel_masked_prior removes those prior entries before prediction — that score is the de-circularized metric. Neither score is a clinical AUROC; no MH patient intensity cohort is bundled.

- Mean **raw** directional accuracy: 1.0
- Mean **panel-masked prior** directional accuracy: 0.9249999999999999

## Panel diseases

- **schizophrenia**: raw=1.0 (10/10); masked=0.9 (9/10); removed_prior=10
- **major_depressive_disorder**: raw=1.0 (8/8); masked=0.875 (7/8); removed_prior=8
- **bipolar**: raw=1.0 (2/2); masked=1.0 (2/2); removed_prior=2

## Thin-evidence conditions

- **anxiety**: no measured panel · `atlas_source=stress_proximal_hypothesis_thin_breath_lit`
- **ptsd**: no measured panel · `atlas_source=trauma_stress_proximal_hypothesis_thin_breath_lit`
- **adhd**: no measured panel · `atlas_source=literature_thin_ef_adjacent_not_adhd_cohort`
- **autism_spectrum_disorder**: no measured panel · `atlas_source=gut_microbiome_voc_proxy_not_exhaled_panel`
