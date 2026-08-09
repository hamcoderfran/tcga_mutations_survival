# Industry benchmark pack

Generated: 2026-08-09T19:06:41.915107+00:00

> Cutting-edge open research-enablement stack for exhaled VOC / GC-MS R&D. Not a clinical diagnostic SOTA claim; not FDA/CE clearance.

Smoke: **8/8** steps OK (PASS)

## Scorecard

| Metric | Value | Standard | Why it matters |
|---|---|---|---|
| Nested patient AUROC (ST000883 hybrid) | 53.3% | TRIPOD+AI / JBR nested CV | Patient-level holdout — unit of clinical inference; beats peak-level CV optimism |
| Optimism gap (non-nested − nested) | 3.5% | JBR 2024 breath-ML leakage audit | Quantifies overfit risk buyers must price into VOC programs |
| Same-feature nested logistic ceiling | 73.3% | Internal ML baseline | Upper bound when labels are fit on this tiny cohort — mechanism signatures are transferable, not maxed here |
| AUPRC / Youden sens-spec | AUPRC=53.1%; sens=94.1%; spec=33.3% | STARD / diagnostic operating point | Discrimination + thresholded clinical operating characteristics |
| Locked split SHA256 | e1b4271ace3b6d550d764e75… | Preregistration / audit trail | Citeable split integrity for partner diligence |
| Sci Data per-sample OVR AUROC | asthma=100.0%; copd=100.0%; bronchiectasis=100.0% | Multi-cohort peak-table ML | Per-sample (not cohort-mean) discrimination across pulmonary cohorts; NO healthy arm |
| Stack vs hybrid directional (public breath) | stack=100.0% hybrid=98.2% | Systems holdout | Multi-head fusion vs hybrid baseline on public panels |
| Stack directional (literature / priority-10) | lit=100.0%; p10=100.0% | Systems holdout | Directional panels — treat as systems evidence; partial prior overlap possible |
| PatientTemplate adversarial break | hard=13/13 soft=5/5 | Adversarial / phenotype stress | Stack must not collapse under comorbidity/smoking template attacks |
| Public breath directional accuracy | 98.2% | Literature/public panel concordance | Hypothesis-generation accuracy on curated elevated/suppressed panels |
| Open-compound coverage vs VOLATILOME universe | 100.0% (777/777) | Open data diligence | Honest ≥99% target on *open* compound universe — not all gated industry atlases |
| Literature concordance (disease suite) | 100.0% | Panel overlay | Directional agreement vs published breath panels (circularity risk flagged) |
| Paper pack completeness | 8 files · zip sha 2b678ccbff26… | Publication / partner diligence pack | One zip: ROC + Methods + overlay + split hash — hours→minutes for BD packs |
| Clinical diagnostic SOTA claim | NONE — research enablement only | Regulatory honesty | Required for credible enterprise sales; overclaim destroys diligence |

## Step timings

| Step | OK | Seconds |
|---|---|---:|
| patient_diagnostic_ST000883 | True | 0.958 |
| signature_benchmark_ST000883 | True | 4.643 |
| scidata_per_sample | True | 8.315 |
| stack_holdout | True | 15.767 |
| public_breath | True | 0.807 |
| coverage_audit_offline | True | 0.033 |
| lit_compare | True | 11.012 |
| deposit_scaffolds | True | 1.663 |

## Draft branch prune

See `DRAFT_BRANCH_PRUNE.md` and `BUYER_BRIEF.md`.
