# Locked patient-level truth leaderboard

Generated: 2026-09-14T08:38:52.873281+00:00

> Patient-level nested AUROC leaderboard for research diligence only. Not clinical diagnostic SOTA. Prefer nested_auroc + optimism_gap + split_sha256 over directional literature percentages.

| Study | Disease | Signature | Method | Nested AUROC | Optimism gap | Split SHA256 | n |
|---|---|---|---|---:|---:|---|---:|
| ST000883 | — | literature | dot | 65.0% | 4.6% | `—` | — |
| ST000883 | — | literature | cosine | 63.3% | 2.7% | `—` | — |
| ST000587 | heart_failure | hybrid | cosine | 60.0% | -4.7% | `c56043631bf5…` | 23 |
| ST000883 | — | stack | cosine | 58.3% | -0.8% | `—` | — |
| ST000883 | — | hybrid | dot | 58.3% | 0.2% | `—` | — |
| ST000883 | — | stack | dot | 55.0% | -0.1% | `—` | — |
| ST000883 | malaria | hybrid | cosine | 53.3% | 3.5% | `e1b4271ace3b…` | 35 |

## How to regenerate

```bash
voc eval-leaderboard
voc eval-patient-diagnostic --all --signature hybrid
```
