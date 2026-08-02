# Signature benchmark — ST000883

| Signature | Method | AUROC | Nested AUROC | Nested logistic | AUPRC |
|---|---|---:|---:|---:|---:|
| hybrid | cosine | 56.9% | 53.3% | 73.3% | 53.1% |
| hybrid | dot | 58.5% | 58.3% | 73.3% | 67.6% |
| stack | cosine | 57.5% | 58.3% | 73.3% | 51.8% |
| stack | dot | 54.9% | 55.0% | 73.3% | 60.9% |
| literature | cosine | 66.0% | 63.3% | 73.3% | 70.5% |
| literature | dot | 69.6% | 65.0% | 73.3% | 74.7% |

Data-fit nested logistic is an upper reference on the same mapped VOC features. Mechanism signatures are transferable templates; literature panels can score higher when they encode the same cohort's published directions (partly circular). Report nested AUROC + optimism gap.
