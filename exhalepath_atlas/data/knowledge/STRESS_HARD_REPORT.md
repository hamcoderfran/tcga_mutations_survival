# Adversarial stress suite — before / after hardening

50 difficult cases (`stress_hard_50.json`) designed to break ExhalePath:
empty / garbage inputs, gene-as-disease traps, alias collisions, bogus
locations, comorbidity weight extremes, Magdeburg schizophrenia suppressions,
zero-shot near-healthy, and metastatic / multi-signal pressure.

## Results

| Run | Pass rate | Passed | Failed |
|---|---:|---:|---:|
| Baseline (pre-harden) | **88.0%** | 44 | 6 |
| After harden | **100.0%** | 50 | 0 |

### Baseline failures (what broke)

| Case | Failure mode | Root cause |
|---|---|---|
| `symbol_location` (`???`) | location falsely matched | `_norm("???")` → `""`; empty string is contained in every alias |
| `comorbidity_weight_overflow` (3.5) | `ValidationError` | `DiseaseQuery.comorbidity_weight` capped at 1.5 with no clamp at API boundary |
| `comorbidity_weight_negative` (−1) | `ValidationError` | same |
| `gene_BRCA1_as_disease` | resolved → breast cancer | alias `BRCA` substring-matched inside `brca1` |
| `sz_suppress_acetone` | acetone elevated | schizophrenia prior still had acetone `+0.255` |
| `sz_heart_suppress_tma` | acetone / TMA elevated | weak SZ suppress + heart_disease comorbidity ketone/TMA push |

### Hardenings applied

1. **Location garbage** (`body/tissues.py`): empty normalized keys return unmatched custom tissue (no fuzzy).
2. **Gene-as-disease** (`knowledge/loader.py`): HGNC-like symbols (digit-bearing or pathway seed genes) resolve as unresolved after exact-alias check; fuzzy containment requires whole-token boundaries.
3. **Comorbidity weight** (`biomarker.py`, `nl/ask.py`): clamp to `[0, 1.5]` before schema validation.
4. **Schizophrenia Magdeburg panel** (`disease_voc_priors.json`): suppress acetone / isoprene / trimethylamine / methanol; elevate pentane / ethane.
5. **Autism Mondo hygiene**: `MONDO:0005260` (was colliding with heart failure `MONDO:0005258`).

### Category pass rates (after)

- robustness 7/7 · gene_trap 4/4 · alias 6/6 · location 5/5
- clinical_hard 9/9 · zero_shot 4/4 · conflict 9/9 · metastatic 6/6

### How to re-run

```bash
python scripts/build_stress_hard_50.py
voc eval-stress-hard --out-dir runs/stress_hard
pytest tests/test_stress_hardening.py -q
```
