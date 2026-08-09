# Draft branch prune memo

## Cherry-picked (additive)

- secure_fetch.py + http allowlist (#14)
- coverage_audit.py + COVERAGE_AUDIT artifacts (#14)
- lit_compare.py + lit_compare artifacts (#13) — NOT priors/predict wholesale
- post-blend smoking/age exo fix in predict.py (surgical from #13)
- ds15_alt_breath_sources + catalogs (#17)
- HBDB 60-disease JSON + extract_hbdb_zenodo.py (#16)
- HMDB Wishart mirror constant + breath metabolite extracts (#18)

## Skipped (regressive / conflicting)

- Wholesale cli.py / predict.py from draft tips (would strip zero-shot hardening)
- disease_voc_priors.json replacements (#13/#15/#16)
- Recalibrated voc_calibrator.joblib from draft tips
- model-improve-100disease full merge (#15)
- Large optional Tedlar xlsx / hbdb_eval zip raw dumps (catalogs kept)

**Action:** Close or leave draft PRs #13–#19; do not merge wholesale into main
