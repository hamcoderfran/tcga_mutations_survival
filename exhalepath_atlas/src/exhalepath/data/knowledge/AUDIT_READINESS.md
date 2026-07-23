# Full readiness audit (post-harden)

Generated after fixing P0–P2 findings from branch audit + Bugbot.

## Gate status

| Gate | Result |
|---|---|
| `pytest -q` | **105 passed** |
| `voc eval-stress-hard` | **100%** (50/50) |
| `voc eval-coverage --offline` | **100%** open-compound (777/777), SHA-bound |
| `voc eval-disease100-external` | Held-out lit-bench **100%**; disease100 **100** diseases |
| `voc eval-vision` | Vision fidelity **98.05%** |

## Bugs fixed in this pass

1. **Double-fuse ratchet** — fuse always restarts from `disease_voc_priors.pristine.json` (idempotent).
2. **dry_run wrote catalogs** — identity enrichment gated; dry-run is read-only.
3. **External validation dropped T2D/TB** — prefer diseases with external GT first.
4. **doi.org → S3 SSRF** — storage CDN redirects only from figshare/zenodo/clowder origins; HTTPS-only; hop-by-hop allowlist.
5. **Histology guard silent** — logged in changelog with `histology_guard:*` provenance.
6. **Integrity absolute paths** — relative paths + `paths_are_relative` policy.
7. **MW expanded catalog wipe** — merge curated + previous expanded + new discoveries.
8. **MW non-breath pollution** — require breath/exhal/EBC (or human+VOC); exclude plant/food.
9. **False 100% chem physchem** — split identity vs physchem metrics.
10. **cirrhosis pair miss** — `chronic_liver_disease` ↔ HCC.
11. **Coverage MD CWD path** — uses `KNOWLEDGE_DIR`.
12. **Circular metrics honesty** — held-out is primary; in-sample marked `circular_with_fuse`.
13. **Prior backups** — `disease_voc_priors.bak-<ts>.json` (gitignored).
14. **request_json body cap** — 32 MiB max.

## Remaining (blocked on user data, not code)

- HBDB Zenodo SQL (~1140 compounds)
- Owlstone Breath Biopsy VOC Atlas (gated)
- PubChem physchem live fill (CIDs present; mw/xlogp often null offline)

Calibrator `*.joblib` untouched.
