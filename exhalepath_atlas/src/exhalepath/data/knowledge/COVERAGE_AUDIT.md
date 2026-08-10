# VOC data coverage audit (secure open corpus)

This audit does NOT claim 99% of all VOC content on the internet. It targets 99% of the defined open machine-readable breath-VOC compound universe (EPA VOLATILOME / HBDB backbone). HBDB live HTML, restricted SQL, Owlstone gated Atlas, and exhaustive PDF scraping remain out of scope.

**Open-compound coverage: 100.0%** (777/777) · meets ≥99% target: **True**

## What “99%” means here

Denominator = EPA VOLATILOME open breath-chemical list (HBDB literature backbone), secured locally with SHA-256 provenance. Not included: paywalled PDFs, Cloudflare-blocked HBDB HTML, gated Owlstone Atlas.

## Metrics

- Atlas prediction VOCs: 50
- Extended catalog mapped to atlas IDs: 44
- Atlas chem enrichment (PubChem/HMDB): 100.0%
- Curated public studies: 9
- Expanded MW study catalog size: 25 (+16 discovered)
- Literature panels: 16 · panel DOIs: 42
- Europe PMC breath-VOC metadata records: 100 (DOIs: 97) — metadata only, no PDF scrape

## Security / anti-poisoning

- host allowlist (SSRF)
- max download size
- executable magic-byte reject
- shebang reject
- HTML refuse by default
- SHA-256 provenance sidecars
- JSON depth/required-key validation
- quarantine-then-promote writes

## Remaining gaps

- **hbdb_full_1140**: HBDB claims ~1140 compounds / 2766 refs; live DB blocked — using VOLATILOME 777 as open substitute
- **owlstone_atlas**: Breath Biopsy VOC Atlas is gated; not bulk-downloadable
- **prediction_panel_50**: Mechanistic predictor still uses 50 VOCs; extended catalog is inventory/coverage, not full retrain
- **pdf_fulltext**: Europe PMC harvest is metadata-only; no automated full-text VOC table extraction yet

## Artifacts

- `data/datasources/coverage_universe.json`
- `data/knowledge/voc_extended_catalog.json`
- `data/datasources/INTEGRITY_MANIFEST.json`
- `data/datasources/literature/europepmc_breath_voc_metadata.json`
- `data/datasources/metabolomics/breath_study_catalog_expanded.json`

Generated: 2026-08-10T02:13:21.385096+00:00
