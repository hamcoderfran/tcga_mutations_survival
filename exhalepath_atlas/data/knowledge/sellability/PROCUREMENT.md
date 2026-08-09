# Procurement pack — ExhalePath / voc-breath

**Lead line:** We cut the cost of being wrong about VOCs.

**Not selling:** clinical diagnostic SOTA, FDA/CE clearance, Sci Data 100% OVR as accuracy, or n≈35 AUROC as product proof.

## Product shape for buyers

| Offer | What they get | Typical ACV |
|---|---|---|
| Pharma biomarker seat | `demo-close`, ledger-cited hypotheses, paper zip, locked splits | $150–600k/yr |
| Platform OEM embed | `exhalepath.oem` library + `voc serve-api` inside their LIMS | $0.5–3M/yr + milestones |
| CRO white-label | Partner OMNI ingest + LOSO diligence + MetaboLights scaffolds | $100–500k/yr |

## One-command diligence

```bash
# Closing demo (note → VOCs + optimism gap + paper zip)
voc demo-close --disease malaria

# Partner multi-site LOSO (bundled OMNI-style fixture or your CSVs)
voc diligence-loso
voc diligence-loso --site-a siteA.csv --site-b siteB.csv

# One-click import of their feature table
voc import-breathvoc partner_table.csv --format omni --disease-id malaria

# OEM API (customer VPC)
voc serve-api --host 127.0.0.1 --port 8787
curl -s localhost:8787/v1/security | jq .
```

## Security / SOC2-ish story

See `SECURITY.md` in this folder and `GET /v1/security` on the OEM API.

- Egress allowlist via `secure_fetch` (SSRF guard)
- No PHI store in default mode (research vignette fields only)
- CSV/JSON feature-table import only
- Customer VPC / air-gapped wheel deploy for regulated buyers
- **Not** a completed SOC 2 Type II attestation — architecture checklist for auditors

## BreathVOC interchange

Partners ship `BreathVOC-1.1` JSON (schema next to exports) or OMNI-style sample×compound CSV.  
`voc import-breathvoc` normalizes both into the same patient×VOC matrix used by LOSO and paper packs.
