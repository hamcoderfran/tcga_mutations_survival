# Buyer brief — selling VOC R&D enablement (honest)

**One-liner:** We cut the cost of being wrong about VOCs: locked patient-level eval, mechanism hypotheses for unseen diseases, and deposit-ready paper packs — so your VOC program spends money on real signal, not optimistic CV.

**Never lead with:** clinical AUROC, Sci Data 100% OVR, or n≈35 malaria numbers as product proof.

**Positioning:** ExhalePath / voc-breath is an open **VOC R&D operating system**: mechanism-aware disease→VOC hypothesis generation + patient-level GC-MS eval harness + deposit/reporting scaffolds. Buyers are breath-biopsy platforms, pharma biomarker teams, and metabolomics CROs entering VOCs — not hospital IVD.

## Sellability pack (how to close)

```bash
voc demo-close --disease malaria          # note → ledger-cited VOCs + optimism gap + paper zip
voc diligence-loso                        # partner OMNI-style multi-site LOSO diligence slide
voc import-breathvoc their_table.csv      # one-click BreathVOC / OMNI import
voc serve-api --port 8787                 # OEM / procurement JSON API
```

Docs: `data/knowledge/sellability/` (`PROCUREMENT.md`, `OEM_EMBED.md`, `SECURITY.md`).

## What is actually revolutionary (evidence-backed)

- Patient-level nested AUROC + locked split SHA256 on public GC-MS (rare in open stacks)
- Optimism-gap reporting that prevents $1–10M wasted follow-ups on peak-level CV mirages
- Sci Data per-sample adapters (not just cohort means) with age/sex strata + paper zip
- 20-head mechanism stack with anti-dilution fusion + epistemic UQ for zero-shot diseases
- MetaboLights mzTab-M/ISA-Tab scaffolds + TRIPOD stub — BD diligence in one command
- Secure open-coverage audit with host allowlist (supply-chain hygiene for enterprise IT)

## What is NOT revolutionary (do not overclaim)

- Not multi-site prospective clinical AUROC SOTA
- Not FDA/CE cleared; not a medical device
- Sci Data OVR AUROCs can look 'perfect' because cohorts are chemically distinct pulmonary diseases without healthy controls — do not sell as clinical accuracy
- Literature directional % can partly overlap atlas priors (circularity)

## Pain → savings levers

| Pain | Cost of status quo | Our lever |
|---|---|---|
| Peak-level CV / leakage → false greenlight of VOC panels | 0.5M–8M per failed biomarker program (assay, cohort, BD time) | Nested patient AUROC + optimism gap + locked splits |
| Months assembling Methods / TRIPOD / MetaboLights deposit packs | 50k–400k FTE + delayed partnership diligence | voc export-paper-pack + export-metabolights in minutes |
| Zero-shot rare disease VOC hypotheses require PhD weeks | 80k–250k per indication scout | stack / zero-shot mechanism heads + literature overlay |
| Gated Owlstone-class atlases block early scouting | license + delay; opportunity cost | open HBDB/Sci Data/MW/HMDB alt-source matrix (ds15) for pre-license triage |

## Credible deal math (not fake $B from n=35 AUROC)

Billions of enterprise value require either (a) platform adoption across many pharma programs or (b) clinical IVD clearance with reimbursable use. Today's evidence supports (a)-style **R&D software / partnership** economics, not (b). Do not pitch a $B valuation on ST000883 n≈35 AUROC.

### Near-term ACV

- **single_pharma_biomarker_seat**: 150k–600k / year
- **breath_platform_oem_embed**: 0.5M–3M / year + milestone
- **CRO white-label eval suite**: 100k–500k / year

### Path to larger outcomes

- Become default pre-clinical VOC hypothesis + eval layer for 3–5 breath platforms
- Accumulate multi-site locked-split benchmarks partners cannot ignore
- Optional later: IVD partnership where *their* prospective data carry clinical claims

TAM context: Breath biopsy / exhaled biomarker tooling sits inside a multi-billion non-invasive diagnostics + pharma biomarker adjacent market; software attach rates are a small slice — sell acceleration and de-risking, not 'we are the diagnostic'.

## Live scorecard snapshot

| Metric | Value |
|---|---|
| Nested patient AUROC (ST000883 hybrid) | 53.3% |
| Optimism gap (non-nested − nested) | 3.5% |
| Same-feature nested logistic ceiling | 73.3% |
| AUPRC / Youden sens-spec | AUPRC=53.1%; sens=94.1%; spec=33.3% |
| Locked split SHA256 | e1b4271ace3b6d550d764e75… |
| Sci Data per-sample OVR AUROC | asthma=100.0%; copd=100.0%; bronchiectasis=100.0% |
| Stack vs hybrid directional (public breath) | stack=100.0% hybrid=98.2% |
| Stack directional (literature / priority-10) | lit=100.0%; p10=100.0% |
| PatientTemplate adversarial break | hard=13/13 soft=5/5 |
| Public breath directional accuracy | 98.2% |
| Open-compound coverage vs VOLATILOME universe | 100.0% (777/777) |
| Literature concordance (disease suite) | 100.0% |
| Paper pack completeness | 8 files · zip sha 2b678ccbff26… |
| Clinical diagnostic SOTA claim | NONE — research enablement only |

## Diligence commands

```bash
voc eval-industry-pack
voc eval-patient-diagnostic --study ST000883 --signature hybrid
voc eval-scidata-samples --all
voc eval-coverage --offline
voc export-metabolights --study ST000883
```
