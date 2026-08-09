# Security & compliance story (SOC2-ish)

**Scope:** R&D software for VOC hypothesis + diligence. Not a medical device. Not a PHI system of record.

## Controls we ship

| Control | Implementation |
|---|---|
| Egress allowlist | `exhalepath.ingest.secure_fetch.ALLOWED_HOSTS` — scientific hosts only |
| SSRF / redirect guard | Rejects redirects off allowlist (storage CDNs only when explicitly permitted) |
| Upload hygiene | Feature tables are CSV/JSON; secure_fetch rejects executable extensions |
| No central PHI | Default `PatientTemplate` is vignette/research fields; deploy air-gapped |
| Deterministic audit | Locked-split `content_sha256`; paper-pack zip sha256 |
| Supply chain | Pip wheel + optional host allowlist for open coverage audits |

## SOC 2 mapping (customer-deployed)

| Trust criteria | Customer responsibility | Our contribution |
|---|---|---|
| Security | VPC, IdP, mTLS in front of `voc serve-api` | Allowlisted egress, no arbitrary code exec on import |
| Availability | Their HA / backups of `runs/` | Stateless API + local wheel |
| Processing integrity | Validate imports; pin versions | BreathVOC schema validation; split verify |
| Confidentiality | Volume encryption; no PHI into vignettes | No cloud PHI store in default mode |
| Privacy | DPA / BAA if they add PHI | Out of default scope — research vignettes |

## Explicit gaps

- Not a completed SOC 2 Type II report  
- Not HIPAA BAA by default  
- Not FDA 21 CFR Part 11 validated  

Use this checklist with your auditor; do not claim “SOC2 certified” from this document alone.
