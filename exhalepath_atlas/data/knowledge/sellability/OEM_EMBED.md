# OEM embed kit — LIMS / platform integration

Recurring ACV beats one-off licenses. Embed these three calls.

## Python (in-process)

```python
from exhalepath.oem import score_sample, verify_locked_split, import_feature_table

# 1) Score a measured VOC vector (ledger-cited)
out = score_sample(
    "malaria",
    {"pentane": 0.4, "hexanal": 0.5, "acetone": -0.2, "isoprene": -0.3},
    signature_source="hybrid",
)
# out["top_vocs"][*].evidence_grade / .doi

# 2) Verify a preregistered locked split
chk = verify_locked_split("data/knowledge/gcms_splits/ST000883_split_v1.json")
assert chk["ok"]

# 3) Import partner OMNI / BreathVOC table
imp = import_feature_table("partner_omni.csv", fmt="omni", disease_id="malaria")
matrix = imp["matrix"]
```

## HTTP API

```bash
voc serve-api --port 8787

curl -s localhost:8787/health
curl -s localhost:8787/v1/security
curl -s -X POST localhost:8787/v1/score-sample \
  -H 'content-type: application/json' \
  -d '{"disease_id":"malaria","vocs":{"pentane":0.4,"hexanal":0.5,"acetone":-0.2}}'
```

## Integration checklist

1. Pin `voc-breath` version; verify wheel hash  
2. Put TLS + IdP in front of `serve-api`  
3. Persist imported matrix content hashes in your LIMS audit log  
4. Surface **optimism gap** + **evidence_grade** in the UI — never raw AUROC alone  
5. Keep clinical claims on *your* prospective data; ExhalePath is the diligence layer  

## Positioning in the product UI

> Cut the cost of wrong VOC panels — locked nested eval and ledger-cited hypotheses before assay spend.
