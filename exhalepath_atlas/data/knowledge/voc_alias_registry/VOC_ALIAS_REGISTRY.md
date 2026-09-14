# VOC alias / m/z registry

Generated: 2026-09-14T08:38:52.878522+00:00

> Unified name/m/z → atlas voc_id map for partner ingest diligence. Conflicts and uncertain Magdeburg channels are retained; do not treat as chemical identity certification.

- aliases: **95**
- Magdeburg m/z rows: **8**

## Confounder-sensitive seeds

- `smoking`: `acetonitrile`, `toluene`, `benzene`, `xylene`, `styrene`, `furan`, `2_butanone`
- `oral_microbiome`: `hydrogen_sulfide`, `methyl_mercaptan`, `indole`, `dimethyl_disulfide`, `dms`
- `diet_gut`: `acetone`, `isoprene`, `ethanol`, `methanol`, `trimethylamine`, `butyric_acid`, `acetic_acid`
- `ambient_exogenous`: `toluene`, `benzene`, `ethylbenzene`, `limonene`, `hexane`

## Regenerate

```bash
voc audit-voc-aliases
```
