# tcga_mutations_survival

## ExhalePath Atlas (primary deliverable)

**`exhalepath_atlas/`** — fully integrated exhaled VOC biomarker platform:

- Priority datasources **1–12** (metabolomics repos → HMDB → λ → mVOC → PubChem → Reactome → GTEx → Open Targets → GDC → BindingDB → blood proxy → NIST RI)
- Whole-body model for **100+ diseases**
- Top-50 VOC Δppb + mechanism WHY (pathways, Census cells, genes)

```bash
cd exhalepath_atlas
pip install -e ".[dev]"
python -m exhalepath integrate-datasources
python -m exhalepath biomarker "lung adenocarcinoma" --location lung --genes KRAS,TP53 --top 50
```

See [`exhalepath_atlas/README.md`](exhalepath_atlas/README.md).

## Other folders

- `exhalepath/` — earlier iterative package (still maintained on this branch)
- `project1_tcga_mutations_survival/` — original TCGA mutation–survival analysis
