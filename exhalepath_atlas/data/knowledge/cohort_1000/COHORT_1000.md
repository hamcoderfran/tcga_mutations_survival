# 1000-patient diversity cohort evaluation

**OK rate: 100.0%** (1000/1000)

## Resources engaged

- diseases=110 VOCs=50 pathways=17 tissues=73 cell_states=18 chains=19
- datasources: bindingdb, blood_proxy, gdc, gtex, hbdb, hmdb, kegg, metabolomics, mvoc, nist, opentargets, partition, pubchem, reactome
- chembl_loaded=True · OT offline=True · GDC=True

## Breakdown flags

- none

## Focus triad VOC alignment (centroid L2)

- `chronic_bronchitis__vs__copd`: 2.471
- `chronic_bronchitis__vs__lung_adenocarcinoma`: 4.395
- `chronic_bronchitis__vs__lung_squamous_cell_carcinoma`: 2.463
- `copd__vs__lung_adenocarcinoma`: 4.462
- `copd__vs__lung_squamous_cell_carcinoma`: 2.450
- `lung_adenocarcinoma__vs__lung_squamous_cell_carcinoma`: 2.623

## Smoking effect (|mean Δlog2fc| current vs never)

- **copd**: n_never=24 n_current=68 mean_abs_delta=0.06307170718430313
- **lung_adenocarcinoma**: n_never=32 n_current=47 mean_abs_delta=0.039641774333837775
- **chronic_bronchitis**: n_never=27 n_current=39 mean_abs_delta=0.04449676569597619

## Embeddings

- **voc**: n=1000 PC1=0.20397223050873947 PC2=0.13848151457392027 silhouette=0.024194189341762696 umap_error=None
- **gene_shift**: n=871 PC1=0.25589652594193946 PC2=0.15721077942899792 silhouette=0.14465866607079625 umap_error=None

## Outputs

- `patients.csv` — covariates + summary metrics
- `patients_dense.csv` — covariates + all VOC log2fc + pathway scores
- `patient_voc_matrix.csv` / `patient_voc_long.csv`
- `patient_gene_shift_matrix.csv` — genetic Δlog2fc vectors
- `patient_pathway_matrix.csv`
- `breakdown_flags.csv`
- `embeddings_voc.csv` / `embeddings_gene_shift.csv`
- `figures/voc_pca.png`, `voc_umap.png`, `gene_shift_*.png`, `*_focus_pulmonary.png`
