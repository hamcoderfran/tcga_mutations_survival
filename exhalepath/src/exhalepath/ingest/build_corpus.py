from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import DEFAULT_GDC_PROJECTS, PROCESSED_DIR
from ..knowledge.loader import KnowledgeBase
from ..utils import stage_ordinal
from .gdc_client import GDCClient
from .reactome_client import ReactomeClient


def _stage_to_num(stage: object) -> float | None:
    if not isinstance(stage, str) or not stage.strip():
        return None
    return stage_ordinal(stage, default=float("nan"))


def _build_offline_demo_tables(kb: KnowledgeBase, n_cases_per_project: int = 400) -> dict[str, pd.DataFrame]:
    """Synthetic multi-project corpus for offline training / CI (still large row counts)."""
    rng = np.random.default_rng(0)
    projects = [
        "TCGA-LUAD",
        "TCGA-BRCA",
        "TCGA-PAAD",
        "TCGA-COAD",
        "TCGA-LIHC",
        "TCGA-OV",
    ]
    sites = {
        "TCGA-LUAD": "Lung",
        "TCGA-BRCA": "Breast",
        "TCGA-PAAD": "Pancreas",
        "TCGA-COAD": "Colon",
        "TCGA-LIHC": "Liver",
        "TCGA-OV": "Ovary",
    }
    gene_sets = kb.pathway_gene_universe()
    clin_rows = []
    edge_rows = []
    top_rows = []
    for pid in projects:
        for i in range(n_cases_per_project):
            case_id = f"{pid.lower()}-case-{i:05d}"
            stage = rng.choice(["Stage I", "Stage II", "Stage III", "Stage IV"], p=[0.25, 0.3, 0.3, 0.15])
            clin_rows.append(
                {
                    "project_id": pid,
                    "case_id": case_id,
                    "submitter_id": case_id,
                    "primary_site": sites[pid],
                    "disease_type": "Adenomas and Adenocarcinomas",
                    "gender": rng.choice(["female", "male"]),
                    "vital_status": rng.choice(["Alive", "Dead"], p=[0.65, 0.35]),
                    "days_to_death": rng.integers(30, 3000),
                    "age_at_diagnosis_days": int(rng.normal(60, 12) * 365.25),
                    "ajcc_pathologic_stage": stage,
                    "primary_diagnosis": "Adenocarcinoma, NOS",
                    "tissue_or_organ_of_origin": sites[pid],
                    "site_of_resection_or_biopsy": sites[pid],
                    "morphology": "8140/3",
                    "classification_of_tumor": "primary",
                }
            )
            # mutate a random subset of pathway genes
            for pid_pw, genes in gene_sets.items():
                for g in genes:
                    if rng.random() < 0.08:
                        edge_rows.append(
                            {
                                "project_id": pid,
                                "case_id": case_id,
                                "gene_symbol": g,
                                "mutated": 1,
                            }
                        )
        for rank, g in enumerate(sorted(kb.all_seed_genes())[:40]):
            top_rows.append(
                {
                    "project_id": pid,
                    "gene_symbol": g,
                    "ssm_occurrence_count": int(1000 - 20 * rank + rng.integers(0, 50)),
                }
            )
    return {
        "clinical": pd.DataFrame(clin_rows),
        "mutation_edges": pd.DataFrame(edge_rows),
        "top_genes": pd.DataFrame(top_rows),
    }


def build_training_corpus(
    *,
    projects: list[str] | None = None,
    out_dir: Path | None = None,
    expand_reactome: bool = True,
    max_projects: int | None = None,
    offline_demo: bool = False,
    demo_cases_per_project: int = 400,
) -> dict[str, Path]:
    """
    Build a multi-million-scale analytical corpus:
      clinical cases × pathway mutation hits × synthetic VOC targets from priors.
    """
    kb = KnowledgeBase()
    out_dir = Path(out_dir or PROCESSED_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    projects = list(projects or DEFAULT_GDC_PROJECTS)
    if max_projects is not None:
        projects = projects[:max_projects]

    pathway_genes = sorted(kb.all_seed_genes())
    gene_sets = kb.pathway_gene_universe()

    if offline_demo:
        print("[corpus] offline demo mode — generating large synthetic multi-project matrix")
        tables = _build_offline_demo_tables(kb, n_cases_per_project=demo_cases_per_project)
    else:
        if expand_reactome:
            try:
                rx = ReactomeClient()
                expanded = rx.expand_pathway_map(kb.pathways)
                expanded_path = out_dir / "pathway_gene_sets.csv"
                expanded.to_csv(expanded_path, index=False)
                gene_sets = {
                    pid: set(g.gene_symbol.tolist())
                    for pid, g in expanded.groupby("pathway_id")
                }
                pathway_genes = sorted({g for s in gene_sets.values() for g in s})
            except Exception as e:  # noqa: BLE001
                print(f"[warn] Reactome expansion skipped: {e}")

        gdc = GDCClient()
        tables = gdc.harvest_projects(
            projects=projects,
            top_n=50,
            include_edges=True,
            pathway_genes=pathway_genes,
        )

    clinical = tables["clinical"].copy()
    edges = tables.get("mutation_edges", pd.DataFrame())
    top_genes = tables["top_genes"]

    clinical["stage_num"] = clinical["ajcc_pathologic_stage"].map(_stage_to_num)
    clinical["age_years"] = clinical["age_at_diagnosis_days"] / 365.25

    # Case-level pathway hit counts / fractions
    if not edges.empty:
        case_gene = edges.drop_duplicates(["case_id", "gene_symbol"])
        pathway_rows = []
        for pid, genes in gene_sets.items():
            sub = case_gene[case_gene["gene_symbol"].isin(genes)]
            hits = sub.groupby("case_id")["gene_symbol"].nunique().rename(f"hits_{pid}")
            pathway_rows.append(hits)
        pathway_mat = pd.concat(pathway_rows, axis=1).fillna(0.0)
        for pid, genes in gene_sets.items():
            n = max(len(genes), 1)
            pathway_mat[f"score_{pid}"] = pathway_mat.get(f"hits_{pid}", 0.0) / n
        feature_df = clinical.merge(pathway_mat, left_on="case_id", right_index=True, how="left")
    else:
        feature_df = clinical.copy()

    score_cols = [c for c in feature_df.columns if c.startswith("score_")]
    for c in score_cols:
        feature_df[c] = feature_df[c].fillna(0.0)

    # Synthetic VOC regression targets from mechanistic priors + noise
    # This creates a large labeled matrix for calibrating the emission model.
    voc_rows = []
    rng = np.random.default_rng(42)
    for _, case in feature_df.iterrows():
        disease = kb.resolve_disease(case["project_id"])
        for voc_id, voc in kb.vocs.items():
            log2fc = float(disease.get("voc_log2fc_prior", {}).get(voc_id, 0.0))
            # pathway contribution
            for pid, p in kb.pathways.items():
                score = float(case.get(f"score_{pid}", 0.0))
                bias = float(disease.get("pathway_bias", {}).get(pid, 1.0))
                coef = float(p.get("voc_effects", {}).get(voc_id, 0.0))
                log2fc += coef * score * bias
            # stage / metastasis style modifiers
            stage = case.get("stage_num")
            if stage is not None and not np.isnan(stage):
                log2fc += 0.08 * (stage - 1.5)
            log2fc += rng.normal(0, 0.15)
            healthy = float(voc["healthy_ppb_median"])
            pred = healthy * (2**log2fc)
            voc_rows.append(
                {
                    "case_id": case["case_id"],
                    "project_id": case["project_id"],
                    "voc_id": voc_id,
                    "log2_fold_change": log2fc,
                    "target_ppb": pred,
                    "healthy_ppb": healthy,
                    "primary_site": case.get("primary_site"),
                    "primary_diagnosis": case.get("primary_diagnosis"),
                    "ajcc_pathologic_stage": case.get("ajcc_pathologic_stage"),
                    "stage_num": case.get("stage_num"),
                }
            )

    voc_df = pd.DataFrame(voc_rows)

    paths = {
        "clinical": out_dir / "clinical_cases.csv",
        "top_genes": out_dir / "top_genes.csv",
        "mutation_edges": out_dir / "mutation_edges.csv",
        "case_features": out_dir / "case_pathway_features.csv",
        "voc_targets": out_dir / "voc_training_targets.csv",
        "manifest": out_dir / "corpus_manifest.json",
    }
    clinical.to_csv(paths["clinical"], index=False)
    top_genes.to_csv(paths["top_genes"], index=False)
    if not edges.empty:
        edges.to_csv(paths["mutation_edges"], index=False)
    feature_df.to_csv(paths["case_features"], index=False)
    voc_df.to_csv(paths["voc_targets"], index=False)

    n_edges = 0 if edges.empty else len(edges)
    n_points = len(voc_df) + len(feature_df) + n_edges + len(top_genes)
    used_projects = sorted(clinical["project_id"].dropna().astype(str).unique().tolist())
    manifest = {
        "projects": used_projects,
        "requested_projects": projects,
        "offline_demo": bool(offline_demo),
        "n_cases": int(len(clinical)),
        "n_mutation_edges": int(n_edges),
        "n_pathway_genes": int(len(pathway_genes)),
        "n_voc_targets": int(len(voc_df)),
        "approx_datapoints": int(n_points),
        "note": (
            "VOC targets are mechanistic-prior synthetic labels for model calibration; "
            "replace/augment with measured breath GC-MS cohorts when available."
        ),
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))
    return paths
