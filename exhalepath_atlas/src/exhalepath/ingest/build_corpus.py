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
    DEPRECATED for VOC label generation.

    The historical path built *synthetic* VOC targets (mechanistic prior + noise).
    Production training now uses ``build_real_training_corpus`` only.

    Passing ``offline_demo=True`` raises unless EXPLICITLY overridden via
    env ``VOC_ALLOW_SYNTHETIC=1`` (CI emergency only).
    """
    import os

    if offline_demo and os.environ.get("VOC_ALLOW_SYNTHETIC") != "1":
        raise RuntimeError(
            "Synthetic VOC training corpora are disabled. "
            "Use: python -m exhalepath build-real-corpus && python -m exhalepath train\n"
            "Or set VOC_ALLOW_SYNTHETIC=1 only for emergency/CI debugging."
        )

    kb = KnowledgeBase()
    out_dir = Path(out_dir or PROCESSED_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Prefer real corpus builder when available
    from .real_breath_corpus import build_real_training_corpus

    print("[corpus] redirecting to REAL breath corpus (no synthetic VOC labels)")
    return build_real_training_corpus(out_dir=out_dir, remove_synthetic=True)
