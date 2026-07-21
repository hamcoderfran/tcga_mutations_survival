from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from ..config import CACHE_DIR, DEFAULT_GDC_PROJECTS, GDC_API
from .http import request_json


class GDCClient:
    """High-throughput GDC pulls for mutation, clinical, and disease context."""

    def __init__(self, cache_dir: Path | None = None):
        self.cache_dir = Path(cache_dir or CACHE_DIR / "gdc")
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def top_mutated_genes(self, project_id: str, top_n: int = 50) -> pd.DataFrame:
        cache = self.cache_dir / f"{project_id}_top_genes.csv"
        if cache.exists():
            return pd.read_csv(cache)

        filters = {
            "op": "and",
            "content": [
                {
                    "op": "in",
                    "content": {"field": "case.project.project_id", "value": [project_id]},
                }
            ],
        }
        params = {
            "filters": json.dumps(filters),
            "size": 0,
            "facets": "ssm.consequence.transcript.gene.symbol",
            "facet_size": max(top_n, 200),
        }
        data = request_json("GET", f"{GDC_API}/ssm_occurrences", params=params)
        buckets = data["data"]["aggregations"]["ssm.consequence.transcript.gene.symbol"][
            "buckets"
        ]
        buckets = sorted(buckets, key=lambda b: b["doc_count"], reverse=True)[:top_n]
        df = pd.DataFrame(
            [
                {
                    "project_id": project_id,
                    "gene_symbol": b["key"].upper(),
                    "ssm_occurrence_count": b["doc_count"],
                }
                for b in buckets
            ]
        )
        df.to_csv(cache, index=False)
        return df

    def clinical_cases(self, project_id: str, max_cases: int = 5000) -> pd.DataFrame:
        cache = self.cache_dir / f"{project_id}_clinical.csv"
        if cache.exists():
            return pd.read_csv(cache)

        filters = {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "project.project_id", "value": [project_id]}}
            ],
        }
        fields = [
            "case_id",
            "submitter_id",
            "primary_site",
            "disease_type",
            "demographic.gender",
            "demographic.vital_status",
            "demographic.days_to_death",
            "diagnoses.age_at_diagnosis",
            "diagnoses.ajcc_pathologic_stage",
            "diagnoses.tumor_stage",
            "diagnoses.primary_diagnosis",
            "diagnoses.site_of_resection_or_biopsy",
            "diagnoses.tissue_or_organ_of_origin",
            "diagnoses.morphology",
            "diagnoses.classification_of_tumor",
        ]
        rows: list[dict[str, Any]] = []
        size = 100
        from_ = 0
        while from_ < max_cases:
            params = {
                "filters": json.dumps(filters),
                "fields": ",".join(fields),
                "format": "JSON",
                "size": size,
                "from": from_,
            }
            data = request_json("GET", f"{GDC_API}/cases", params=params)
            hits = data.get("data", {}).get("hits", [])
            if not hits:
                break
            for h in hits:
                demo = h.get("demographic") or {}
                diags = h.get("diagnoses") or [{}]
                d0 = diags[0] if diags else {}
                rows.append(
                    {
                        "project_id": project_id,
                        "case_id": h.get("case_id"),
                        "submitter_id": h.get("submitter_id"),
                        "primary_site": h.get("primary_site"),
                        "disease_type": h.get("disease_type"),
                        "gender": demo.get("gender"),
                        "vital_status": demo.get("vital_status"),
                        "days_to_death": demo.get("days_to_death"),
                        "age_at_diagnosis_days": d0.get("age_at_diagnosis"),
                        "ajcc_pathologic_stage": d0.get("ajcc_pathologic_stage")
                        or d0.get("tumor_stage"),
                        "primary_diagnosis": d0.get("primary_diagnosis"),
                        "tissue_or_organ_of_origin": d0.get("tissue_or_organ_of_origin"),
                        "site_of_resection_or_biopsy": d0.get("site_of_resection_or_biopsy"),
                        "morphology": d0.get("morphology"),
                        "classification_of_tumor": d0.get("classification_of_tumor"),
                    }
                )
            from_ += size
            if len(hits) < size:
                break

        df = pd.DataFrame(rows)
        df.to_csv(cache, index=False)
        return df

    def mutation_case_gene_edges(
        self,
        project_id: str,
        genes: Iterable[str],
        page_size: int = 2000,
        max_rows: int = 500_000,
        gene_batch_size: int = 40,
    ) -> pd.DataFrame:
        """Pull case↔gene mutation edges for pathway genes (scales to large matrices)."""
        genes = sorted({g.upper() for g in genes})
        cache_csv = self.cache_dir / f"{project_id}_mut_edges.csv"
        if cache_csv.exists():
            return pd.read_csv(cache_csv)

        rows: list[dict[str, Any]] = []
        for i in range(0, len(genes), gene_batch_size):
            batch = genes[i : i + gene_batch_size]
            filters = {
                "op": "and",
                "content": [
                    {
                        "op": "in",
                        "content": {
                            "field": "case.project.project_id",
                            "value": [project_id],
                        },
                    },
                    {
                        "op": "in",
                        "content": {
                            "field": "ssm.consequence.transcript.gene.symbol",
                            "value": batch,
                        },
                    },
                ],
            }
            from_ = 0
            while from_ < max_rows:
                params = {
                    "filters": json.dumps(filters),
                    "fields": "case.case_id,ssm.consequence.transcript.gene.symbol",
                    "format": "JSON",
                    "size": page_size,
                    "from": from_,
                }
                data = request_json("GET", f"{GDC_API}/ssm_occurrences", params=params)
                hits = data.get("data", {}).get("hits", [])
                if not hits:
                    break
                for h in hits:
                    case = (h.get("case") or {}).get("case_id")
                    cons = h.get("ssm", {}).get("consequence") or []
                    gene_sym = None
                    for c in cons:
                        gene_sym = (
                            ((c.get("transcript") or {}).get("gene") or {}).get("symbol")
                        )
                        if gene_sym:
                            gene_sym = gene_sym.upper()
                            break
                    if not gene_sym or gene_sym not in batch:
                        # fall back: keep row tagged with first batch gene only if parse fails
                        continue
                    rows.append(
                        {
                            "project_id": project_id,
                            "case_id": case,
                            "gene_symbol": gene_sym,
                            "mutated": 1,
                        }
                    )
                from_ += page_size
                if len(hits) < page_size:
                    break
                if len(rows) >= max_rows:
                    break

        df = pd.DataFrame(rows).drop_duplicates()
        df.to_csv(cache_csv, index=False)
        return df

    def harvest_projects(
        self,
        projects: list[str] | None = None,
        top_n: int = 40,
        include_edges: bool = True,
        pathway_genes: Iterable[str] | None = None,
    ) -> dict[str, pd.DataFrame]:
        projects = projects or DEFAULT_GDC_PROJECTS
        genes = list(pathway_genes) if pathway_genes is not None else []
        top_frames = []
        clin_frames = []
        edge_frames = []
        for pid in projects:
            print(f"[GDC] {pid}: top genes + clinical")
            top_frames.append(self.top_mutated_genes(pid, top_n=top_n))
            clin_frames.append(self.clinical_cases(pid))
            if include_edges and genes:
                print(f"[GDC] {pid}: mutation edges for {len(genes)} pathway genes")
                edge_frames.append(self.mutation_case_gene_edges(pid, genes))

        out = {
            "top_genes": pd.concat(top_frames, ignore_index=True),
            "clinical": pd.concat(clin_frames, ignore_index=True),
        }
        if edge_frames:
            out["mutation_edges"] = pd.concat(edge_frames, ignore_index=True)
        return out
