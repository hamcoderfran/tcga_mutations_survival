"""
Project 1: TCGA mutation + survival data download
====================================================

Downloads, for each of four TCGA cancer cohorts (LUAD, BRCA, PAAD, COAD):
  1. The most frequently mutated genes (by number of simple somatic
     mutation occurrences) via the GDC ssm_occurrences aggregation endpoint.
  2. For the top mutated genes, the list of case IDs that carry a mutation
     in that gene (via paginated ssm_occurrences queries).
  3. Clinical / survival data (vital status, days to death, days to last
     follow-up, age, sex, tumor stage) for every case in the cohort, via
     the GDC cases endpoint.

All data comes directly from the official NCI Genomic Data Commons (GDC)
public REST API (https://api.gdc.cancer.gov) -- no authentication, no file
downloads of arbitrary binaries, JSON only.

Outputs (per project), written to ../data/<project>/:
  - top_genes.csv             top mutated genes with occurrence counts
  - mutated_cases_<GENE>.csv  case IDs carrying a mutation in GENE
  - clinical.csv              per-case clinical/survival data
"""

import json
import os
import time

import pandas as pd
import requests

GDC_API = "https://api.gdc.cancer.gov"
PROJECTS = ["TCGA-LUAD", "TCGA-BRCA", "TCGA-PAAD", "TCGA-COAD"]

# Number of top mutated genes to pull per-case mutation status for
N_TOP_GENES = 10

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")

session = requests.Session()
session.headers.update({"User-Agent": "tcga-mutation-survival-research/1.0"})


def gdc_get(endpoint, params, retries=5):
    url = f"{GDC_API}/{endpoint}"
    for attempt in range(retries):
        try:
            r = session.get(url, params=params, timeout=60)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, json.JSONDecodeError) as e:
            wait = 2 ** attempt
            print(f"  [retry {attempt+1}] {endpoint} failed: {e} -> sleeping {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"GDC API call to {endpoint} failed after {retries} retries")


def top_mutated_genes(project_id, top_n=25):
    """Top genes by number of ssm_occurrences (mutation x case rows) in a project."""
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "case.project.project_id", "value": [project_id]}}
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "size": 0,
        "facets": "ssm.consequence.transcript.gene.symbol",
        "facet_size": 200,
    }
    d = gdc_get("ssm_occurrences", params)
    buckets = d["data"]["aggregations"]["ssm.consequence.transcript.gene.symbol"]["buckets"]
    buckets = sorted(buckets, key=lambda b: b["doc_count"], reverse=True)[:top_n]
    return pd.DataFrame(
        [{"gene_symbol": b["key"].upper(), "ssm_occurrence_count": b["doc_count"]} for b in buckets]
    )


def cases_with_mutation(project_id, gene_symbol):
    """Return the set of case_ids in `project_id` carrying >=1 SSM in `gene_symbol`."""
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "case.project.project_id", "value": [project_id]}},
            {
                "op": "in",
                "content": {
                    "field": "ssm.consequence.transcript.gene.symbol",
                    "value": [gene_symbol.lower()],
                },
            },
        ],
    }
    case_ids = set()
    size = 100
    frm = 0
    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": "case.case_id",
            "size": size,
            "from": frm,
        }
        d = gdc_get("ssm_occurrences", params)
        hits = d["data"]["hits"]
        for h in hits:
            cid = h.get("case", {}).get("case_id")
            if cid:
                case_ids.add(cid)
        total = d["data"]["pagination"]["total"]
        frm += size
        if frm >= total:
            break
    return case_ids


def clinical_data(project_id):
    """Pull per-case clinical / survival fields for every case in a project."""
    filters = {
        "op": "in",
        "content": {"field": "project.project_id", "value": [project_id]},
    }
    fields = ",".join(
        [
            "case_id",
            "submitter_id",
            "demographic.vital_status",
            "demographic.days_to_death",
            "demographic.gender",
            "demographic.race",
            "demographic.ethnicity",
            "diagnoses.age_at_diagnosis",
            "diagnoses.days_to_last_follow_up",
            "diagnoses.ajcc_pathologic_stage",
            "diagnoses.primary_diagnosis",
            "diagnoses.tumor_grade",
        ]
    )
    rows = []
    size = 100
    frm = 0
    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": fields,
            "size": size,
            "from": frm,
        }
        d = gdc_get("cases", params)
        hits = d["data"]["hits"]
        for h in hits:
            demo = h.get("demographic", {}) or {}
            diags = h.get("diagnoses", []) or [{}]
            diag = diags[0] if diags else {}
            rows.append(
                {
                    "case_id": h.get("case_id"),
                    "submitter_id": h.get("submitter_id"),
                    "vital_status": demo.get("vital_status"),
                    "days_to_death": demo.get("days_to_death"),
                    "gender": demo.get("gender"),
                    "race": demo.get("race"),
                    "ethnicity": demo.get("ethnicity"),
                    "age_at_diagnosis_days": diag.get("age_at_diagnosis"),
                    "days_to_last_follow_up": diag.get("days_to_last_follow_up"),
                    "ajcc_pathologic_stage": diag.get("ajcc_pathologic_stage"),
                    "primary_diagnosis": diag.get("primary_diagnosis"),
                    "tumor_grade": diag.get("tumor_grade"),
                }
            )
        total = d["data"]["pagination"]["total"]
        frm += size
        if frm >= total:
            break
    return pd.DataFrame(rows)


def main():
    for project_id in PROJECTS:
        print(f"\n=== {project_id} ===")
        proj_dir = os.path.join(DATA_DIR, project_id)
        os.makedirs(proj_dir, exist_ok=True)

        print("  fetching top mutated genes...")
        top_genes = top_mutated_genes(project_id, top_n=25)
        top_genes.to_csv(os.path.join(proj_dir, "top_genes.csv"), index=False)
        print(top_genes.head(N_TOP_GENES).to_string(index=False))

        print(f"  fetching mutated-case lists for top {N_TOP_GENES} genes...")
        for gene in top_genes["gene_symbol"].head(N_TOP_GENES):
            cases = cases_with_mutation(project_id, gene)
            out = os.path.join(proj_dir, f"mutated_cases_{gene}.csv")
            pd.DataFrame({"case_id": sorted(cases)}).to_csv(out, index=False)
            print(f"    {gene}: {len(cases)} cases")

        print("  fetching clinical/survival data...")
        clin = clinical_data(project_id)
        clin.to_csv(os.path.join(proj_dir, "clinical.csv"), index=False)
        print(f"    {len(clin)} cases total")

    print("\nDone.")


if __name__ == "__main__":
    main()
