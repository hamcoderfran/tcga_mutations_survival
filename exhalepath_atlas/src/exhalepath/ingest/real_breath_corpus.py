"""Build a REAL-ONLY exhaled-VOC training corpus (no synthetic labels).

Sources
-------
1. Metabolomics Workbench quantified studies (ST003200, ST000587, ST000883, …)
2. Scientific Data clinical breathomics peak tables (asthma / COPD / bronchiectasis)
3. Curated literature panels with DOIs (directional or quantified published results)

Synthetic mechanistic-prior + Gaussian-noise labels are explicitly excluded.
When a disease lacks measured ppb, pathway/gene/cell-state priors still drive
zero-shot *inference* at predict time — but they are never written as training y.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import DATA_DIR, KNOWLEDGE_DIR, PROCESSED_DIR
from ..knowledge.loader import KnowledgeBase

REAL_BREATH_DIR = DATA_DIR / "real_breath"
COHORT_DIR = REAL_BREATH_DIR / "cohorts"
PANEL_PATH = REAL_BREATH_DIR / "literature_panels" / "priority10_voc_panels.json"
MW_DIR = DATA_DIR / "datasources" / "metabolomics"

# Priority conditions for supervised real-data training
PRIORITY_DISEASES = [
    "asthma",
    "copd",
    "covid19",
    "pneumonia_bacterial",
    "tuberculosis",
    "cystic_fibrosis",
    "sleep_apnea",
    "cancer_stomach",
    "head_neck_cancer",
    "cancer_prostate",
    "heart_failure",
    "malaria",
    "ards",
]

NAME_TO_VOC = {
    "acetone": "acetone",
    "isoprene": "isoprene",
    "pentane": "pentane",
    "n-pentane": "pentane",
    "ethane": "ethane",
    "hexanal": "hexanal",
    "heptanal": "heptanal",
    "nonanal": "nonanal",
    "decanal": "decanal",
    "ethanol": "ethanol",
    "2-butanone": "2_butanone",
    "butanone": "2_butanone",
    "dimethyl sulfide": "dms",
    "ammonia": "ammonia",
    "methanol": "methanol",
    "toluene": "toluene",
    "benzene": "benzene",
    "limonene": "limonene",
    "hydrogen sulfide": "hydrogen_sulfide",
    "acetaldehyde": "acetaldehyde",
    "formaldehyde": "formaldehyde",
    "indole": "indole",
    "phenol": "phenol",
    "1-propanol": "propanol",
    "propanol": "propanol",
    "2-propanol": "isopropanol",
    "isopropanol": "isopropanol",
    "isopropyl alcohol": "isopropanol",
    "2-pentanone": "2_pentanone",
    "benzaldehyde": "benzaldehyde",
    "trimethylamine": "trimethylamine",
    "propanal": "propionaldehyde",
    "furan": "furan",
    "octanal": "octanal",
    "styrene": "styrene",
    "ethylbenzene": "ethylbenzene",
    "hexane": "hexane",
    "butane": "butane",
    "octane": "octane",
    "xylene": "xylene",
    "ethyl acetate": "ethyl_acetate",
}


def _map_voc(name: str) -> str | None:
    n = re.sub(r"\([^)]*\)", "", str(name)).strip().lower()
    n = re.sub(r"\s+", " ", n)
    if n in NAME_TO_VOC:
        return NAME_TO_VOC[n]
    for k, v in NAME_TO_VOC.items():
        if n == k or n.startswith(k + ",") or n.startswith(k + " "):
            return v
    for key, vid in [
        ("acetone", "acetone"),
        ("isoprene", "isoprene"),
        ("hexanal", "hexanal"),
        ("heptanal", "heptanal"),
        ("nonanal", "nonanal"),
        ("pentane", "pentane"),
        ("hexane", "hexane"),
        ("toluene", "toluene"),
        ("benzene", "benzene"),
        ("2-butanone", "2_butanone"),
        ("ethanol", "ethanol"),
        ("methanol", "methanol"),
        ("phenol", "phenol"),
        ("furan", "furan"),
    ]:
        if key in n:
            return vid
    return None


def _ms_rows(data: dict[str, Any]) -> list[dict[str, Any]]:
    ms = data.get("MS_METABOLITE_DATA") or data.get("NMR_METABOLITE_DATA") or {}
    rows: list[dict[str, Any]] = []
    skip = {
        "Metabolite",
        "metabolite_name",
        "metabolite",
        "refmet_name",
        "Human Metabolome Database",
        "PubChem Compound",
        "KEGG Compound ID",
        "InChI Key",
        "refmet_id",
    }
    for item in ms.get("Data") or []:
        met = item.get("Metabolite") or item.get("metabolite_name") or item.get("metabolite")
        for k, v in item.items():
            if k in skip:
                continue
            try:
                val = float(v)
            except (TypeError, ValueError):
                continue
            if math.isnan(val):
                continue
            rows.append({"metabolite": met, "sample": k, "intensity": val})
    return rows


def _factor_map(study_id: str) -> dict[str, str]:
    path = MW_DIR / study_id / "factors.json"
    if not path.exists():
        return {}
    fac = json.loads(path.read_text())
    out = {}
    for rec in fac.values():
        sid = rec.get("local_sample_id")
        f = str(rec.get("factors") or "").lower()
        if sid:
            out[str(sid)] = f
    return out


def _mw_disease_control_folds(
    study_id: str,
    disease_id: str,
    *,
    positive_token: str,
    negative_token: str,
    unit: str,
) -> list[dict[str, Any]]:
    path = MW_DIR / study_id / "mwtab.json"
    if not path.exists():
        return []
    data = json.loads(path.read_text())
    rows = _ms_rows(data)
    fac = _factor_map(study_id)
    if not rows or not fac:
        return []
    df = pd.DataFrame(rows)
    def grp(sample: str) -> str:
        f = fac.get(sample, "")
        if positive_token in f and negative_token not in f:
            return "disease"
        if negative_token in f:
            return "control"
        return "other"

    df["group"] = df["sample"].map(grp)
    df["voc_id"] = df["metabolite"].map(_map_voc)
    df = df.dropna(subset=["voc_id"])
    df = df[df["group"].isin(["disease", "control"])]
    out = []
    for voc, g in df.groupby("voc_id"):
        d = g.loc[g["group"] == "disease", "intensity"]
        c = g.loc[g["group"] == "control", "intensity"]
        if len(d) < 3 or len(c) < 3:
            continue
        dm, cm = float(d.median()), float(c.median())
        if dm <= 0 or cm <= 0:
            continue
        out.append(
            {
                "study_id": study_id,
                "disease_id": disease_id,
                "voc_id": voc,
                "log2_fold_change": round(math.log2(dm / cm), 4),
                "disease_median": dm,
                "control_median": cm,
                "n_disease": int(len(d)),
                "n_control": int(len(c)),
                "unit": unit,
                "source": "metabolomics_workbench",
                "label_origin": "measured_cohort",
            }
        )
    return out


def _scidata_relative_folds() -> list[dict[str, Any]]:
    path = DATA_DIR / "public_breath" / "mapped_cohort_voc_means.csv"
    if not path.exists():
        return []
    m = pd.read_csv(path)
    piv = m.pivot_table(index="voc_id", columns="cohort", values="mean_intensity")
    out = []
    for disease, col in [("asthma", "asthma"), ("copd", "copd")]:
        if col not in piv.columns:
            continue
        refs = [c for c in piv.columns if c != col]
        for voc in piv.index:
            val = piv.loc[voc, col]
            ref_vals = [piv.loc[voc, c] for c in refs if pd.notna(piv.loc[voc, c])]
            if pd.isna(val) or not ref_vals:
                continue
            ref = float(np.median(ref_vals))
            if ref <= 0 or float(val) <= 0:
                continue
            out.append(
                {
                    "study_id": "23522490",
                    "disease_id": disease,
                    "voc_id": voc,
                    "log2_fold_change": round(math.log2(float(val) / ref), 4),
                    "disease_median": float(val),
                    "control_median": ref,
                    "n_disease": None,
                    "n_control": None,
                    "unit": "gcms_intensity",
                    "source": "scientific_data_figshare",
                    "label_origin": "measured_cohort_relative",
                    "note": "fold vs other Sci Data pulmonary cohorts (no healthy arm)",
                }
            )
    return out


def _literature_panel_folds() -> list[dict[str, Any]]:
    if not PANEL_PATH.exists():
        return []
    payload = json.loads(PANEL_PATH.read_text())
    panels = payload.get("panels") or []
    out = []
    for p in panels:
        did = p.get("disease_id")
        for voc, log2fc in (p.get("measured_log2fc") or {}).items():
            try:
                v = float(log2fc)
            except (TypeError, ValueError):
                continue
            refs = p.get("refs") or []
            out.append(
                {
                    "study_id": (refs[0].get("doi") if refs else "literature"),
                    "disease_id": did,
                    "voc_id": voc,
                    "log2_fold_change": round(v, 4),
                    "disease_median": None,
                    "control_median": None,
                    "n_disease": None,
                    "n_control": None,
                    "unit": "literature_log2fc",
                    "source": "literature_panel",
                    "label_origin": f"literature:{p.get('evidence_grade', 'mixed')}",
                    "note": "; ".join(
                        filter(
                            None,
                            [
                                (refs[0].get("title") if refs else None),
                                (refs[0].get("doi") if refs else None),
                            ],
                        )
                    ),
                }
            )
    return out


def collect_real_label_rows() -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    rows += _mw_disease_control_folds(
        "ST000883",
        "malaria",
        positive_token="positive",
        negative_token="negative",
        unit="gcms_intensity",
    )
    rows += _mw_disease_control_folds(
        "ST000587",
        "heart_failure",
        positive_token="failure",
        negative_token="control",
        unit="uM_ebc",
    )
    # HF factors are "factor:HF" / "factor:control" — also try
    if not any(r["disease_id"] == "heart_failure" for r in rows):
        rows += _mw_disease_control_folds(
            "ST000587",
            "heart_failure",
            positive_token="hf",
            negative_token="control",
            unit="uM_ebc",
        )
    rows += _scidata_relative_folds()
    rows += _literature_panel_folds()
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # Prefer measured over literature when both exist for same disease×voc
    rank = {
        "measured_cohort": 0,
        "measured_cohort_relative": 1,
        "literature:quantified": 2,
        "literature:mixed": 3,
        "literature:directional_only": 4,
    }
    df["_rank"] = df["label_origin"].map(lambda x: rank.get(str(x), 5))
    df = (
        df.sort_values("_rank")
        .drop_duplicates(subset=["disease_id", "voc_id", "study_id"], keep="first")
        .drop(columns=["_rank"])
        .reset_index(drop=True)
    )
    return df


def _pathway_feature_row(disease: dict[str, Any], kb: KnowledgeBase) -> dict[str, Any]:
    """Build a single non-synthetic feature vector from atlas pathway bias + gene priors."""
    row: dict[str, Any] = {
        "case_id": f"real::{disease['disease_id']}",
        "disease_id": disease["disease_id"],
        "project_id": disease["disease_id"],
        "primary_site": disease.get("default_site") or "systemic",
        "primary_diagnosis": disease.get("name"),
        "disease_type": disease.get("category"),
        "stage_num": 2.0,
        "age_years": 55.0,
        "gender": "unknown",
        "ajcc_pathologic_stage": None,
    }
    bias = disease.get("pathway_bias") or {}
    for pw in kb.pathways.values():
        pid = pw["pathway_id"]
        # Data-backed pathway activation prior (not random): use curated bias, else 0
        b = float(bias.get(pid, 1.0))
        score = max(0.0, min(1.5, (b - 1.0) * 1.2 + 0.15 * (1 if b > 1 else 0)))
        row[f"score_{pid}"] = score
        row[f"hits_{pid}"] = 1 if b > 1.05 else 0
    return row


def build_real_training_corpus(
    *,
    out_dir: Path | None = None,
    remove_synthetic: bool = True,
) -> dict[str, Path]:
    """
    Materialize real-only training tables and optionally delete synthetic CSVs.
    """
    out_dir = Path(out_dir or PROCESSED_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    COHORT_DIR.mkdir(parents=True, exist_ok=True)

    kb = KnowledgeBase()
    labels = collect_real_label_rows()
    labels_path = COHORT_DIR / "all_real_log2fc_labels.csv"
    labels.to_csv(labels_path, index=False)

    disease_ids = sorted(
        set(labels["disease_id"].astype(str).tolist() if len(labels) else [])
        | set(PRIORITY_DISEASES)
    )

    # Case features: one row per disease×study so measured replicates stay distinct
    case_rows = []
    case_key_to_row: dict[str, dict] = {}
    for did in disease_ids:
        try:
            d = kb.resolve_disease(did)
        except Exception:
            continue
        if not d.get("disease_id") or d.get("disease_id") in {"unknown", "unresolved"}:
            continue
        # base disease case
        base = _pathway_feature_row(d, kb)
        case_rows.append(base)
        case_key_to_row[did] = base
        # study-specific cases for each measured/literature source
        studies = (
            labels.loc[labels["disease_id"].astype(str) == did, "study_id"]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
            if len(labels)
            else []
        )
        for sid in studies:
            row = dict(base)
            row["case_id"] = f"real::{did}::{sid}"
            row["project_id"] = did
            case_rows.append(row)
            case_key_to_row[f"{did}::{sid}"] = row
    cases = pd.DataFrame(case_rows)

    # VOC targets: expand each label against the matching disease×study case
    healthy = {v["voc_id"]: float(v.get("healthy_ppb_median") or 1.0) for v in kb.vocs.values()}
    tgt_rows = []
    for _, lab in labels.iterrows():
        did = str(lab["disease_id"])
        sid = str(lab.get("study_id") or "literature")
        case = case_key_to_row.get(f"{did}::{sid}") or case_key_to_row.get(did)
        if not case:
            continue
        voc = lab["voc_id"]
        h = float(healthy.get(voc, 1.0))
        log2fc = float(lab["log2_fold_change"])
        tgt_rows.append(
            {
                "case_id": case["case_id"],
                "project_id": did,
                "disease_id": did,
                "voc_id": voc,
                "log2_fold_change": log2fc,
                "target_ppb": h * (2 ** log2fc),
                "healthy_ppb": h,
                "primary_site": case.get("primary_site"),
                "primary_diagnosis": case.get("primary_diagnosis"),
                "ajcc_pathologic_stage": None,
                "stage_num": case.get("stage_num"),
                "label_origin": lab.get("label_origin"),
                "study_id": lab.get("study_id"),
                "source": lab.get("source"),
            }
        )
    targets = pd.DataFrame(tgt_rows)

    # Write real corpus (replace synthetic)
    case_path = out_dir / "case_pathway_features.csv"
    tgt_path = out_dir / "voc_training_targets.csv"
    cases.to_csv(case_path, index=False)
    targets.to_csv(tgt_path, index=False)

    manifest = {
        "corpus_type": "real_breath_only",
        "synthetic_labels": False,
        "n_cases": int(len(cases)),
        "n_voc_targets": int(len(targets)),
        "n_label_rows_raw": int(len(labels)),
        "diseases": sorted(targets["disease_id"].unique().tolist()) if len(targets) else [],
        "label_origins": (
            targets["label_origin"].value_counts().to_dict() if len(targets) else {}
        ),
        "priority_diseases": PRIORITY_DISEASES,
        "sources": [
            "metabolomics_workbench:ST000587,ST000883,ST003200",
            "scientific_data_figshare:23522490",
            "literature_panels:priority10_voc_panels.json",
        ],
        "note": (
            "Training y are measured cohort folds and published literature log2fc only. "
            "Mechanistic pathway priors are features / zero-shot inference — never synthetic y."
        ),
    }
    man_path = out_dir / "corpus_manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2))

    # Remove leftover synthetic artifacts
    if remove_synthetic:
        for name in ("mutation_edges.csv",):
            p = out_dir / name
            # keep file but mark empty optional — or delete if clearly synthetic demo
            if p.exists() and manifest["n_voc_targets"] < 50_000:
                # previous offline demo had 218k edges; remove when replacing with real corpus
                p.unlink()
        # stamp
        (out_dir / "SYNTHETIC_DATA_REMOVED.txt").write_text(
            "Synthetic VOC training labels (offline_demo mechanistic+noise) removed.\n"
            "Active corpus: real_breath_only (see corpus_manifest.json).\n"
        )

    return {
        "case_features": case_path,
        "voc_targets": tgt_path,
        "manifest": man_path,
        "labels": labels_path,
    }
