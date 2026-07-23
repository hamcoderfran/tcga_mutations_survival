"""
Download and normalize public exhaled-breath VOC datasets for external validation.

Primary open dataset
--------------------
Kuo et al. Scientific Data 2024 — clinical breathomics (asthma / COPD / bronchiectasis)
Figshare: https://doi.org/10.6084/m9.figshare.23522490.v6
Paper:    https://doi.org/10.1038/s41597-024-03052-2

Peak tables are within-disease GC-MS intensities (no healthy controls). We derive
cross-disease differential markers (e.g. higher in Asthma than COPD/Bronchiectasis)
and map compounds onto the ExhalePath 50-VOC catalog for top-k ranking tests.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import DATA_DIR, KNOWLEDGE_DIR
from ..knowledge.loader import KnowledgeBase

PUBLIC_BREATH_DIR = DATA_DIR / "public_breath"

FIGSHARE_FILES = {
    "Asthma_peaktable_ver3.csv": "https://ndownloader.figshare.com/files/43994397",
    "COPD_peaktable_ver3.csv": "https://ndownloader.figshare.com/files/43994391",
    "Bronchi_peaktable_ver3.csv": "https://ndownloader.figshare.com/files/43994394",
    "CBD_metadata_for_ver3.xlsx": "https://ndownloader.figshare.com/files/43994400",
    "intersection_of_detected_compunds.xlsx": "https://ndownloader.figshare.com/files/44414426",
}

# PubChem CIDs for ExhalePath catalog VOCs (where known)
VOC_PUBCHEM_CID: dict[str, int] = {
    "acetone": 180,
    "isoprene": 6557,
    "pentane": 8003,
    "ethane": 6324,
    "hexanal": 6184,
    "heptanal": 8130,
    "nonanal": 31289,
    "decanal": 8175,
    "ethanol": 702,
    "2_butanone": 6560,
    "dms": 1068,
    "ammonia": 222,
    "methanol": 887,
    "toluene": 1140,
    "benzene": 241,
    "limonene": 22311,
    "hydrogen_sulfide": 402,
    "acetaldehyde": 177,
    "formaldehyde": 712,
    "indole": 798,
    "phenol": 996,
    "dimethyl_disulfide": 12232,
    "propanol": 1031,
    "isopropanol": 3776,
    "butane": 7843,
    "hexane": 8058,
    "octane": 356,
    "styrene": 7501,
    "ethylbenzene": 7500,
    "xylene": 7809,
    "cyclohexane": 8078,
    "methyl_acetate": 6584,
    "ethyl_acetate": 8857,
    "acetonitrile": 6342,
    "furan": 8029,
    "2_pentanone": 7895,
    "3_methylbutanal": 11552,
    "benzaldehyde": 240,
    "octanal": 454,
    "undecane": 14257,
    "dodecane": 8182,
    "carbon_disulfide": 6348,
    "methyl_mercaptan": 878,
    "trimethylamine": 1146,
    "dimethyl_amine": 674,
    "pyrrole": 8027,
    "propionaldehyde": 527,
    "crotonaldehyde": 447466,
}


def _norm_name(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def download_scientific_data_breathomics(out_dir: Path | None = None) -> dict[str, Path]:
    """Fetch Figshare peak tables for the Scientific Data clinical breathomics set."""
    from .secure_fetch import secure_fetch

    out_dir = Path(out_dir or PUBLIC_BREATH_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, url in FIGSHARE_FILES.items():
        dest = out_dir / name
        if dest.exists() and dest.stat().st_size > 1000:
            paths[name] = dest
            continue
        secure_fetch(
            url,
            dest=dest,
            quarantine_dir=out_dir / ".quarantine",
            max_bytes=80 * 1024 * 1024,
            timeout=120,
        )
        paths[name] = dest
    return paths


def _load_peak_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    # normalize column names
    cols = {c: str(c).strip() for c in df.columns}
    df = df.rename(columns=cols)
    name_col = "IUPAC Name" if "IUPAC Name" in df.columns else df.columns[1]
    cid_col = "pubchem_CID" if "pubchem_CID" in df.columns else df.columns[0]
    sample_cols = [c for c in df.columns if c not in {name_col, cid_col}]
    long = df.melt(
        id_vars=[cid_col, name_col],
        value_vars=sample_cols,
        var_name="sample_id",
        value_name="intensity",
    )
    long = long.rename(columns={cid_col: "pubchem_cid", name_col: "iupac_name"})
    long["intensity"] = pd.to_numeric(long["intensity"], errors="coerce")
    long = long.dropna(subset=["intensity"])
    return long


def map_compound_to_voc(
    *,
    pubchem_cid: Any,
    iupac_name: str,
    kb: KnowledgeBase | None = None,
) -> str | None:
    """Map a peak-table compound to an ExhalePath voc_id when possible."""
    kb = kb or KnowledgeBase()
    try:
        cid = int(float(pubchem_cid))
    except (TypeError, ValueError):
        cid = None
    if cid is not None:
        for voc_id, pcid in VOC_PUBCHEM_CID.items():
            if pcid == cid and voc_id in kb.vocs:
                return voc_id
    name = _norm_name(iupac_name)
    # direct / substring match against catalog names and ids
    for voc_id, voc in kb.vocs.items():
        aliases = {_norm_name(voc_id), _norm_name(voc.get("name") or "")}
        for a in aliases:
            if not a:
                continue
            if name == a or a in name or name in a:
                return voc_id
    # common IUPAC → catalog
    aliases_extra = {
        "benzaldehyde": "benzaldehyde",
        "hexanal": "hexanal",
        "octanal": "octanal",
        "nonanal": "nonanal",
        "propan 2 ol": "isopropanol",
        "propan 1 ol": "propanol",
        "butan 2 one": "2_butanone",
        "pentan 2 one": "2_pentanone",
        "1 phenylethanone": None,
        "styrene": "styrene",
        "ethylbenzene": "ethylbenzene",
        "methylbenzene": "toluene",
        "phenol": "phenol",
        "indole": "indole",
    }
    for key, voc_id in aliases_extra.items():
        if key in name and voc_id and voc_id in kb.vocs:
            return voc_id
    return None


def build_public_breath_benchmark(
    *,
    out_dir: Path | None = None,
    fold_threshold: float = 1.25,
) -> dict[str, Path]:
    """
    Build cross-disease differential VOC expectations from Scientific Data peak tables
    + curated literature public-validation cases.
    """
    out_dir = Path(out_dir or PUBLIC_BREATH_DIR)
    paths = download_scientific_data_breathomics(out_dir)
    kb = KnowledgeBase()

    cohort_files = {
        "asthma": paths["Asthma_peaktable_ver3.csv"],
        "copd": paths["COPD_peaktable_ver3.csv"],
        "bronchiectasis": paths["Bronchi_peaktable_ver3.csv"],
    }
    # disease_id used by ExhalePath atlas
    cohort_to_disease = {
        "asthma": "asthma",
        "copd": "copd",
        "bronchiectasis": "copd",  # atlas proxy (no dedicated bronchiectasis prior yet)
    }

    means = {}
    mapped_rows = []
    for cohort, path in cohort_files.items():
        long = _load_peak_table(path)
        long["cohort"] = cohort
        long["voc_id"] = [
            map_compound_to_voc(pubchem_cid=r.pubchem_cid, iupac_name=r.iupac_name, kb=kb)
            for r in long.itertuples()
        ]
        g = (
            long.dropna(subset=["voc_id"])
            .groupby("voc_id")["intensity"]
            .mean()
        )
        means[cohort] = g
        for voc_id, inten in g.items():
            mapped_rows.append(
                {"cohort": cohort, "voc_id": voc_id, "mean_intensity": float(inten)}
            )

    mean_df = pd.DataFrame(mapped_rows)
    cases = []
    for cohort, g in means.items():
        others = [means[c] for c in means if c != cohort]
        if not others:
            continue
        other = pd.concat(others, axis=1).mean(axis=1)
        elevated = []
        suppressed = []
        for voc_id, val in g.items():
            base = float(other.get(voc_id, np.nan))
            if not np.isfinite(base) or base <= 0:
                continue
            ratio = float(val) / base
            if ratio >= fold_threshold:
                elevated.append(voc_id)
            elif ratio <= 1.0 / fold_threshold:
                suppressed.append(voc_id)
        did = cohort_to_disease[cohort]
        # prefer atlas resolve
        resolved = kb.resolve_disease(did)
        cases.append(
            {
                "case_id": f"sci_data_2024_{cohort}",
                "source": "scientific_data_2024_figshare_23522490",
                "doi": "10.1038/s41597-024-03052-2",
                "cohort": cohort,
                "disease": resolved.get("name") or did,
                "disease_id": resolved.get("disease_id") or did,
                "location": "lung",
                "expect_elevated": sorted(set(elevated)),
                "expect_suppressed": sorted(set(suppressed)),
                "n_mapped_vocs": int(len(g)),
                "notes": (
                    "Cross-cohort differential from GC-MS peak intensities "
                    "(no healthy controls in source tables)."
                ),
            }
        )

    # Curated open-literature directional cases (publicly reported panels)
    curated = [
        {
            "case_id": "lit_t2d_acetone",
            "source": "literature_public",
            "doi": "10.1016/j.jchromb.2012.12.008",
            "disease": "type 2 diabetes",
            "disease_id": "type_2_diabetes",
            "location": "systemic",
            "expect_elevated": ["acetone", "isopropanol", "2_butanone"],
            "expect_suppressed": [],
            "notes": "Classic ketosis / breath acetone axis",
        },
        {
            "case_id": "lit_luad_aldehydes",
            "source": "literature_public",
            "doi": "10.1038/bjc.2011.277",
            "disease": "lung adenocarcinoma",
            "disease_id": "lung_adenocarcinoma",
            "location": "lung",
            "expect_elevated": ["hexanal", "heptanal", "nonanal", "pentane", "2_butanone"],
            "expect_suppressed": ["isoprene"],
            "notes": "Lung cancer breath aldehyde / alkane panel",
        },
        {
            "case_id": "lit_ad_oxidative",
            "source": "literature_public",
            "doi": "10.1016/j.neurobiolaging.2014.06.015",
            "disease": "Alzheimer disease",
            "disease_id": "alzheimer_disease",
            "location": "brain",
            "expect_elevated": ["hexanal", "pentane", "ethane"],
            "expect_suppressed": [],
            "notes": "Neurodegeneration oxidative-stress VOCs",
        },
        {
            "case_id": "lit_ibd_microbiome",
            "source": "literature_public",
            "doi": "10.1136/gutjnl-2013-305799",
            "disease": "inflammatory bowel disease",
            "disease_id": "inflammatory_bowel_disease",
            "location": "gut",
            "expect_elevated": ["indole", "phenol", "hydrogen_sulfide", "trimethylamine"],
            "expect_suppressed": [],
            "notes": "Gut dysbiosis volatiles in IBD",
        },
    ]
    for c in curated:
        c["n_mapped_vocs"] = len(c["expect_elevated"]) + len(c["expect_suppressed"])
        cases.append(c)

    bench = {
        "version": "1.0.0",
        "description": (
            "Public exhaled-breath validation cases: Scientific Data 2024 clinical "
            "breathomics cross-cohort differentials + curated open literature panels."
        ),
        "cases": cases,
    }
    bench_path = KNOWLEDGE_DIR / "public_breath_benchmarks.json"
    mapped_path = out_dir / "mapped_cohort_voc_means.csv"
    bench_path.write_text(json.dumps(bench, indent=2))
    mean_df.to_csv(mapped_path, index=False)
    manifest = {
        "n_cases": len(cases),
        "n_sci_data_cases": sum(1 for c in cases if c["source"].startswith("scientific_data")),
        "paths": {"benchmark": str(bench_path), "mapped_means": str(mapped_path)},
    }
    man_path = out_dir / "public_breath_manifest.json"
    man_path.write_text(json.dumps(manifest, indent=2))
    return {"benchmark": bench_path, "mapped_means": mapped_path, "manifest": man_path}
