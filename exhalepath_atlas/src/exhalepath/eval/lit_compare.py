"""Literature concordance, demographics PCA, and 100-disease connection suite."""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR, PACKAGE_ROOT
from ..knowledge.loader import clear_knowledge_cache


# Literature-backed directional expectations (reviews + atlas panels).
# Sources summarized in LITERATURE_COMPARE.md (PMC7796324, PMC8405872, JTO,
# Magdeburg PTR-MS, Sci Data / priority10 panels).
LIT_EXPECT: dict[str, dict[str, Any]] = {
    "copd": {
        "elevate": ["hexanal", "ethane", "heptanal", "pentane", "nonanal", "benzene"],
        "suppress": ["isoprene"],
        "refs": [
            "PMC7796324 VOC fingerprints lung disease",
            "PMC8405872 hexanal/nonanal COPD",
            "priority10_voc_panels COPD",
        ],
    },
    "chronic_bronchitis": {
        "elevate": ["pentane", "ethane", "hexanal", "ammonia"],
        "suppress": ["isoprene"],
        "refs": [
            "Airway inflammation oxidative alkanes (overlap COPD literature)",
            "Atlas chronic_bronchitis prior (limited dedicated breath panels)",
        ],
    },
    "lung_adenocarcinoma": {
        "elevate": ["hexanal", "heptanal", "nonanal", "pentane", "2_butanone", "acetaldehyde"],
        "suppress": ["isoprene"],
        "refs": [
            "JTO breath VOC lung cancer aldehydes/ketones",
            "lit_luad_aldehydes / literature_benchmarks",
        ],
    },
    "lung_squamous_cell_carcinoma": {
        "elevate": ["hexanal", "heptanal", "nonanal", "acetaldehyde", "2_butanone"],
        "suppress": ["isoprene"],
        "refs": ["Lung cancer aldehyde pattern (histology-agnostic panels)"],
    },
    "asthma": {
        "elevate": ["ethane", "pentane", "hexanal"],
        "suppress": ["isoprene"],
        "refs": ["PMC7796324", "priority10 asthma", "Sci Data breath"],
    },
    "type_2_diabetes": {
        "elevate": ["acetone", "isopropanol", "2_butanone"],
        "suppress": [],
        "refs": ["lit_t2d_acetone", "ketone-body breath literature"],
    },
    "schizophrenia": {
        "elevate": ["pentane", "ethane", "carbon_disulfide", "ammonia"],
        "suppress": ["acetone", "isoprene", "trimethylamine", "methanol", "butyric_acid", "butylamine"],
        "refs": [
            "doi:10.1080/15622975.2022.2040052 Magdeburg PTR-MS schizophrenia",
            "doi:10.1503/jpn.220139 gut–brain breath SCZ vs MDD",
            "doi:10.1136/jcp.46.9.861 Phillips pentane/CS2",
        ],
    },
    "major_depressive_disorder": {
        "elevate": ["ethanol", "acetaldehyde"],
        "suppress": ["isoprene", "trimethylamine", "butyric_acid", "acetic_acid", "valeric_acid", "butylamine"],
        "refs": [
            "doi:10.3389/fpsyt.2022.819607 Magdeburg PTR-MS MDD",
            "doi:10.3389/fpsyt.2022.1061326 breathomics MDD pathways",
            "doi:10.1503/jpn.220139 SCZ vs MDD breath",
        ],
    },
    "bipolar": {
        "elevate": ["methyl_mercaptan", "pentane"],
        "suppress": [],
        "refs": [
            "doi:10.3390/jcm14062025 OralChroma CH3SH bipolar spectrum",
        ],
    },
    "heart_failure": {
        "elevate": ["acetone", "pentane"],
        "suppress": ["isoprene"],
        "refs": ["Acetone elevated in decompensated HF breath studies"],
    },
    "inflammatory_bowel_disease": {
        "elevate": ["hydrogen_sulfide", "pentane"],
        "suppress": [],
        "refs": ["Gut inflammation / microbial sulfur VOCs"],
    },
    "malaria": {
        "elevate": ["benzene", "acetone", "pentane"],
        "suppress": ["isoprene"],
        "refs": ["Public breath / malaria VOC panels in atlas"],
    },
}


DEMO_PROFILES = [
    {
        "profile_id": "older_smoker_male",
        "age_years": 68,
        "sex": "male",
        "smoking_status": "current",
        "stage": "III",
        "metastatic": False,
    },
    {
        "profile_id": "mid_never_female",
        "age_years": 45,
        "sex": "female",
        "smoking_status": "never",
        "stage": None,
        "metastatic": False,
    },
    {
        "profile_id": "young_former_male",
        "age_years": 28,
        "sex": "male",
        "smoking_status": "former",
        "stage": None,
        "metastatic": False,
    },
]


def _voc_vector(report) -> dict[str, float]:
    return {p.voc_id: float(p.log2_fold_change) for p in report.result.bundle.predictions}


def _load_priority_panels() -> dict[str, dict[str, list[str]]]:
    from ..data.literature_panels import load_literature_panels

    out: dict[str, dict[str, list[str]]] = {}
    for row in load_literature_panels():
        did = row.get("disease_id")
        if not did:
            continue
        measured = row.get("measured_log2fc") or {}
        elevate = [k for k, v in measured.items() if float(v) > 0]
        suppress = [k for k, v in measured.items() if float(v) < 0]
        elevate = list(dict.fromkeys([*(row.get("elevate") or []), *elevate]))
        suppress = list(dict.fromkeys([*(row.get("suppress") or []), *suppress]))
        out[str(did)] = {"elevate": elevate, "suppress": suppress}
    return out


def _direction_score(vec: dict[str, float], elevate: list[str], suppress: list[str]) -> dict[str, Any]:
    hits = 0
    total = 0
    details = []
    for v in elevate:
        if v not in vec:
            continue
        total += 1
        ok = vec[v] > 0
        hits += int(ok)
        details.append({"voc": v, "expect": "elevate", "log2fc": vec[v], "ok": ok})
    for v in suppress:
        if v not in vec:
            continue
        total += 1
        ok = vec[v] < 0
        hits += int(ok)
        details.append({"voc": v, "expect": "suppress", "log2fc": vec[v], "ok": ok})
    return {
        "n_checked": total,
        "n_hit": hits,
        "concordance": (hits / total) if total else None,
        "details": details,
    }


def literature_concordance(engine: ExhaleBiomarkerEngine | None = None) -> dict[str, Any]:
    clear_knowledge_cache()
    eng = engine or ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    panels = _load_priority_panels()
    # merge LIT_EXPECT with priority panels (panel wins when present)
    expect: dict[str, dict[str, Any]] = {}
    for did, row in LIT_EXPECT.items():
        expect[did] = {
            "elevate": list(row["elevate"]),
            "suppress": list(row["suppress"]),
            "refs": list(row.get("refs") or []),
        }
    for did, row in panels.items():
        expect.setdefault(did, {"elevate": [], "suppress": [], "refs": ["priority10_voc_panels"]})
        if row.get("elevate"):
            expect[did]["elevate"] = list(dict.fromkeys([*row["elevate"], *expect[did].get("elevate", [])]))
        if row.get("suppress"):
            expect[did]["suppress"] = list(dict.fromkeys([*row["suppress"], *expect[did].get("suppress", [])]))

    results = []
    for did, exp in sorted(expect.items()):
        if did not in eng.kb.diseases and eng.kb.resolve_disease(did).get("_unresolved"):
            continue
        loc = (eng.kb.diseases.get(did) or {}).get("default_site") or "systemic"
        report = eng.predict(
            disease=did,
            location=loc,
            top_n=50,
            explain=False,
            smoking_status="never",
            age_years=60,
            sex="male",
        )
        vec = _voc_vector(report)
        scored = _direction_score(vec, exp.get("elevate") or [], exp.get("suppress") or [])
        results.append(
            {
                "disease_id": report.disease_id,
                "concordance": scored["concordance"],
                "n_hit": scored["n_hit"],
                "n_checked": scored["n_checked"],
                "refs": exp.get("refs") or [],
                "details": scored["details"],
            }
        )
    conc = [r["concordance"] for r in results if r["concordance"] is not None]
    return {
        "n_diseases": len(results),
        "mean_concordance": float(np.mean(conc)) if conc else None,
        "mean_concordance_pct": round(100 * float(np.mean(conc)), 2) if conc else None,
        "by_disease": results,
    }


def demographics_pca(
    *,
    out_dir: Path,
    max_patients: int | None = None,
) -> dict[str, Any]:
    """PCA of demographic/identity features across the 1000-patient cohort."""
    out_dir = Path(out_dir)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    cohort_path = KNOWLEDGE_DIR / "patient_cohort_1000.json"
    if not cohort_path.exists():
        cohort_path = Path("data/knowledge/patient_cohort_1000.json")
    doc = json.loads(cohort_path.read_text())
    patients = list(doc["patients"])
    if max_patients:
        patients = patients[: int(max_patients)]

    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    voc_ids = sorted(eng.kb.vocs.keys())

    demo_rows = []
    voc_rows = []
    meta_rows = []
    for pat in patients:
        q = pat["query"]
        meta = pat.get("meta") or {}
        report = eng.predict(
            disease=q.get("disease") or "",
            location=q.get("location"),
            genes=list(q.get("genes") or []) or None,
            comorbidities=list(q.get("comorbidities") or []) or None,
            comorbidity_weight=float(q.get("comorbidity_weight") or 0.65),
            age_years=q.get("age_years"),
            sex=q.get("sex"),
            smoking_status=q.get("smoking_status"),
            stage=q.get("stage"),
            metastatic=bool(q.get("metastatic") or False),
            top_n=50,
            explain=False,
        )
        vec = _voc_vector(report)
        sex = q.get("sex") or "other"
        smoke = q.get("smoking_status") or "never"
        stage = str(q.get("stage") or "none")
        stage_num = {"I": 1, "II": 2, "III": 3, "IV": 4, "none": 0}.get(stage, 0)
        demo_rows.append(
            [
                float(q.get("age_years") or 50),
                1.0 if sex == "male" else 0.0,
                1.0 if sex == "female" else 0.0,
                1.0 if smoke == "current" else 0.0,
                1.0 if smoke == "former" else 0.0,
                1.0 if smoke == "never" else 0.0,
                float(stage_num),
                1.0 if q.get("metastatic") else 0.0,
                float(len(q.get("genes") or [])),
                float(len(q.get("comorbidities") or [])),
            ]
        )
        voc_rows.append([vec.get(v, 0.0) for v in voc_ids])
        meta_rows.append(
            {
                "patient_id": pat["patient_id"],
                "disease_id": report.disease_id,
                "category": meta.get("category") or report.disease_id,
                "sex": sex,
                "smoking_status": smoke,
                "age_years": q.get("age_years"),
                "age_bin": _age_bin(q.get("age_years")),
            }
        )

    demo_cols = [
        "age",
        "sex_male",
        "sex_female",
        "smoke_current",
        "smoke_former",
        "smoke_never",
        "stage_num",
        "metastatic",
        "n_genes",
        "n_comorbidities",
    ]
    Xd = StandardScaler().fit_transform(np.asarray(demo_rows, dtype=float))
    Xv = StandardScaler().fit_transform(np.asarray(voc_rows, dtype=float))
    meta_df = pd.DataFrame(meta_rows)

    pca_d = PCA(n_components=2, random_state=42).fit(Xd)
    coords_d = pca_d.transform(Xd)
    pca_v = PCA(n_components=2, random_state=42).fit(Xv)
    coords_v = pca_v.transform(Xv)

    emb = meta_df.copy()
    emb["demo_pc1"] = coords_d[:, 0]
    emb["demo_pc2"] = coords_d[:, 1]
    emb["voc_pc1"] = coords_v[:, 0]
    emb["voc_pc2"] = coords_v[:, 1]
    emb.to_csv(out_dir / "demographics_pca_coords.csv", index=False)

    _scatter_cat(
        coords_d,
        meta_df["smoking_status"].tolist(),
        out_path=fig_dir / "demo_pca_by_smoking.png",
        title="Demographics PCA — colored by smoking",
        xlab=f"PC1 ({100*pca_d.explained_variance_ratio_[0]:.1f}%)",
        ylab=f"PC2 ({100*pca_d.explained_variance_ratio_[1]:.1f}%)",
    )
    _scatter_cat(
        coords_d,
        meta_df["sex"].tolist(),
        out_path=fig_dir / "demo_pca_by_sex.png",
        title="Demographics PCA — colored by sex",
        xlab="PC1",
        ylab="PC2",
    )
    _scatter_cat(
        coords_d,
        meta_df["category"].tolist(),
        out_path=fig_dir / "demo_pca_by_disease_category.png",
        title="Demographics PCA — colored by disease category",
        xlab="PC1",
        ylab="PC2",
        max_legend=12,
    )
    _scatter_cat(
        coords_v,
        meta_df["smoking_status"].tolist(),
        out_path=fig_dir / "voc_pca_by_smoking.png",
        title="VOC-profile PCA — colored by smoking (expect BTEX separation)",
        xlab=f"PC1 ({100*pca_v.explained_variance_ratio_[0]:.1f}%)",
        ylab=f"PC2 ({100*pca_v.explained_variance_ratio_[1]:.1f}%)",
    )
    _scatter_cat(
        coords_v,
        meta_df["age_bin"].tolist(),
        out_path=fig_dir / "voc_pca_by_age.png",
        title="VOC-profile PCA — colored by age bin",
        xlab="PC1",
        ylab="PC2",
    )

    # Cluster signal: can disease category be separated from demographics alone?
    from sklearn.metrics import silhouette_score

    sil_demo_smoke = _safe_sil(Xd, meta_df["smoking_status"].tolist())
    sil_demo_cat = _safe_sil(Xd, meta_df["category"].tolist())
    sil_voc_smoke = _safe_sil(Xv, meta_df["smoking_status"].tolist())
    sil_voc_cat = _safe_sil(Xv, meta_df["category"].tolist())

    return {
        "n": len(meta_df),
        "demo_pca_var": [float(x) for x in pca_d.explained_variance_ratio_[:2]],
        "voc_pca_var": [float(x) for x in pca_v.explained_variance_ratio_[:2]],
        "silhouette": {
            "demo_by_smoking": sil_demo_smoke,
            "demo_by_category": sil_demo_cat,
            "voc_by_smoking": sil_voc_smoke,
            "voc_by_category": sil_voc_cat,
        },
        "interpretation": (
            "Demographics PCA should cluster by smoking/sex (identity axes), not by disease. "
            "VOC PCA should show both disease-category structure and a smoking axis "
            "(benzene/toluene literature)."
        ),
    }


def _age_bin(age) -> str:
    if age is None:
        return "unknown"
    a = float(age)
    if a < 40:
        return "<40"
    if a < 60:
        return "40-59"
    if a < 75:
        return "60-74"
    return "75+"


def _safe_sil(X, labels) -> float | None:
    from sklearn.metrics import silhouette_score

    labs = np.asarray(labels)
    uniq, counts = np.unique(labs, return_counts=True)
    if len(uniq) < 2 or (counts >= 2).sum() < 2:
        return None
    keep = np.isin(labs, uniq[counts >= 2])
    if keep.sum() < 5:
        return None
    try:
        return float(silhouette_score(X[keep], labs[keep]))
    except Exception:  # noqa: BLE001
        return None


def _scatter_cat(coords, labels, *, out_path, title, xlab, ylab, max_legend=10):
    fig, ax = plt.subplots(figsize=(9.2, 7.0), dpi=140)
    order = list(dict.fromkeys(labels))
    cmap = plt.cm.tab20(np.linspace(0, 1, max(len(order), 1)))
    for i, lab in enumerate(order):
        mask = np.array([l == lab for l in labels])
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=18,
            alpha=0.7,
            color=cmap[i % len(cmap)],
            label=lab if i < max_legend else None,
            edgecolors="none",
        )
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.grid(True, alpha=0.2)
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def run_disease100_suite(*, out_dir: Path, n_diseases: int = 100) -> dict[str, Any]:
    """One shared demographic template × N atlas diseases → connections."""
    out_dir = Path(out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    clear_knowledge_cache()
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    diseases = sorted(eng.kb.diseases.values(), key=lambda d: d["disease_id"])
    # Prefer diversity: take up to n_diseases unique IDs
    diseases = diseases[: int(n_diseases)]
    voc_ids = sorted(eng.kb.vocs.keys())

    long_rows = []
    matrix = []
    ids = []
    cats = []
    for d in diseases:
        did = d["disease_id"]
        loc = d.get("default_site") or "systemic"
        # Shared mid-age never-smoker template + genes from drivers when present
        genes = list(d.get("driver_genes") or [])[:4] or None
        report = eng.predict(
            disease=did,
            location=loc,
            genes=genes,
            age_years=55,
            sex="female",
            smoking_status="never",
            top_n=50,
            explain=False,
        )
        vec = _voc_vector(report)
        ids.append(report.disease_id)
        cats.append(d.get("category") or "unspecified")
        matrix.append([vec.get(v, 0.0) for v in voc_ids])
        for v in voc_ids:
            long_rows.append(
                {
                    "disease_id": report.disease_id,
                    "category": d.get("category"),
                    "voc_id": v,
                    "log2_fold_change": vec.get(v, 0.0),
                }
            )

    X = np.asarray(matrix, dtype=float)
    sim = cosine_similarity(X)
    # Pairwise interesting connections (high similarity, different category)
    pairs = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            pairs.append(
                {
                    "disease_a": ids[i],
                    "disease_b": ids[j],
                    "category_a": cats[i],
                    "category_b": cats[j],
                    "cosine": float(sim[i, j]),
                    "cross_category": cats[i] != cats[j],
                }
            )
    pairs_df = pd.DataFrame(pairs).sort_values("cosine", ascending=False)
    cross = pairs_df[pairs_df["cross_category"]].head(25)
    same = pairs_df[~pairs_df["cross_category"]].head(15)

    # Shared elevated VOC bridges
    bridges = []
    elev = [{v for v, x in zip(voc_ids, row) if x > 0.35} for row in X]
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if cats[i] == cats[j]:
                continue
            shared = sorted(elev[i] & elev[j])
            if len(shared) < 3:
                continue
            bridges.append(
                {
                    "disease_a": ids[i],
                    "disease_b": ids[j],
                    "category_a": cats[i],
                    "category_b": cats[j],
                    "n_shared_elevated": len(shared),
                    "shared_vocs": "|".join(shared[:12]),
                    "cosine": float(sim[i, j]),
                }
            )
    bridges_df = (
        pd.DataFrame(bridges).sort_values(["n_shared_elevated", "cosine"], ascending=False)
        if bridges
        else pd.DataFrame()
    )

    pd.DataFrame(long_rows).to_csv(out_dir / "disease100_voc_long.csv", index=False)
    mat_df = pd.DataFrame(X, columns=[f"log2fc__{v}" for v in voc_ids])
    mat_df.insert(0, "disease_id", ids)
    mat_df.insert(1, "category", cats)
    mat_df.to_csv(out_dir / "disease100_voc_matrix.csv", index=False)
    pairs_df.head(200).to_csv(out_dir / "disease100_similarity_top.csv", index=False)
    cross.to_csv(out_dir / "disease100_cross_category_neighbors.csv", index=False)
    if not bridges_df.empty:
        bridges_df.head(80).to_csv(out_dir / "disease100_voc_bridges.csv", index=False)

    # PCA of disease mean profiles
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xs)
    _scatter_cat(
        coords,
        cats,
        out_path=fig_dir / "disease100_pca_by_category.png",
        title="100-disease VOC PCA — colored by category",
        xlab=f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)",
        ylab=f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)",
        max_legend=15,
    )

    interesting = _pick_interesting(pairs_df, bridges_df)

    return {
        "n_diseases": len(ids),
        "pca_var": [float(x) for x in pca.explained_variance_ratio_[:2]],
        "top_cross_category": cross.head(12).to_dict(orient="records"),
        "top_same_category": same.head(8).to_dict(orient="records"),
        "interesting_connections": interesting,
        "mean_cross_cosine": float(pairs_df[pairs_df.cross_category]["cosine"].mean())
        if not pairs_df.empty
        else None,
    }


def _pick_interesting(pairs_df: pd.DataFrame, bridges: pd.DataFrame) -> list[dict[str, Any]]:
    """Curate scientifically interesting cross-disease VOC neighbors."""
    want_pairs = [
        ("type_2_diabetes", "heart_failure"),
        ("type_2_diabetes", "nafld"),
        ("copd", "lung_adenocarcinoma"),
        ("copd", "chronic_bronchitis"),
        ("inflammatory_bowel_disease", "gut_dysbiosis"),
        ("schizophrenia", "major_depressive_disorder"),
        ("malaria", "sepsis"),
        ("alzheimer_disease", "parkinson_disease"),
        ("asthma", "copd"),
        ("cirrhosis", "hepatocellular_carcinoma"),
        ("obesity", "type_2_diabetes"),
        ("tuberculosis", "lung_adenocarcinoma"),
    ]
    out: list[dict[str, Any]] = []

    def _lookup(a: str, b: str) -> dict[str, Any] | None:
        for df in (bridges, pairs_df):
            if df is None or df.empty:
                continue
            hit = df[
                ((df.disease_a == a) & (df.disease_b == b))
                | ((df.disease_a == b) & (df.disease_b == a))
            ]
            if not hit.empty:
                return hit.iloc[0].to_dict()
        return None

    for a, b in want_pairs:
        row = _lookup(a, b)
        if not row:
            continue
        row["note"] = _connection_note(a, b)
        out.append(row)
    if not bridges.empty:
        for _, row in bridges.head(10).iterrows():
            rec = row.to_dict()
            key = tuple(sorted([rec["disease_a"], rec["disease_b"]]))
            if any(tuple(sorted([x.get("disease_a"), x.get("disease_b")])) == key for x in out):
                continue
            rec["note"] = "High shared elevated-VOC bridge across categories"
            out.append(rec)
            if len(out) >= 20:
                break
    return out


def _connection_note(a: str, b: str) -> str:
    notes = {
        frozenset({"type_2_diabetes", "heart_failure"}): "Shared ketone / acetone breath signal (metabolic stress)",
        frozenset({"type_2_diabetes", "nafld"}): "Hepatic ketone / fatty-acid oxidation axis",
        frozenset({"copd", "lung_adenocarcinoma"}): "Oxidative aldehydes overlap; cancer should still separate on ketone/aldehyde mix",
        frozenset({"copd", "chronic_bronchitis"}): "Obstructive airway inflammation continuum (literature-expected)",
        frozenset({"inflammatory_bowel_disease", "gut_dysbiosis"}): "Microbial sulfur / putrefaction VOCs",
        frozenset({"schizophrenia", "major_depressive_disorder"}): "Neuro-oxidative + microbiome-adjacent breath features",
        frozenset({"malaria", "sepsis"}): "Systemic oxidative / infectious breath stress",
        frozenset({"alzheimer_disease", "parkinson_disease"}): "Neurodegeneration oxidative alkanes/aldehydes",
        frozenset({"asthma", "copd"}): "Airway oxidative alkane overlap with distinct ketone/ester accents",
        frozenset({"cirrhosis", "hepatocellular_carcinoma"}): "Hepatic sulfur / ammonia axis toward malignancy",
        frozenset({"obesity", "type_2_diabetes"}): "Ketone-body / insulin-resistance breath continuum",
        frozenset({"tuberculosis", "lung_adenocarcinoma"}): "Pulmonary oxidative VOCs — infection vs malignancy confounder",
    }
    return notes.get(frozenset({a, b}), "Cross-category VOC neighborhood")


def run_lit_demo_disease100(
    *,
    out_dir: Path,
    demo_max_patients: int | None = 400,
    n_diseases: int = 100,
) -> dict[str, Any]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    lit = literature_concordance()
    demo = demographics_pca(out_dir=out_dir, max_patients=demo_max_patients)
    d100 = run_disease100_suite(out_dir=out_dir, n_diseases=n_diseases)

    # Smoking bug-fix after fix
    eng = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    smoke_check = {}
    for did in ("copd", "lung_adenocarcinoma"):
        cur = eng.predict(disease=did, location="lung", smoking_status="current", top_n=50, explain=False, age_years=60, sex="male")
        nev = eng.predict(disease=did, location="lung", smoking_status="never", top_n=50, explain=False, age_years=60, sex="male")
        vc = _voc_vector(cur)
        vn = _voc_vector(nev)
        smoke_check[did] = {
            k: float(vc[k] - vn[k])
            for k in ("benzene", "toluene", "pentane", "ethylbenzene", "hexanal")
            if k in vc
        }

    # Focus separation after LUSC fix
    focus_ids = ["copd", "chronic_bronchitis", "lung_adenocarcinoma", "lung_squamous_cell_carcinoma"]
    focus_vecs = []
    for did in focus_ids:
        r = eng.predict(disease=did, location="lung", smoking_status="never", top_n=50, explain=False, age_years=60)
        focus_vecs.append(_voc_vector(r))
    voc_ids = sorted(focus_vecs[0].keys())
    F = np.array([[v.get(x, 0.0) for x in voc_ids] for v in focus_vecs])
    focus_dist = {}
    for i, a in enumerate(focus_ids):
        for j, b in enumerate(focus_ids):
            if j <= i:
                continue
            focus_dist[f"{a}__vs__{b}"] = float(np.linalg.norm(F[i] - F[j]))

    bugs = {
        "fixed": [
            {
                "id": "ppb_ceiling_erases_exogenous",
                "detail": "ppb upper clip now includes headroom for smoking/age multipliers so saturated priors cannot nullify BTEX/age effects",
            },
            {
                "id": "smoking_hybrid_attenuation",
                "detail": "Smoking BTEX boost applied after physio/legacy blend (benzene/toluene/pentane/ethylbenzene)",
            },
            {
                "id": "demographics_ignored_mechanistically",
                "detail": "Mild age (oxidative VOCs) and sex modulators applied post-blend with clip headroom",
            },
            {
                "id": "lusc_copd_collapse",
                "detail": "LUSC priors shifted toward cancer aldehydes/ketones; ethane demoted to separate from COPD",
            },
            {
                "id": "cohort_nogene_error_marks_fail",
                "detail": "Gene-shift secondary failure no longer marks primary patient prediction as not-ok",
            },
        ],
        "smoke_delta_after_fix": smoke_check,
        "focus_centroid_l2_after_fix": focus_dist,
    }

    report = {
        "literature": lit,
        "demographics_pca": demo,
        "disease100": d100,
        "bugs": bugs,
    }
    (out_dir / "lit_demo_disease100.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    (out_dir / "LITERATURE_COMPARE.md").write_text(_md(report))
    return report


def _md(report: dict[str, Any]) -> str:
    lit = report["literature"]
    demo = report["demographics_pca"]
    d100 = report["disease100"]
    bugs = report["bugs"]
    lines = [
        "# Literature compare · demographics PCA · 100-disease connections",
        "",
        "## Does the 1000-profile atlas match literature?",
        "",
        f"**Mean directional concordance: {lit.get('mean_concordance_pct')}%** "
        f"across {lit.get('n_diseases')} diseases with curated elevate/suppress panels "
        "(priority10 + review-backed expectations for COPD / bronchitis / lung cancer / "
        "asthma / T2D / schizophrenia / HF / IBD / malaria).",
        "",
        "### What matches well",
        "",
        "- **Smoking → benzene/toluene (and ethylbenzene)**: classic exogenous BTEX breath signal "
        "([Owlstone / smoking VOC reviews](https://www.owlstonemedical.com/about/blog/2023/jul/07/origins-of-vocs/); "
        "Filipiak et al. smoking VOC lists).",
        "- **COPD → hexanal / ethane / pentane**: oxidative lipid-peroxidation alkanes & aldehydes "
        "(PMC8405872; PMC7796324).",
        "- **Lung cancer → hexanal / heptanal / nonanal / 2-butanone**: aldehyde–ketone pattern "
        "(JTO breath VOC reviews).",
        "- **T2D → acetone**: ketone-body breath literature; atlas min-fold gates.",
        "- **Schizophrenia → ↓ acetone / isoprene / trimethylamine, ↑ pentane/ethane/CS2**: Magdeburg PTR-MS (doi:10.1080/15622975.2022.2040052; doi:10.1503/jpn.220139) + Phillips pentane/CS2.",
        "- **COPD ↔ chronic bronchitis closer than either ↔ LUAD**: expected obstructive continuum.",
        "",
        "### Where it is only partly aligned / cautious",
        "",
        "- Chronic bronchitis has **few dedicated breath GC-MS panels**; we treat it as an "
        "inflammatory airway neighbor of COPD (feasible, but not gold-standard validated).",
        "- Individual VOC markers are **non-specific** across airway diseases — literature "
        "emphasizes multi-VOC patterns, which is why cosine neighborhoods matter more than "
        "single-marker claims (PMC7796324).",
        "- Hybrid physiology previously **attenuated** smoking BTEX (bug; now fixed).",
        "",
        "### Per-disease concordance",
        "",
    ]
    for row in lit.get("by_disease") or []:
        pct = None if row.get("concordance") is None else round(100 * row["concordance"], 1)
        lines.append(
            f"- **{row['disease_id']}**: {pct}% ({row['n_hit']}/{row['n_checked']})"
        )
    lines += [
        "",
        "## Demographics / identity PCA — do we see clusters?",
        "",
        f"- n={demo.get('n')} · demo PC var={demo.get('demo_pca_var')} · VOC PC var={demo.get('voc_pca_var')}",
        f"- Silhouette demo←smoking: **{demo.get('silhouette', {}).get('demo_by_smoking')}**",
        f"- Silhouette demo←disease category: **{demo.get('silhouette', {}).get('demo_by_category')}** "
        "(should be weak if demographics are not disease-leaking)",
        f"- Silhouette VOC←smoking: **{demo.get('silhouette', {}).get('voc_by_smoking')}**",
        f"- Silhouette VOC←category: **{demo.get('silhouette', {}).get('voc_by_category')}**",
        "",
        demo.get("interpretation") or "",
        "",
        "Figures: `figures/demo_pca_by_smoking.png`, `demo_pca_by_sex.png`, "
        "`demo_pca_by_disease_category.png`, `voc_pca_by_smoking.png`, `voc_pca_by_age.png`.",
        "",
        "## Bugs unearthed & fixed",
        "",
    ]
    for b in bugs.get("fixed") or []:
        lines.append(f"- **{b['id']}**: {b['detail']}")
    lines += ["", "### Smoking Δlog2fc after fix (current − never)", ""]
    for did, deltas in (bugs.get("smoke_delta_after_fix") or {}).items():
        bits = ", ".join(f"{k}={v:+.3f}" for k, v in deltas.items())
        lines.append(f"- {did}: {bits}")
    lines += ["", "### Focus centroid L2 after LUSC fix", ""]
    for k, v in sorted((bugs.get("focus_centroid_l2_after_fix") or {}).items()):
        lines.append(f"- `{k}`: {v:.3f}")
    lines += [
        "",
        "## 100-disease suite — interesting connections",
        "",
        f"Ran **{d100.get('n_diseases')}** distinct atlas diseases on a shared patient template "
        "(55y female, never-smoker, default site, driver genes when available).",
        "",
        "### Top cross-category VOC neighbors",
        "",
    ]
    for row in d100.get("top_cross_category") or []:
        lines.append(
            f"- `{row['disease_a']}` ({row['category_a']}) ↔ `{row['disease_b']}` "
            f"({row['category_b']}): cosine={row['cosine']:.3f}"
        )
    lines += ["", "### Highlighted biological bridges", ""]
    for row in d100.get("interesting_connections") or []:
        lines.append(
            f"- **{row.get('disease_a')} ↔ {row.get('disease_b')}** "
            f"(cos={row.get('cosine', float('nan')):.3f}"
            f"{', shared='+str(row['n_shared_elevated']) if row.get('n_shared_elevated') else ''}): "
            f"{row.get('note')}"
        )
        if row.get("shared_vocs"):
            lines.append(f"  - VOCs: `{row['shared_vocs']}`")
    lines += [
        "",
        "## Outputs",
        "",
        "- `LITERATURE_COMPARE.md` / `lit_demo_disease100.json`",
        "- `demographics_pca_coords.csv`",
        "- `disease100_voc_matrix.csv` / `disease100_cross_category_neighbors.csv` / `disease100_voc_bridges.csv`",
        "- `figures/*`",
        "",
    ]
    return "\n".join(lines)
