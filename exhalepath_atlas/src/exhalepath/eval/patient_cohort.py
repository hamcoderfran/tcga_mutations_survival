"""1000-patient diversity cohort evaluation: dense CSV + VOC/genetic PCA–UMAP."""

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
from sklearn.metrics import pairwise_distances, silhouette_score
from sklearn.preprocessing import StandardScaler

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR
from ..knowledge.loader import clear_knowledge_cache


def _load_cohort(path: Path | None = None) -> dict[str, Any]:
    for p in [
        path,
        KNOWLEDGE_DIR / "patient_cohort_1000.json",
        Path("data/knowledge/patient_cohort_1000.json"),
    ]:
        if p and Path(p).exists():
            return json.loads(Path(p).read_text())
    raise FileNotFoundError(
        "Missing patient_cohort_1000.json — run scripts/build_patient_cohort_1000.py"
    )


def _finite(v: float) -> bool:
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _voc_maps(report) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    log2: dict[str, float] = {}
    delta: dict[str, float] = {}
    fold: dict[str, float] = {}
    for p in report.result.bundle.predictions:
        log2[p.voc_id] = float(p.log2_fold_change)
        delta[p.voc_id] = float(p.delta_ppb)
        fold[p.voc_id] = float(p.fold_change)
    return log2, delta, fold


def _pathway_map(report) -> dict[str, float]:
    return {p.pathway_id: float(p.score) for p in report.result.bundle.pathway_scores}


def _predict_pair(engine: ExhaleBiomarkerEngine, q: dict[str, Any]) -> tuple[Any, Any | None, str | None]:
    """Predict with genes; optionally without genes for genetic-shift vector."""
    kw = dict(
        disease=q.get("disease") or "",
        location=q.get("location"),
        top_n=50,
        mode="hybrid",
        explain=False,
        genes=list(q.get("genes") or []) or None,
        comorbidities=list(q.get("comorbidities") or []) or None,
        comorbidity_weight=float(q.get("comorbidity_weight") or 0.65),
        age_years=q.get("age_years"),
        sex=q.get("sex"),
        smoking_status=q.get("smoking_status"),
        stage=q.get("stage"),
        metastatic=bool(q.get("metastatic") or False),
    )
    try:
        with_genes = engine.predict(**kw)
    except Exception as e:  # noqa: BLE001
        return None, None, f"{type(e).__name__}: {e}"

    no_genes = None
    if kw.get("genes"):
        try:
            kw2 = dict(kw)
            kw2["genes"] = None
            no_genes = engine.predict(**kw2)
        except Exception as e:  # noqa: BLE001
            return with_genes, None, f"nogene:{type(e).__name__}: {e}"
    return with_genes, no_genes, None


def _resource_inventory(engine: ExhaleBiomarkerEngine) -> dict[str, Any]:
    kb = engine.kb
    chembl = bool(getattr(kb, "chembl_priors", None))
    return {
        "n_diseases": len(kb.diseases),
        "n_vocs": len(kb.vocs),
        "n_pathways": len(kb.pathways),
        "n_tissues": len(kb.tissues),
        "n_cell_states": len(kb.cell_states),
        "n_pathway_chains": len(kb.pathway_chains),
        "datasources": sorted(kb.datasources.keys()),
        "chembl_loaded": chembl,
        "atlas_capability": bool(getattr(kb, "atlas_capability", None)),
        "mode": "hybrid",
        "opentargets_live": False,
        "opentargets_offline_priors": "opentargets" in kb.datasources,
        "gdc_offline_priors": "gdc" in kb.datasources,
    }


def run_cohort(
    *,
    out_dir: Path,
    cohort_path: Path | None = None,
    max_patients: int | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    clear_knowledge_cache()
    doc = _load_cohort(cohort_path)
    patients = list(doc["patients"])
    if max_patients is not None:
        patients = patients[: int(max_patients)]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    resources = _resource_inventory(engine)
    voc_ids = sorted(engine.kb.vocs.keys())
    pathway_ids = sorted(engine.kb.pathways.keys())

    patient_rows: list[dict[str, Any]] = []
    long_rows: list[dict[str, Any]] = []
    voc_matrix: list[list[float]] = []
    gene_shift_matrix: list[list[float]] = []
    pathway_matrix: list[list[float]] = []
    breakdown_flags: list[dict[str, Any]] = []
    pid_order: list[str] = []
    labels_disease: list[str] = []
    labels_group: list[str] = []
    has_gene_shift: list[bool] = []

    focus = set(doc.get("focus_triad") or ["copd", "chronic_bronchitis", "lung_adenocarcinoma"])

    for pat in patients:
        pid = pat["patient_id"]
        q = dict(pat.get("query") or {})
        meta = dict(pat.get("meta") or {})
        report, report_ng, err = _predict_pair(engine, q)

        row: dict[str, Any] = {
            "patient_id": pid,
            "stratum": pat.get("stratum"),
            "site_match": pat.get("site_match"),
            "disease_query": q.get("disease"),
            "location_query": q.get("location"),
            "genes": "|".join(q.get("genes") or []),
            "n_genes": len(q.get("genes") or []),
            "comorbidities": "|".join(q.get("comorbidities") or []),
            "n_comorbidities": len(q.get("comorbidities") or []),
            "comorbidity_weight": q.get("comorbidity_weight"),
            "age_years": q.get("age_years"),
            "sex": q.get("sex"),
            "smoking_status": q.get("smoking_status"),
            "stage": q.get("stage"),
            "metastatic": bool(q.get("metastatic") or False),
            "category": meta.get("category"),
            "error": err,
            "ok": err is None and report is not None,
        }

        if report is None:
            patient_rows.append(row)
            breakdown_flags.append(
                {"patient_id": pid, "flag": "exception", "detail": err, "disease": q.get("disease")}
            )
            continue

        did = report.disease_id
        unresolved = did.startswith("custom::")
        loc = report.location or {}
        log2, delta, fold = _voc_maps(report)
        pw = _pathway_map(report)

        max_abs = max((abs(v) for v in log2.values()), default=0.0)
        n_nonfinite = sum(1 for v in log2.values() if not _finite(v))
        n_nonpos = sum(
            1
            for p in report.result.bundle.predictions
            if not (p.predicted_ppb > 0 and p.healthy_ppb > 0)
        )

        row.update(
            {
                "resolved_disease_id": did,
                "resolved_disease_name": report.disease_name,
                "unresolved": unresolved,
                "location_tissue_id": loc.get("tissue_id"),
                "location_matched": bool(loc.get("matched")),
                "n_vocs_modeled": report.n_vocs_modeled,
                "model_version": report.model_version,
                "max_abs_log2fc": max_abs,
                "mean_abs_log2fc": float(np.mean([abs(v) for v in log2.values()])) if log2 else 0.0,
                "top_voc": max(delta, key=lambda k: abs(delta[k])) if delta else None,
                "top_voc_delta_ppb": max((abs(v) for v in delta.values()), default=0.0),
                "n_nonfinite_voc": n_nonfinite,
                "n_nonpositive_ppb": n_nonpos,
                "focus_triad": did in focus or (q.get("disease") in focus),
            }
        )

        # Breakdown heuristics
        if unresolved:
            breakdown_flags.append(
                {"patient_id": pid, "flag": "unresolved_disease", "detail": did, "disease": q.get("disease")}
            )
        if n_nonfinite:
            breakdown_flags.append(
                {"patient_id": pid, "flag": "nonfinite_voc", "detail": n_nonfinite, "disease": did}
            )
        if n_nonpos:
            breakdown_flags.append(
                {"patient_id": pid, "flag": "nonpositive_ppb", "detail": n_nonpos, "disease": did}
            )
        if unresolved and max_abs > 0.45:
            breakdown_flags.append(
                {
                    "patient_id": pid,
                    "flag": "unresolved_not_near_healthy",
                    "detail": max_abs,
                    "disease": did,
                }
            )
        if not loc.get("matched") and pat.get("site_match") == "matched":
            breakdown_flags.append(
                {
                    "patient_id": pid,
                    "flag": "expected_site_unmatched",
                    "detail": loc.get("tissue_id"),
                    "disease": did,
                }
            )

        # Dense long VOC rows
        for vid in voc_ids:
            long_rows.append(
                {
                    "patient_id": pid,
                    "disease_id": did,
                    "category": meta.get("category"),
                    "smoking_status": q.get("smoking_status"),
                    "stage": q.get("stage"),
                    "voc_id": vid,
                    "log2_fold_change": log2.get(vid, 0.0),
                    "delta_ppb": delta.get(vid, 0.0),
                    "fold_change": fold.get(vid, 1.0),
                }
            )

        voc_vec = [log2.get(vid, 0.0) for vid in voc_ids]
        pw_vec = [pw.get(pid_, 0.0) for pid_ in pathway_ids]

        # Genetic shift: with-genes − without-genes log2fc
        if report_ng is not None:
            log2_ng, _, _ = _voc_maps(report_ng)
            shift = [log2.get(vid, 0.0) - log2_ng.get(vid, 0.0) for vid in voc_ids]
            has_shift = True
        else:
            shift = [0.0] * len(voc_ids)
            has_shift = False

        shift_l2 = float(np.linalg.norm(shift))
        row["gene_shift_l2"] = shift_l2
        row["has_gene_shift"] = has_shift
        if has_shift and shift_l2 < 1e-9 and row["n_genes"] > 0:
            breakdown_flags.append(
                {
                    "patient_id": pid,
                    "flag": "genes_no_effect",
                    "detail": "|".join(q.get("genes") or []),
                    "disease": did,
                }
            )

        patient_rows.append(row)
        pid_order.append(pid)
        labels_disease.append(did)
        # Group label for pulmonary triad vs other
        if did in focus or did == "lung_squamous_cell_carcinoma":
            labels_group.append(did)
        else:
            labels_group.append("other")
        voc_matrix.append(voc_vec)
        gene_shift_matrix.append(shift)
        pathway_matrix.append(pw_vec)
        has_gene_shift.append(has_shift)

    # --- DataFrames / dense CSVs ---
    patients_df = pd.DataFrame(patient_rows)
    long_df = pd.DataFrame(long_rows)
    break_df = pd.DataFrame(breakdown_flags)

    voc_df = pd.DataFrame(voc_matrix, columns=[f"log2fc__{v}" for v in voc_ids])
    voc_df.insert(0, "patient_id", pid_order)
    voc_df.insert(1, "disease_id", labels_disease)
    voc_df.insert(2, "group", labels_group)

    shift_df = pd.DataFrame(gene_shift_matrix, columns=[f"dlog2fc__{v}" for v in voc_ids])
    shift_df.insert(0, "patient_id", pid_order)
    shift_df.insert(1, "disease_id", labels_disease)
    shift_df.insert(2, "group", labels_group)
    shift_df.insert(3, "has_gene_shift", has_gene_shift)

    pw_df = pd.DataFrame(pathway_matrix, columns=[f"pw__{p}" for p in pathway_ids])
    pw_df.insert(0, "patient_id", pid_order)
    pw_df.insert(1, "disease_id", labels_disease)

    patients_df.to_csv(out_dir / "patients.csv", index=False)
    long_df.to_csv(out_dir / "patient_voc_long.csv", index=False)
    voc_df.to_csv(out_dir / "patient_voc_matrix.csv", index=False)
    shift_df.to_csv(out_dir / "patient_gene_shift_matrix.csv", index=False)
    pw_df.to_csv(out_dir / "patient_pathway_matrix.csv", index=False)
    break_df.to_csv(out_dir / "breakdown_flags.csv", index=False)

    # Wide dense mega-table: covariates + all VOC log2fc + pathway scores
    dense = patients_df.merge(voc_df.drop(columns=["disease_id", "group"], errors="ignore"), on="patient_id", how="left")
    dense = dense.merge(pw_df.drop(columns=["disease_id"], errors="ignore"), on="patient_id", how="left")
    dense.to_csv(out_dir / "patients_dense.csv", index=False)

    # --- Embeddings ---
    embed_summary: dict[str, Any] = {}
    X = np.asarray(voc_matrix, dtype=float)
    if len(X) >= 5:
        embed_summary["voc"] = _embed_and_plot(
            X,
            labels=labels_group,
            disease_ids=labels_disease,
            patient_ids=pid_order,
            out_dir=fig_dir,
            prefix="voc",
            title="VOC log2FC profiles (PCA / UMAP)",
            focus=focus | {"lung_squamous_cell_carcinoma"},
        )
        # Write embedding coords CSV
        coords = pd.DataFrame(
            {
                "patient_id": pid_order,
                "disease_id": labels_disease,
                "group": labels_group,
                "pca1": embed_summary["voc"]["pca_coords"][:, 0],
                "pca2": embed_summary["voc"]["pca_coords"][:, 1],
            }
        )
        if embed_summary["voc"].get("umap_coords") is not None:
            coords["umap1"] = embed_summary["voc"]["umap_coords"][:, 0]
            coords["umap2"] = embed_summary["voc"]["umap_coords"][:, 1]
        coords.to_csv(out_dir / "embeddings_voc.csv", index=False)

    Xs = np.asarray(gene_shift_matrix, dtype=float)
    mask = np.array(has_gene_shift, dtype=bool)
    if mask.sum() >= 5:
        embed_summary["gene_shift"] = _embed_and_plot(
            Xs[mask],
            labels=[labels_group[i] for i, m in enumerate(mask) if m],
            disease_ids=[labels_disease[i] for i, m in enumerate(mask) if m],
            patient_ids=[pid_order[i] for i, m in enumerate(mask) if m],
            out_dir=fig_dir,
            prefix="gene_shift",
            title="Genetic shift → VOC Δlog2FC (with genes − without)",
            focus=focus | {"lung_squamous_cell_carcinoma"},
        )
        gcoords = pd.DataFrame(
            {
                "patient_id": [pid_order[i] for i, m in enumerate(mask) if m],
                "disease_id": [labels_disease[i] for i, m in enumerate(mask) if m],
                "group": [labels_group[i] for i, m in enumerate(mask) if m],
                "pca1": embed_summary["gene_shift"]["pca_coords"][:, 0],
                "pca2": embed_summary["gene_shift"]["pca_coords"][:, 1],
            }
        )
        if embed_summary["gene_shift"].get("umap_coords") is not None:
            gcoords["umap1"] = embed_summary["gene_shift"]["umap_coords"][:, 0]
            gcoords["umap2"] = embed_summary["gene_shift"]["umap_coords"][:, 1]
        gcoords.to_csv(out_dir / "embeddings_gene_shift.csv", index=False)

    # Focus triad alignment distances
    alignment = _focus_alignment(X, labels_disease, focus | {"lung_squamous_cell_carcinoma"})
    smoking_effect = _smoking_effect(patients_df, voc_df, voc_ids)

    overall = _summarize(
        patients_df,
        break_df,
        resources,
        embed_summary,
        alignment,
        smoking_effect,
        n=len(patients),
    )
    report = {
        "overall": overall,
        "resources": resources,
        "alignment": alignment,
        "smoking_effect": smoking_effect,
        "embed": {
            k: {kk: vv for kk, vv in v.items() if kk not in {"pca_coords", "umap_coords"}}
            for k, v in embed_summary.items()
        },
        "cohort_version": doc.get("version"),
        "n_patients_run": len(patients),
    }
    (out_dir / "cohort_eval.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    (out_dir / "COHORT_1000.md").write_text(_md(report, break_df, patients_df))
    return report


def _embed_and_plot(
    X: np.ndarray,
    *,
    labels: list[str],
    disease_ids: list[str],
    patient_ids: list[str],
    out_dir: Path,
    prefix: str,
    title: str,
    focus: set[str],
) -> dict[str, Any]:
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    pca = PCA(n_components=min(10, Xs.shape[0], Xs.shape[1]), random_state=42)
    pcs = pca.fit_transform(Xs)
    pca_coords = pcs[:, :2]

    umap_coords = None
    umap_error = None
    try:
        import umap  # type: ignore

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reducer = umap.UMAP(
                n_neighbors=min(30, max(5, len(Xs) // 10)),
                min_dist=0.15,
                metric="euclidean",
                random_state=42,
            )
            umap_coords = reducer.fit_transform(Xs)
    except Exception as e:  # noqa: BLE001
        umap_error = f"{type(e).__name__}: {e}"

    # Silhouette on focus+other labels (need ≥2 clusters with ≥2 samples)
    sil = None
    uniq, counts = np.unique(labels, return_counts=True)
    if len(uniq) >= 2 and (counts >= 2).sum() >= 2:
        try:
            sil = float(silhouette_score(Xs, labels, metric="euclidean"))
        except Exception:  # noqa: BLE001
            sil = None

    _scatter(
        pca_coords,
        labels,
        out_path=out_dir / f"{prefix}_pca.png",
        title=f"{title} — PCA",
        xlab=f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}%)",
        ylab=f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}%)",
        focus=focus,
    )
    if umap_coords is not None:
        _scatter(
            umap_coords,
            labels,
            out_path=out_dir / f"{prefix}_umap.png",
            title=f"{title} — UMAP",
            xlab="UMAP1",
            ylab="UMAP2",
            focus=focus,
        )

    # Focus-only panel
    focus_mask = np.array([d in focus for d in disease_ids])
    if focus_mask.sum() >= 5:
        _scatter(
            pca_coords[focus_mask],
            [disease_ids[i] for i, m in enumerate(focus_mask) if m],
            out_path=out_dir / f"{prefix}_pca_focus_pulmonary.png",
            title=f"{title} — COPD / bronchitis / lung cancer (PCA)",
            xlab="PC1",
            ylab="PC2",
            focus=focus,
        )
        if umap_coords is not None:
            _scatter(
                umap_coords[focus_mask],
                [disease_ids[i] for i, m in enumerate(focus_mask) if m],
                out_path=out_dir / f"{prefix}_umap_focus_pulmonary.png",
                title=f"{title} — COPD / bronchitis / lung cancer (UMAP)",
                xlab="UMAP1",
                ylab="UMAP2",
                focus=focus,
            )

    return {
        "n": int(len(X)),
        "pca_var_pc1": float(pca.explained_variance_ratio_[0]),
        "pca_var_pc2": float(pca.explained_variance_ratio_[1]),
        "pca_var_top5": [float(x) for x in pca.explained_variance_ratio_[:5]],
        "silhouette": sil,
        "umap_error": umap_error,
        "pca_coords": pca_coords,
        "umap_coords": umap_coords,
        "loadings_pc1_top": _top_loadings(pca, [f"v{i}" for i in range(X.shape[1])], 0, 8),
    }


def _top_loadings(pca: PCA, names: list[str], comp: int, k: int) -> list[dict[str, float]]:
    if comp >= pca.components_.shape[0]:
        return []
    vec = pca.components_[comp]
    idx = np.argsort(np.abs(vec))[::-1][:k]
    return [{"feature_index": int(i), "loading": float(vec[i])} for i in idx]


_COLOR = {
    "copd": "#c45c26",
    "chronic_bronchitis": "#2a9d8f",
    "lung_adenocarcinoma": "#264653",
    "lung_squamous_cell_carcinoma": "#6d597a",
    "other": "#b0b0b0",
}


def _scatter(
    coords: np.ndarray,
    labels: list[str],
    *,
    out_path: Path,
    title: str,
    xlab: str,
    ylab: str,
    focus: set[str],
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 7.2), dpi=140)
    order = list(dict.fromkeys([*{x for x in labels if x in focus}, "other", *labels]))
    for lab in order:
        mask = np.array([l == lab for l in labels])
        if not mask.any():
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=28 if lab != "other" else 12,
            alpha=0.85 if lab != "other" else 0.25,
            c=_COLOR.get(lab, "#4a5568"),
            label=lab,
            edgecolors="none",
        )
    ax.set_title(title)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _focus_alignment(
    X: np.ndarray, disease_ids: list[str], focus: set[str]
) -> dict[str, Any]:
    centroids: dict[str, np.ndarray] = {}
    for did in sorted(focus):
        idx = [i for i, d in enumerate(disease_ids) if d == did]
        if len(idx) >= 2:
            centroids[did] = X[idx].mean(axis=0)

    dist: dict[str, float] = {}
    keys = list(centroids.keys())
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            dist[f"{a}__vs__{b}"] = float(np.linalg.norm(centroids[a] - centroids[b]))

    # Within-disease mean pairwise distance
    within: dict[str, float] = {}
    for did, cen in centroids.items():
        idx = [i for i, d in enumerate(disease_ids) if d == did]
        if len(idx) < 3:
            continue
        D = pairwise_distances(X[idx])
        within[did] = float(D[np.triu_indices_from(D, k=1)].mean())

    return {
        "centroid_l2_distances": dist,
        "within_disease_mean_pairwise_l2": within,
        "n_per_focus": {d: sum(1 for x in disease_ids if x == d) for d in sorted(focus)},
        "interpretation": (
            "Smaller centroid distances ⇒ more similar VOC profiles. "
            "COPD vs chronic_bronchitis should be closer than either vs lung cancer "
            "if the atlas separates obstructive inflammation from malignancy."
        ),
    }


def _smoking_effect(
    patients_df: pd.DataFrame, voc_df: pd.DataFrame, voc_ids: list[str]
) -> dict[str, Any]:
    """Mean |Δlog2fc| between current vs never smokers within COPD/lung cancer."""
    merged = patients_df.merge(voc_df, on="patient_id", suffixes=("", "_y"))
    out: dict[str, Any] = {}
    for did in ("copd", "lung_adenocarcinoma", "chronic_bronchitis"):
        sub = merged[merged["resolved_disease_id"] == did]
        if sub.empty and "disease_id" in merged.columns:
            sub = merged[merged["disease_id"] == did]
        never = sub[sub["smoking_status"] == "never"]
        current = sub[sub["smoking_status"] == "current"]
        if len(never) < 3 or len(current) < 3:
            out[did] = {"n_never": len(never), "n_current": len(current), "mean_abs_delta": None}
            continue
        cols = [f"log2fc__{v}" for v in voc_ids]
        delta = current[cols].mean().to_numpy() - never[cols].mean().to_numpy()
        top = sorted(
            zip(voc_ids, delta),
            key=lambda kv: abs(float(kv[1])),
            reverse=True,
        )[:8]
        out[did] = {
            "n_never": int(len(never)),
            "n_current": int(len(current)),
            "mean_abs_delta": float(np.mean(np.abs(delta))),
            "top_voc_deltas": [{"voc_id": v, "delta_log2fc": float(d)} for v, d in top],
        }
    return out


def _summarize(
    patients_df: pd.DataFrame,
    break_df: pd.DataFrame,
    resources: dict[str, Any],
    embed_summary: dict[str, Any],
    alignment: dict[str, Any],
    smoking_effect: dict[str, Any],
    *,
    n: int,
) -> dict[str, Any]:
    ok = int(patients_df["ok"].sum()) if "ok" in patients_df.columns else 0
    flag_counts: dict[str, int] = {}
    if not break_df.empty and "flag" in break_df.columns:
        flag_counts = {str(k): int(v) for k, v in break_df["flag"].value_counts().items()}
    return {
        "n_patients": n,
        "n_ok": ok,
        "n_failed": n - ok,
        "ok_pct": round(100.0 * ok / n, 2) if n else None,
        "n_unresolved": int(patients_df.get("unresolved", pd.Series(dtype=bool)).fillna(False).sum())
        if "unresolved" in patients_df.columns
        else 0,
        "breakdown_flag_counts": flag_counts,
        "n_resources_datasources": len(resources.get("datasources") or []),
        "voc_silhouette": (embed_summary.get("voc") or {}).get("silhouette"),
        "gene_shift_silhouette": (embed_summary.get("gene_shift") or {}).get("silhouette"),
        "focus_centroid_distances": alignment.get("centroid_l2_distances"),
        "smoking_effect_mean_abs": {
            k: (v or {}).get("mean_abs_delta") for k, v in smoking_effect.items()
        },
    }


def _md(report: dict[str, Any], break_df: pd.DataFrame, patients_df: pd.DataFrame) -> str:
    o = report["overall"]
    lines = [
        "# 1000-patient diversity cohort evaluation",
        "",
        f"**OK rate: {o['ok_pct']}%** ({o['n_ok']}/{o['n_patients']})",
        "",
        "## Resources engaged",
        "",
    ]
    res = report.get("resources") or {}
    lines.append(
        f"- diseases={res.get('n_diseases')} VOCs={res.get('n_vocs')} pathways={res.get('n_pathways')} "
        f"tissues={res.get('n_tissues')} cell_states={res.get('n_cell_states')} chains={res.get('n_pathway_chains')}"
    )
    lines.append(f"- datasources: {', '.join(res.get('datasources') or [])}")
    lines.append(f"- chembl_loaded={res.get('chembl_loaded')} · OT offline={res.get('opentargets_offline_priors')} · GDC={res.get('gdc_offline_priors')}")
    lines += ["", "## Breakdown flags", ""]
    for k, v in sorted((o.get("breakdown_flag_counts") or {}).items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{k}`: {v}")
    if not (o.get("breakdown_flag_counts") or {}):
        lines.append("- none")
    lines += ["", "## Focus triad VOC alignment (centroid L2)", ""]
    for k, v in sorted(
        ((report.get("alignment") or {}).get("centroid_l2_distances") or {}).items()
    ):
        lines.append(f"- `{k}`: {v:.3f}")
    lines += ["", "## Smoking effect (|mean Δlog2fc| current vs never)", ""]
    for did, row in (report.get("smoking_effect") or {}).items():
        lines.append(
            f"- **{did}**: n_never={row.get('n_never')} n_current={row.get('n_current')} "
            f"mean_abs_delta={row.get('mean_abs_delta')}"
        )
    lines += ["", "## Embeddings", ""]
    for name, emb in (report.get("embed") or {}).items():
        lines.append(
            f"- **{name}**: n={emb.get('n')} PC1={emb.get('pca_var_pc1')} "
            f"PC2={emb.get('pca_var_pc2')} silhouette={emb.get('silhouette')} "
            f"umap_error={emb.get('umap_error')}"
        )
    lines += ["", "## Outputs", ""]
    lines += [
        "- `patients.csv` — covariates + summary metrics",
        "- `patients_dense.csv` — covariates + all VOC log2fc + pathway scores",
        "- `patient_voc_matrix.csv` / `patient_voc_long.csv`",
        "- `patient_gene_shift_matrix.csv` — genetic Δlog2fc vectors",
        "- `patient_pathway_matrix.csv`",
        "- `breakdown_flags.csv`",
        "- `embeddings_voc.csv` / `embeddings_gene_shift.csv`",
        "- `figures/voc_pca.png`, `voc_umap.png`, `gene_shift_*.png`, `*_focus_pulmonary.png`",
        "",
    ]
    return "\n".join(lines)
