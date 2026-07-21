from __future__ import annotations

import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import DATA_DIR, KNOWLEDGE_DIR

CENSUS_DIR = DATA_DIR / "census"
CENSUS_VERSION = "2025-11-08"

# Map Census cell_type strings → ExhalePath cell_state_id (substring match, longest wins)
CELLTYPE_TO_STATE = [
    (r"hepatocyte", "hepatocyte_ketogenic"),
    (r"kupffer", "hepatocyte_sulfur"),
    (r"malignant|cancer cell|tumor cell|neoplastic", "tumor_epithelial_warburg"),
    (r"epithelial cell|enterocyte|goblet|pneumocyte|keratinocyte", "oxidative_stress_cell"),
    (r"microglial", "microglia_activated"),
    (r"neuron|neuroblast|glutamatergic|gabaergic|medium spiny", "neuron_stressed"),
    (r"astrocyte|oligodendrocyte|opc", "microglia_activated"),
    (r"macrophage|monocyte|neutrophil|dendritic", "oxidative_stress_cell"),
    (r"t cell|b cell|plasma cell|nk cell", "oxidative_stress_cell"),
    (r"adipocyte|fat cell", "adipocyte_lipolytic"),
    (r"enterocyte|paneth|tuft|intestinal", "gut_fermentative_microbe"),
    (r"fibroblast|stromal|pericyte|endothelial", "oxidative_stress_cell"),
]


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip().lower())


def _load_top100(path: Path | None = None) -> list[dict[str, Any]]:
    p = Path(path or KNOWLEDGE_DIR / "top100_us_diseases.json")
    return json.loads(p.read_text())["diseases"]


def match_census_diseases(
    census_disease_counts: pd.Series,
    top100: list[dict[str, Any]],
) -> pd.DataFrame:
    """Match top-100 US disease intents to Census disease labels by alias coverage."""
    rows = []
    census_names = list(census_disease_counts.index.astype(str))
    used = set()
    for d in top100:
        aliases = [_norm(a) for a in d.get("aliases", [])] + [_norm(d["name"]), _norm(d["id"])]
        hits = []
        for name in census_names:
            nn = _norm(name)
            if nn == "normal" and d["id"] != "normal_baseline":
                continue
            for a in aliases:
                if not a:
                    continue
                if nn == a or a in nn or nn in a:
                    hits.append(name)
                    break
        # Dedup preserve order
        seen = set()
        uniq = []
        for h in hits:
            if h not in seen:
                seen.add(h)
                uniq.append(h)
        n_cells = int(sum(int(census_disease_counts.get(h, 0)) for h in uniq))
        rows.append(
            {
                "us_disease_id": d["id"],
                "us_disease_name": d["name"],
                "census_diseases": "|".join(uniq),
                "n_census_labels": len(uniq),
                "n_cells": n_cells,
                "matched": len(uniq) > 0,
            }
        )
        used.update(uniq)

    df = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    # Also keep unmatched high-count census diseases not in top100 (data-rich opportunistic harvest)
    leftover = []
    for name, n in census_disease_counts.items():
        if name in used or _norm(str(name)) == "normal":
            continue
        if int(n) < 50_000:
            continue
        leftover.append(
            {
                "us_disease_id": f"census::{_norm(name).replace(' ', '_')[:60]}",
                "us_disease_name": str(name),
                "census_diseases": str(name),
                "n_census_labels": 1,
                "n_cells": int(n),
                "matched": False,
                "opportunistic": True,
            }
        )
    if leftover:
        df = pd.concat([df, pd.DataFrame(leftover)], ignore_index=True)
        df = df.sort_values("n_cells", ascending=False)
    return df


def map_celltype_to_state(cell_type: str) -> str | None:
    ct = _norm(cell_type)
    best = None
    best_len = -1
    for pat, state in CELLTYPE_TO_STATE:
        if re.search(pat, ct):
            if len(pat) > best_len:
                best = state
                best_len = len(pat)
    return best


def harvest_census_compositions(
    *,
    out_dir: Path | None = None,
    census_version: str = CENSUS_VERSION,
    min_disease_cells: int = 5_000,
    top_n_diseases: int = 100,
    include_opportunistic: bool = True,
) -> dict[str, Path]:
    """
    Download CELLxGENE Census cell metadata summaries for:
      - healthy/normal cells across all tissue_general
      - top US diseases (re-ranked by Census cell abundance)
      - opportunistic high-N Census diseases

    Writes compact composition tables (not raw count matrices).
    """
    import cellxgene_census

    out_dir = Path(out_dir or CENSUS_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[census] opening SOMA ({census_version}) …")
    with cellxgene_census.open_soma(census_version=census_version) as census:
        print("[census] loading primary human obs (disease, tissue, cell_type) …")
        obs = cellxgene_census.get_obs(
            census,
            "homo_sapiens",
            value_filter="is_primary_data == True",
            column_names=["disease", "tissue_general", "cell_type", "assay"],
        )

    print(f"[census] loaded {len(obs):,} cells; aggregating …")
    # Keep categoricals — never cast 96M-row columns to wide unicode strings
    disease_counts = obs["disease"].value_counts()
    tissue_counts = obs["tissue_general"].value_counts()

    top100 = _load_top100()
    match_df = match_census_diseases(disease_counts, top100)
    match_path = out_dir / "top100_disease_census_match.csv"
    match_df.to_csv(match_path, index=False)

    # Healthy / normal baseline across all tissues
    disease_cats = set(map(str, disease_counts.index.tolist()))
    normal_label = "normal" if "normal" in disease_cats else next(
        (c for c in disease_cats if str(c).lower() == "normal"), None
    )
    if normal_label is None:
        raise RuntimeError("No 'normal' disease label found in Census obs")
    normal = obs[obs["disease"] == normal_label]
    healthy_comp = (
        normal.groupby(["tissue_general", "cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    healthy_comp["fraction_in_tissue"] = healthy_comp.groupby(
        "tissue_general", observed=True
    )["n_cells"].transform(lambda s: s / s.sum())
    healthy_path = out_dir / "healthy_tissue_celltype_composition.csv"
    healthy_comp.to_csv(healthy_path, index=False)

    tissue_summary = (
        normal.groupby("tissue_general", observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    tissue_summary_path = out_dir / "healthy_tissue_summary.csv"
    tissue_summary.to_csv(tissue_summary_path, index=False)

    # Select diseases to harvest: matched top100 with cells + opportunistic rich diseases
    selected = match_df[match_df["n_cells"] >= min_disease_cells].copy()
    if not include_opportunistic:
        selected = selected[selected.get("matched", True) == True]  # noqa: E712
    # Prefer matched top100 first, then fill with opportunistic by n_cells
    matched = selected[selected["matched"] == True].head(top_n_diseases)  # noqa: E712
    if include_opportunistic and len(matched) < top_n_diseases:
        opp = selected[selected["matched"] == False]  # noqa: E712
        need = top_n_diseases - len(matched)
        selected = pd.concat([matched, opp.head(need)], ignore_index=True)
    else:
        selected = matched

    # Build disease×tissue×cell_type for selected census disease labels
    census_labels = set()
    for labels in selected["census_diseases"]:
        if not labels:
            continue
        census_labels.update([x for x in str(labels).split("|") if x])

    disease_obs = obs[obs["disease"].isin(census_labels)]
    disease_comp = (
        disease_obs.groupby(["disease", "tissue_general", "cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    disease_comp["fraction_in_disease_tissue"] = disease_comp.groupby(
        ["disease", "tissue_general"], observed=True
    )["n_cells"].transform(lambda s: s / max(s.sum(), 1))

    disease_path = out_dir / "disease_tissue_celltype_composition.csv"
    disease_comp.to_csv(disease_path, index=False)

    # Map to ExhalePath states
    def add_state(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["exhalepath_state"] = out["cell_type"].map(map_celltype_to_state)
        return out

    healthy_state = add_state(healthy_comp)
    disease_state = add_state(disease_comp)

    healthy_state_path = out_dir / "healthy_state_composition.csv"
    disease_state_path = out_dir / "disease_state_composition.csv"
    healthy_state.to_csv(healthy_state_path, index=False)
    disease_state.to_csv(disease_state_path, index=False)

    # Per-disease state fractions (pooled across tissues, weighted by n_cells)
    state_rows = []
    for _, sel in selected.iterrows():
        labels = [x for x in str(sel["census_diseases"]).split("|") if x]
        sub = disease_state[disease_state["disease"].isin(labels)]
        if sub.empty:
            continue
        mapped = sub.dropna(subset=["exhalepath_state"])
        if mapped.empty:
            continue
        g = mapped.groupby("exhalepath_state")["n_cells"].sum()
        total = float(g.sum()) or 1.0
        for state_id, n in g.items():
            state_rows.append(
                {
                    "us_disease_id": sel["us_disease_id"],
                    "us_disease_name": sel["us_disease_name"],
                    "exhalepath_state": state_id,
                    "n_cells": int(n),
                    "fraction": float(n) / total,
                    "n_cells_disease_total": int(sel["n_cells"]),
                }
            )
    state_frac = pd.DataFrame(state_rows)
    state_frac_path = out_dir / "disease_exhalepath_state_fractions.csv"
    state_frac.to_csv(state_frac_path, index=False)

    # Healthy state fractions by tissue + global
    hs = healthy_state.dropna(subset=["exhalepath_state"])
    healthy_global = hs.groupby("exhalepath_state")["n_cells"].sum()
    hg_total = float(healthy_global.sum()) or 1.0
    healthy_global_df = healthy_global.reset_index(name="n_cells")
    healthy_global_df["fraction"] = healthy_global_df["n_cells"] / hg_total
    healthy_global_path = out_dir / "healthy_exhalepath_state_fractions.csv"
    healthy_global_df.to_csv(healthy_global_path, index=False)

    # Assay diversity (data-quality signal)
    assay_div = (
        obs.groupby(["disease", "assay"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    assay_path = out_dir / "disease_assay_diversity.csv"
    assay_div.to_csv(assay_path, index=False)

    manifest = {
        "census_version": census_version,
        "n_primary_cells": int(len(obs)),
        "n_normal_cells": int(len(normal)),
        "n_diseases_in_census": int(disease_counts.shape[0]),
        "n_tissues": int(tissue_counts.shape[0]),
        "n_cell_types": int(obs["cell_type"].nunique()),
        "n_top100_matched_with_cells": int(
            ((match_df["matched"] == True) & (match_df["n_cells"] > 0)).sum()  # noqa: E712
        ),
        "n_diseases_harvested": int(len(selected)),
        "min_disease_cells": min_disease_cells,
        "top_harvested": selected.head(30)[
            ["us_disease_id", "us_disease_name", "n_cells"]
        ].to_dict(orient="records"),
        "largest_healthy_tissues": tissue_summary.head(25).to_dict(orient="records"),
        "note": (
            "Summaries are cell-type compositions from CELLxGENE Census primary human cells. "
            "Raw count matrices are not downloaded; fractions calibrate ExhalePath cell-state densities."
        ),
    }
    manifest_path = out_dir / "census_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2)[:2000])

    # Free memory hint
    del obs, normal, disease_obs

    return {
        "manifest": manifest_path,
        "match": match_path,
        "healthy_composition": healthy_path,
        "healthy_tissue_summary": tissue_summary_path,
        "disease_composition": disease_path,
        "healthy_state": healthy_state_path,
        "disease_state": disease_state_path,
        "disease_state_fractions": state_frac_path,
        "healthy_state_fractions": healthy_global_path,
        "assay_diversity": assay_path,
    }


def calibrate_atlas_from_census(
    *,
    census_dir: Path | None = None,
    atlas_path: Path | None = None,
    out_path: Path | None = None,
) -> Path:
    """
    Update cell_state_atlas.json baseline densities from healthy Census fractions,
    and disease_density_mult from disease/healthy fraction ratios (data-rich diseases).
    """
    census_dir = Path(census_dir or CENSUS_DIR)
    atlas_path = Path(atlas_path or KNOWLEDGE_DIR / "cell_state_atlas.json")
    out_path = Path(out_path or KNOWLEDGE_DIR / "cell_state_atlas.census_calibrated.json")

    atlas = json.loads(atlas_path.read_text())
    healthy = pd.read_csv(census_dir / "healthy_exhalepath_state_fractions.csv")
    disease = pd.read_csv(census_dir / "disease_exhalepath_state_fractions.csv")
    match = pd.read_csv(census_dir / "top100_disease_census_match.csv")

    # Map us_disease_id → exhalepath disease_id where possible
    us_to_exhale = {
        "cancer_lung": "lung_adenocarcinoma",
        "cancer_breast": "breast_invasive_carcinoma",
        "cancer_colorectal": "colon_adenocarcinoma",
        "cancer_pancreas": "pancreatic_adenocarcinoma",
        "cancer_liver": "hepatocellular_carcinoma",
        "cancer_ovary": "ovarian_cancer",
        "glioblastoma": "glioblastoma",
        "alzheimer": "alzheimer_disease",
        "parkinson": "parkinson_disease",
        "schizophrenia": "schizophrenia",
        "depression": "major_depressive_disorder",
        "epilepsy": "epilepsy",
        "multiple_sclerosis": "multiple_sclerosis",
        "type2_diabetes": "type_2_diabetes",
        "ckd": "chronic_kidney_disease",
        "cirrhosis": "chronic_liver_disease",
        "ibd": "inflammatory_bowel_disease",
        "autism": "autism_spectrum_disorder",
        "hpylori": "helicobacter_pylori_infection",
        "cdiff": "clostridioides_difficile_infection",
        "lewy_body": "alzheimer_disease",
    }
    # Prefer disease rows with the most cells when duplicates map to same exhalepath id
    disease = disease.sort_values("n_cells_disease_total", ascending=False)

    healthy_frac = {
        r.exhalepath_state: float(r.fraction) for r in healthy.itertuples()
    }

    # Build disease fraction matrix
    dis_frac: dict[str, dict[str, float]] = defaultdict(dict)
    for r in disease.itertuples():
        eid = us_to_exhale.get(r.us_disease_id)
        if not eid:
            continue
        # Prefer larger n_cells diseases — already one row per state
        dis_frac[eid][r.exhalepath_state] = float(r.fraction)

    for st in atlas["cell_states"]:
        sid = st["state_id"]
        if sid in healthy_frac:
            # Blend atlas prior with census healthy fraction
            st["baseline_density"] = float(
                np.clip(0.35 * float(st.get("baseline_density", 0.1)) + 0.65 * healthy_frac[sid], 0.01, 0.9)
            )
            st["census_healthy_fraction"] = healthy_frac[sid]

        dens_mult = dict(st.get("disease_density_mult") or {})
        for eid, fracs in dis_frac.items():
            if sid not in fracs:
                continue
            h = max(healthy_frac.get(sid, 1e-3), 1e-3)
            ratio = float(np.clip(fracs[sid] / h, 0.25, 8.0))
            # Blend with existing expert prior if present
            prior = float(dens_mult.get(eid, 1.0))
            dens_mult[eid] = float(np.clip(0.4 * prior + 0.6 * ratio, 0.3, 8.0))
        st["disease_density_mult"] = dens_mult
        st["census_calibrated"] = True

    atlas["version"] = str(atlas.get("version", "1.0.0")) + "+census"
    atlas["census_calibration"] = {
        "source": str(census_dir / "census_manifest.json"),
        "n_states_with_healthy": len(healthy_frac),
        "n_diseases_mapped": len(dis_frac),
    }
    out_path.write_text(json.dumps(atlas, indent=2))

    # Also overwrite active atlas used by the tool
    active = KNOWLEDGE_DIR / "cell_state_atlas.json"
    active.write_text(json.dumps(atlas, indent=2))
    print(f"[census] calibrated atlas → {active}")
    return active
