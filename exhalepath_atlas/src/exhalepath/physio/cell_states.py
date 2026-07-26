from __future__ import annotations

from typing import Any, Iterable

from ..schemas import CellStateActivity, TumorContext
from ..utils import stage_ordinal


def _marker_hit_score(markers: list[str], mutated: set[str], associated: dict[str, float]) -> float:
    if not markers:
        return 0.0
    hits = 0.0
    for m in markers:
        mu = m.upper()
        if mu in mutated:
            hits += 1.0
        elif mu in associated:
            hits += min(float(associated[mu]), 1.0)
    return min(hits / len(markers), 1.0)


def estimate_cell_states(
    *,
    atlas: list[dict[str, Any]],
    disease: dict[str, Any],
    mutated_genes: Iterable[str],
    associated_genes: dict[str, float] | None = None,
    tumor: TumorContext | None = None,
    cell_state_fractions: dict[str, float] | None = None,
    cell_state_activity: dict[str, float] | None = None,
) -> list[CellStateActivity]:
    """
    Predict affected cell-state density and activity.

    Density ≈ baseline × disease_mult × (1 + marker_expansion) × tumor_burden
    Activity ≈ baseline × disease_mult × (1 + marker_activation)
    Optional user overrides from single-cell fractions / activity.
    """
    mutated = {g.upper() for g in mutated_genes if g}
    associated_genes = {str(k).upper(): float(v) for k, v in (associated_genes or {}).items()}
    fractions = cell_state_fractions or {}
    activities = cell_state_activity or {}
    did = disease.get("disease_id", "")

    burden = 1.0
    if tumor and tumor.tumor_burden_proxy is not None:
        burden = 0.7 + 0.8 * tumor.tumor_burden_proxy
    elif tumor and tumor.stage:
        # Map stage ordinal 0–4 → ~0.85–1.45
        burden = 0.85 + 0.15 * stage_ordinal(tumor.stage)
    if tumor and tumor.metastatic:
        burden *= 1.2

    out: list[CellStateActivity] = []
    for st in atlas:
        sid = st["state_id"]
        markers = list(st.get("markers") or [])
        marker_score = _marker_hit_score(markers, mutated, associated_genes)

        dens_mult = float(st.get("disease_density_mult", {}).get(did, 1.0))
        act_mult = float(st.get("disease_activity_mult", {}).get(did, 1.0))
        # Soft-max with comorbidity / zero-shot cell-state donor ids
        donor_ids = list(disease.get("_comorbid_ids") or [])
        donor_ids.extend(list(disease.get("_cell_state_donor_ids") or []))
        for cid in donor_ids:
            dens_mult = max(
                dens_mult, float(st.get("disease_density_mult", {}).get(cid, 1.0))
            )
            act_mult = max(
                act_mult, float(st.get("disease_activity_mult", {}).get(cid, 1.0))
            )
        disease_modulated = dens_mult != 1.0 or act_mult != 1.0 or marker_score >= 0.08

        density = float(st.get("baseline_density", 0.05)) * dens_mult * (1.0 + 0.8 * marker_score)
        activity = float(st.get("baseline_activity", 0.2)) * act_mult * (1.0 + 0.6 * marker_score)

        # Tumor epithelial states scale with burden; brain/gut less so
        tissue = (st.get("tissue") or "").lower()
        if tissue in {"tumor", "lung", "multi"} or did.endswith("carcinoma") or "cancer" in (disease.get("category") or ""):
            if sid.startswith("tumor_") or tissue == "tumor":
                density *= burden

        if sid in fractions:
            density = max(float(fractions[sid]), 0.0)
        if sid in activities:
            activity = max(float(activities[sid]), 0.0)

        # Effective source strength
        effective = density * activity
        out.append(
            CellStateActivity(
                state_id=sid,
                name=st.get("name", sid),
                tissue=st.get("tissue", "unknown"),
                density=float(density),
                activity=float(activity),
                effective_source=float(effective),
                marker_hit_score=float(marker_score),
                marker_genes_hit=sorted({m.upper() for m in markers if m.upper() in mutated or m.upper() in associated_genes}),
                produces_chains=list(st.get("produces_chains") or []),
                disease_modulated=bool(disease_modulated),
            )
        )

    out.sort(key=lambda x: x.effective_source, reverse=True)
    return out
