from __future__ import annotations

from typing import Any

from ..schemas import CellStateActivity, PhysiologyTrace


def farhi_alveolar_fraction(lambda_ba: float, va: float, q: float) -> float:
    """
    Alveolar release efficiency from mixed-venous blood.

    Uses excretion form: f = 1 / (1 + λ·Q/VA)
    Low-λ gases (pentane/isoprene) are exhaled efficiently; high-λ solvents
    (acetone/ethanol) are more blood-retained per pass.
    """
    lam = max(float(lambda_ba), 1e-6)
    return float(1.0 / (1.0 + (lam * q) / max(va, 1e-9)))


def blood_transfer_rate(
    *,
    tissue_production: float,
    tissue: str,
    voc_id: str,
    physio: dict[str, Any],
) -> dict[str, float]:
    """Map tissue VOC production into mixed-venous blood delivery (relative units)."""
    perfusion = float(physio.get("tissue_perfusion_fraction", {}).get(tissue, 0.1))
    first_pass_table = physio.get("hepatic_first_pass", {})
    first_pass = float(first_pass_table.get(voc_id, first_pass_table.get("default", 0.15)))
    if tissue in {"gut", "stomach", "liver"}:
        first_pass = min(0.85, first_pass * 1.15)
    if tissue == "tumor":
        first_pass = min(first_pass, 0.25)
    delivered = max(tissue_production, 0.0) * perfusion * (1.0 - first_pass)
    return {
        "perfusion_fraction": perfusion,
        "hepatic_first_pass": first_pass,
        "blood_delivery": float(delivered),
    }


def alveolar_release_ppb(
    *,
    blood_delivery: float,
    healthy_blood_delivery: float,
    healthy_ppb: float,
    voc_id: str,
    physio: dict[str, Any],
) -> dict[str, float]:
    """Convert relative blood delivery into exhaled alveolar ppb via Farhi exchange."""
    physiology = physio.get("physiology", {})
    va = float(physiology.get("alveolar_ventilation_l_per_min", 5.0))
    q = float(physiology.get("cardiac_output_l_per_min", 5.5))
    lam_table = physio.get("blood_air_partition_lambda", {})
    lam = float(lam_table.get(voc_id, lam_table.get("default", 50.0)))
    alv_frac = farhi_alveolar_fraction(lam, va, q)

    denom = max(healthy_blood_delivery * max(alv_frac, 1e-12), 1e-12)
    scale = healthy_ppb / denom
    predicted = max(blood_delivery, 0.0) * alv_frac * scale
    return {
        "lambda_blood_air": lam,
        "alveolar_fraction": alv_frac,
        "va_l_per_min": va,
        "q_l_per_min": q,
        "predicted_ppb": float(max(predicted, 0.0)),
    }


def integrate_production(
    *,
    voc_id: str,
    chain_fluxes: list[dict[str, Any]],
    cell_states: list[CellStateActivity],
    physio: dict[str, Any],
    healthy_ppb: float,
    disease_prior_log2fc: float = 0.0,
) -> PhysiologyTrace:
    """
    Full stack for one VOC:
      cell density×activity × chain flux → tissue production
      → blood transfer → alveolar release → exhaled ppb
    """
    relevant_chains = [c for c in chain_fluxes if c["voc_id"] == voc_id]
    chain_ids = {c["chain_id"] for c in relevant_chains}
    flux_by_chain = {c["chain_id"]: float(c["flux"]) for c in relevant_chains}

    physiology = physio.get("physiology", {})
    scale = float(physiology.get("production_scale", 6.0))
    healthy_prod = float(physiology.get("healthy_reference_production", 0.08))
    prior_mult = float(2 ** disease_prior_log2fc)

    tissue_prod: dict[str, float] = {}
    contributing_states: list[str] = []
    contributing_chains: list[str] = []

    for st in cell_states:
        for cid in st.produces_chains:
            if cid not in chain_ids:
                continue
            flux = flux_by_chain.get(cid, 0.0)
            # Only disease-modulated states or strong chain flux contribute
            if st.disease_modulated:
                flux_eff = max(flux, 0.12)
            elif flux >= 0.20:
                flux_eff = flux
            elif abs(disease_prior_log2fc) >= 0.25 and flux >= 0.08:
                flux_eff = flux
            else:
                continue
            prod = st.effective_source * flux_eff * scale
            if prod <= 1e-9:
                continue
            tissue_prod[st.tissue] = tissue_prod.get(st.tissue, 0.0) + prod
            contributing_states.append(st.state_id)
            contributing_chains.append(cid)

    if tissue_prod:
        for t in list(tissue_prod):
            tissue_prod[t] *= prior_mult
    else:
        # Healthy-like systemic baseline × disease prior (can be <1 for suppressed VOCs)
        tissue_prod = {"multi": healthy_prod * prior_mult}
        contributing_states = []
        contributing_chains = []

    blood_total = 0.0
    perfusion_used = 0.0
    first_pass_used = 0.0
    for tissue, prod in tissue_prod.items():
        bt = blood_transfer_rate(
            tissue_production=prod, tissue=tissue, voc_id=voc_id, physio=physio
        )
        blood_total += bt["blood_delivery"]
        perfusion_used = max(perfusion_used, bt["perfusion_fraction"])
        first_pass_used = max(first_pass_used, bt["hepatic_first_pass"])

    healthy_bt = blood_transfer_rate(
        tissue_production=healthy_prod, tissue="multi", voc_id=voc_id, physio=physio
    )
    healthy_delivery = max(healthy_bt["blood_delivery"], 1e-9)

    alv = alveolar_release_ppb(
        blood_delivery=blood_total,
        healthy_blood_delivery=healthy_delivery,
        healthy_ppb=healthy_ppb,
        voc_id=voc_id,
        physio=physio,
    )

    return PhysiologyTrace(
        voc_id=voc_id,
        tissue_production=float(sum(tissue_prod.values())),
        blood_delivery=float(blood_total),
        hepatic_first_pass=float(first_pass_used),
        perfusion_fraction=float(perfusion_used),
        lambda_blood_air=float(alv["lambda_blood_air"]),
        alveolar_fraction=float(alv["alveolar_fraction"]),
        alveolar_ventilation_l_per_min=float(alv["va_l_per_min"]),
        cardiac_output_l_per_min=float(alv["q_l_per_min"]),
        predicted_ppb=float(alv["predicted_ppb"]),
        contributing_cell_states=sorted(set(contributing_states)),
        contributing_chains=sorted(set(contributing_chains)),
        tissue_breakdown={k: float(v) for k, v in tissue_prod.items()},
    )
