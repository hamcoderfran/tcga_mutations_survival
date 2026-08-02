"""Cutting-edge adaptive fusion: anti-dilution, mode-aware weights, epistemic UQ."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from .data import load_json
from .types import AspectHit, FusedVOC, ModelOutput, StackQuery

# Family multipliers by prediction regime
_ATLAS_FAMILY_MULT = {
    "exhalepath": 1.35,
    "ml": 1.25,
    "literature": 0.85,
    "qsar": 1.05,
    "pbpk": 1.05,
    "zero_shot": 0.75,
    "gsmm": 0.95,
    "graph": 0.9,
    "genetics": 0.95,
    "enrichment": 0.95,
    "microbiome": 0.85,
    "single_cell": 1.0,
    "clinical": 0.9,
    "pharmacology": 0.7,
    "network": 0.85,
    "meta": 1.15,
    "null": 0.4,
    "phenotype": 0.9,
}

_ZS_FAMILY_MULT = {
    "exhalepath": 0.85,
    "ml": 0.7,
    "literature": 0.25,  # avoid circular VOC prior copy
    "qsar": 1.1,
    "pbpk": 1.1,
    "zero_shot": 1.55,
    "gsmm": 1.4,
    "graph": 1.35,
    "genetics": 1.3,
    "enrichment": 1.35,
    "microbiome": 1.05,
    "single_cell": 1.1,
    "clinical": 0.8,
    "pharmacology": 0.75,
    "network": 1.2,
    "meta": 1.1,
    "null": 0.5,
    "phenotype": 1.35,
}

_ANCHOR_IDS = {
    "exhalepath_hybrid",
    "exhalepath_physiology",
    "exhalepath_calibrator",
}


def detect_zero_shot(model_outputs: list[ModelOutput], query: StackQuery) -> bool:
    """True only for genuinely unresolved / custom diseases.

    Do not treat every zero_shot_evidence theme hit as zero-shot — that would
    flip atlas diseases into the ZS weight regime and dilute hybrid anchors.
    """
    for mo in model_outputs:
        if mo.model_id != "zero_shot_mechanism":
            continue
        if mo.metadata.get("zero_shot"):
            return True
        did = str(mo.metadata.get("disease_id") or "")
        if did.startswith("custom"):
            return True
    if str(query.disease).lower().startswith("custom"):
        return True
    return False


def calibrated_weight_overrides() -> dict[str, float]:
    doc = load_json("fusion_weights_calibrated.json")
    return {str(k): float(v) for k, v in (doc.get("weights") or {}).items()}


def effective_weight(
    mo: ModelOutput,
    *,
    zero_shot: bool,
    overrides: dict[str, float],
) -> float:
    if mo.status == "skipped" or mo.weight <= 0:
        return 0.0
    base = float(overrides.get(mo.model_id, mo.weight))
    fam = _ZS_FAMILY_MULT if zero_shot else _ATLAS_FAMILY_MULT
    return base * float(fam.get(mo.family, 1.0))


def reciprocal_rank_fusion(
    model_outputs: list[ModelOutput],
    *,
    weights: dict[str, float],
    k: int = 60,
) -> dict[str, float]:
    scores: dict[str, float] = defaultdict(float)
    for mo in model_outputs:
        w = weights.get(mo.model_id, 0.0)
        if w <= 0 or mo.status == "skipped":
            continue
        ranked = sorted(mo.voc_signals, key=lambda s: abs(s.log2fc), reverse=True)
        for r, sig in enumerate(ranked, 1):
            scores[sig.voc_id] += w * (1.0 / (k + r))
    return dict(scores)


def fuse_vocs(
    model_outputs: list[ModelOutput],
    *,
    voc_catalog: dict[str, Any],
    top_n: int = 20,
    query: StackQuery | None = None,
    zero_shot: bool | None = None,
) -> list[FusedVOC]:
    q = query or StackQuery(disease="unknown")
    zs = detect_zero_shot(model_outputs, q) if zero_shot is None else zero_shot
    overrides = calibrated_weight_overrides()
    eff = {mo.model_id: effective_weight(mo, zero_shot=zs, overrides=overrides) for mo in model_outputs}
    rrf = reciprocal_rank_fusion(model_outputs, weights=eff)

    by_voc: dict[str, list[tuple[ModelOutput, Any]]] = defaultdict(list)
    for mo in model_outputs:
        if eff.get(mo.model_id, 0.0) <= 0:
            continue
        for sig in mo.voc_signals:
            by_voc[sig.voc_id].append((mo, sig))

    fused: list[FusedVOC] = []
    for vid, pairs in by_voc.items():
        w_sum = 0.0
        num = 0.0
        conf_num = 0.0
        conf_den = 0.0
        delta_num = 0.0
        delta_den = 0.0
        votes: dict[str, float] = {}
        confs: dict[str, float] = {}
        evidence: list[str] = []
        agreeing = 0
        weighted_vals: list[tuple[float, float]] = []

        for mo, sig in pairs:
            w = eff[mo.model_id] * max(0.05, float(sig.confidence))
            fc = float(sig.log2fc)
            num += w * fc
            w_sum += w
            conf_num += eff[mo.model_id] * float(sig.confidence)
            conf_den += eff[mo.model_id]
            votes[mo.model_id] = fc
            confs[mo.model_id] = float(sig.confidence)
            weighted_vals.append((fc, w))
            if sig.delta_ppb is not None:
                delta_num += w * float(sig.delta_ppb)
                delta_den += w
            if abs(fc) >= 0.05:
                agreeing += 1
            for e in sig.evidence[:2]:
                evidence.append(f"[{mo.model_id}] {e}")
        if w_sum <= 0:
            continue

        log2fc = num / w_sum

        # --- anti-dilution anchor: hybrid/physio/calibrator ---
        # Prefer hybrid+calibrator when they agree (physiology can disagree on
        # sparse VOCs). Fall back to full 3-anchor agreement.
        def _agreeing(ids: set[str]) -> list[float]:
            vals = [
                votes[mid]
                for mid in ids
                if mid in votes and abs(votes[mid]) >= 0.03
            ]
            if len(vals) < 2:
                return []
            signs = [1 if v >= 0 else -1 for v in vals]
            return vals if abs(sum(signs)) == len(signs) else []

        primary = _agreeing({"exhalepath_hybrid", "exhalepath_calibrator"})
        fallback = _agreeing(_ANCHOR_IDS)
        anchor_vals = primary or fallback
        if anchor_vals:
            anchor_mean = sum(anchor_vals) / len(anchor_vals)
            # atlas: trust anchors more; zero-shot: softer pull
            # hybrid+calibrator agreement gets a stronger atlas pull
            if primary:
                alpha = 0.5 if zs else 0.82
                tag = "hybrid_calibrator"
            else:
                alpha = 0.45 if zs else 0.72
                tag = "full_anchor"
            log2fc = (1.0 - alpha) * log2fc + alpha * anchor_mean
            evidence.append(f"[fusion] anti_dilution_{tag} alpha={alpha:.2f}")

        # epistemic uncertainty = weighted std of votes
        if sum(w for _, w in weighted_vals) > 0:
            mean = sum(v * w for v, w in weighted_vals) / sum(w for _, w in weighted_vals)
            var = sum(w * (v - mean) ** 2 for v, w in weighted_vals) / sum(
                w for _, w in weighted_vals
            )
            epi = math.sqrt(max(var, 0.0))
        else:
            epi = 0.0

        base_conf = conf_num / max(conf_den, 1e-9)
        agree_bonus = min(0.25, 0.03 * agreeing)
        signs = [1 if v >= 0 else -1 for v in votes.values() if abs(v) >= 0.05]
        sign_agree = abs(sum(signs)) / len(signs) if signs else 0.5
        # shrink confidence when epistemic uncertainty high
        uq_penalty = min(0.35, 0.18 * epi)
        fused_conf = max(
            0.05,
            min(0.99, base_conf * (0.55 + 0.45 * sign_agree) + agree_bonus - uq_penalty),
        )

        name = (voc_catalog.get(vid) or {}).get("name") or vid
        healthy = float((voc_catalog.get(vid) or {}).get("healthy_ppb") or 10.0)
        if delta_den > 0:
            delta = delta_num / delta_den
        else:
            delta = healthy * (2.0**log2fc - 1.0)

        # 90% CI from epistemic std (approximate)
        z = 1.645
        ci_low = log2fc - z * epi
        ci_high = log2fc + z * epi

        fused.append(
            FusedVOC(
                voc_id=vid,
                name=name,
                fused_log2fc=float(log2fc),
                fused_delta_ppb=float(delta),
                fused_confidence=float(fused_conf),
                rrf_score=float(rrf.get(vid, 0.0)),
                n_models_agreeing=agreeing,
                model_votes=votes,
                model_confidences=confs,
                evidence=evidence[:14],
                epistemic_std=float(epi),
                ci_low_log2fc=float(ci_low),
                ci_high_log2fc=float(ci_high),
            )
        )

    def sort_key(f: FusedVOC) -> float:
        # prefer high |effect| * conf / (1+uncertainty), plus RRF
        return (
            0.55 * abs(f.fused_log2fc) * f.fused_confidence / (1.0 + f.epistemic_std)
            + 0.30 * f.rrf_score
            + 0.15 * (f.n_models_agreeing / max(1, len(f.model_votes)))
        )

    fused.sort(key=sort_key, reverse=True)
    return fused[:top_n]


def fuse_aspects(model_outputs: list[ModelOutput]) -> dict[str, list[AspectHit]]:
    buckets: dict[str, dict[str, AspectHit]] = defaultdict(dict)
    for mo in model_outputs:
        if mo.status == "skipped":
            continue
        for hit in mo.aspects:
            kind = hit.kind or mo.aspect
            prev = buckets[kind].get(hit.id)
            score = hit.score * mo.weight
            if prev is None:
                buckets[kind][hit.id] = AspectHit(
                    id=hit.id,
                    name=hit.name,
                    score=score,
                    kind=kind,
                    evidence=list(hit.evidence) + [f"from:{mo.model_id}"],
                )
            else:
                prev.score += score
                prev.evidence.extend(hit.evidence[:2])
                prev.evidence.append(f"from:{mo.model_id}")
    out: dict[str, list[AspectHit]] = {}
    for kind, items in buckets.items():
        ranked = sorted(items.values(), key=lambda h: abs(h.score), reverse=True)
        out[kind] = ranked[:15]
    return out


def fusion_summary(
    query: StackQuery,
    model_outputs: list[ModelOutput],
    fused_vocs: list[FusedVOC],
    *,
    zero_shot: bool = False,
) -> dict[str, Any]:
    ok = [m for m in model_outputs if m.status == "ok"]
    degraded = [m for m in model_outputs if m.status == "degraded"]
    skipped = [m for m in model_outputs if m.status == "skipped"]
    mean_agree = (
        sum(v.n_models_agreeing for v in fused_vocs) / len(fused_vocs) if fused_vocs else 0.0
    )
    mean_epi = (
        sum(v.epistemic_std for v in fused_vocs) / len(fused_vocs) if fused_vocs else 0.0
    )
    return {
        "disease": query.disease,
        "zero_shot_mode": zero_shot,
        "n_models_total": len(model_outputs),
        "n_models_ok": len(ok),
        "n_models_degraded": len(degraded),
        "n_models_skipped": len(skipped),
        "models_ok": [m.model_id for m in ok],
        "models_degraded": [m.model_id for m in degraded],
        "models_skipped": [m.model_id for m in skipped],
        "n_fused_vocs": len(fused_vocs),
        "mean_model_agreement": round(mean_agree, 3),
        "mean_epistemic_std": round(mean_epi, 4),
        "top_voc": fused_vocs[0].voc_id if fused_vocs else None,
        "fusion_method": (
            "adaptive_mode_aware + anti_dilution_anchor + RRF + "
            "sign_agreement + epistemic_uq"
            + (" + calibrated_weights" if calibrated_weight_overrides() else "")
        ),
    }
