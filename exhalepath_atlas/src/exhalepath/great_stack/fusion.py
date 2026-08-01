"""Fuse 10+ model VOC signals + disease aspects into one consensus."""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

from .types import AspectHit, FusedVOC, ModelOutput, StackQuery


def reciprocal_rank_fusion(
    model_outputs: list[ModelOutput], *, k: int = 60
) -> dict[str, float]:
    scores: dict[str, float] = defaultdict(float)
    for mo in model_outputs:
        if mo.status == "skipped" or mo.weight <= 0:
            continue
        ranked = sorted(mo.voc_signals, key=lambda s: abs(s.log2fc), reverse=True)
        for r, sig in enumerate(ranked, 1):
            scores[sig.voc_id] += mo.weight * (1.0 / (k + r))
    return dict(scores)


def fuse_vocs(
    model_outputs: list[ModelOutput],
    *,
    voc_catalog: dict[str, Any],
    top_n: int = 20,
) -> list[FusedVOC]:
    rrf = reciprocal_rank_fusion(model_outputs)
    by_voc: dict[str, list[tuple[ModelOutput, Any]]] = defaultdict(list)
    for mo in model_outputs:
        if mo.status == "skipped" or mo.weight <= 0:
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
        for mo, sig in pairs:
            w = mo.weight * max(0.05, float(sig.confidence))
            num += w * float(sig.log2fc)
            w_sum += w
            conf_num += mo.weight * float(sig.confidence)
            conf_den += mo.weight
            votes[mo.model_id] = float(sig.log2fc)
            confs[mo.model_id] = float(sig.confidence)
            if sig.delta_ppb is not None:
                delta_num += w * float(sig.delta_ppb)
                delta_den += w
            if abs(sig.log2fc) >= 0.05:
                agreeing += 1
            for e in sig.evidence[:2]:
                evidence.append(f"[{mo.model_id}] {e}")
        if w_sum <= 0:
            continue
        log2fc = num / w_sum
        # agreement bonus: more independent models → higher confidence
        base_conf = conf_num / max(conf_den, 1e-9)
        agree_bonus = min(0.25, 0.03 * agreeing)
        # shrink toward 0 when models disagree in sign
        signs = [1 if v >= 0 else -1 for v in votes.values() if abs(v) >= 0.05]
        if signs:
            sign_agree = abs(sum(signs)) / len(signs)
        else:
            sign_agree = 0.5
        fused_conf = max(0.05, min(0.99, base_conf * (0.55 + 0.45 * sign_agree) + agree_bonus))
        name = (voc_catalog.get(vid) or {}).get("name") or vid
        # approximate delta from log2fc if missing (relative to ~10 ppb baseline)
        healthy = float((voc_catalog.get(vid) or {}).get("healthy_ppb") or 10.0)
        if delta_den > 0:
            delta = delta_num / delta_den
        else:
            delta = healthy * (2.0**log2fc - 1.0)
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
                evidence=evidence[:12],
            )
        )

    # final rank: combine |log2fc| * confidence with RRF
    def sort_key(f: FusedVOC) -> float:
        return 0.65 * abs(f.fused_log2fc) * f.fused_confidence + 0.35 * f.rrf_score

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
) -> dict[str, Any]:
    ok = [m for m in model_outputs if m.status == "ok"]
    degraded = [m for m in model_outputs if m.status == "degraded"]
    skipped = [m for m in model_outputs if m.status == "skipped"]
    mean_agree = (
        sum(v.n_models_agreeing for v in fused_vocs) / len(fused_vocs) if fused_vocs else 0.0
    )
    return {
        "disease": query.disease,
        "n_models_total": len(model_outputs),
        "n_models_ok": len(ok),
        "n_models_degraded": len(degraded),
        "n_models_skipped": len(skipped),
        "models_ok": [m.model_id for m in ok],
        "models_degraded": [m.model_id for m in degraded],
        "models_skipped": [m.model_id for m in skipped],
        "n_fused_vocs": len(fused_vocs),
        "mean_model_agreement": round(mean_agree, 3),
        "top_voc": fused_vocs[0].voc_id if fused_vocs else None,
        "fusion_method": "confidence_weighted_log2fc + reciprocal_rank_fusion + sign_agreement",
    }
