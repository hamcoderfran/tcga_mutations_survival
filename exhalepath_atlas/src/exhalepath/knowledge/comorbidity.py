"""Merge primary disease priors with comorbidities for whole-body VOC prediction."""

from __future__ import annotations

from copy import deepcopy
from typing import Any


DEFAULT_COMORBIDITY_WEIGHT = 0.65


def merge_disease_with_comorbidities(
    primary: dict[str, Any],
    comorbid_diseases: list[dict[str, Any]],
    *,
    weight: float = DEFAULT_COMORBIDITY_WEIGHT,
    weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    """
    Fuse comorbidity pathway biases and VOC log2FC priors into the primary disease.

    - pathway_bias: soft-max style — keep primary, boost toward comorbid biases
      ``1 + w * (c - 1)`` then take element-wise max with primary
    - voc_log2fc_prior: ``primary + w * comorbid`` (clipped later by predictor)
    - metadata: ``_comorbidities``, ``_comorbid_ids`` for cell-state / explain layers
    """
    if not comorbid_diseases:
        out = deepcopy(primary)
        out["_comorbidities"] = []
        out["_comorbid_ids"] = []
        return out

    out = deepcopy(primary)
    pb = dict(out.get("pathway_bias") or {})
    vp = dict(out.get("voc_log2fc_prior") or {})
    meta: list[dict[str, Any]] = []
    ids: list[str] = []

    for c in comorbid_diseases:
        cid = c.get("disease_id") or c.get("name") or "unknown"
        ids.append(str(cid))
        w = float((weights or {}).get(str(cid), weight))
        w = max(0.0, min(1.5, w))
        meta.append(
            {
                "disease_id": c.get("disease_id"),
                "name": c.get("name"),
                "category": c.get("category"),
                "weight": w,
                "default_site": c.get("default_site"),
            }
        )
        for pid, bias in (c.get("pathway_bias") or {}).items():
            b = float(bias)
            boosted = 1.0 + w * (b - 1.0)
            pb[pid] = max(float(pb.get(pid, 1.0)), boosted)
        for vid, prior in (c.get("voc_log2fc_prior") or {}).items():
            vp[vid] = float(vp.get(vid, 0.0)) + w * float(prior)

    out["pathway_bias"] = pb
    out["voc_log2fc_prior"] = vp
    out["_comorbidities"] = meta
    out["_comorbid_ids"] = ids
    # Keep primary identity; annotate display name
    names = [m.get("name") or m.get("disease_id") for m in meta]
    out["_comorbidity_label"] = "+".join(str(n) for n in names if n)
    return out


def resolve_comorbid_diseases(
    kb: Any,
    comorbidity_queries: list[str],
) -> list[dict[str, Any]]:
    """Resolve free-text comorbidity names via KnowledgeBase.resolve_disease."""
    out: list[dict[str, Any]] = []
    seen: set[str] = set()
    for q in comorbidity_queries:
        q = (q or "").strip()
        if not q:
            continue
        d = kb.resolve_disease(q)
        did = d.get("disease_id") or q
        if did in seen or d.get("_unresolved"):
            # still allow unresolved with empty priors? skip unresolved
            if d.get("_unresolved"):
                continue
        if did in seen:
            continue
        seen.add(did)
        out.append(d)
    return out
