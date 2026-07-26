"""Offline disease→gene index beyond atlas disease_ids (zero-shot)."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..config import KNOWLEDGE_DIR


def _norm(s: str) -> str:
    s = str(s or "").strip().lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


@lru_cache(maxsize=1)
def load_zero_shot_gene_index(path: str | None = None) -> dict[str, Any]:
    p = Path(path or KNOWLEDGE_DIR / "zero_shot_disease_genes.json")
    if not p.exists():
        return {"diseases": [], "_by_alias": {}}
    doc = json.loads(p.read_text())
    by_alias: dict[str, dict[str, Any]] = {}
    for row in doc.get("diseases") or []:
        keys = [row.get("disease_id"), row.get("name"), row.get("mondo_id")]
        keys.extend(list(row.get("aliases") or []))
        for k in keys:
            if not k:
                continue
            by_alias[_norm(str(k))] = row
    doc["_by_alias"] = by_alias
    return doc


def clear_gene_index_cache() -> None:
    load_zero_shot_gene_index.cache_clear()


def lookup_disease_genes(
    query: str,
    *,
    kb_datasources: dict[str, Any] | None = None,
    atlas_disease_id: str | None = None,
) -> dict[str, float]:
    """
    Resolve free-text disease → gene scores from:
      1) zero_shot_disease_genes.json (rare / unseen names)
      2) fused Open Targets rows (atlas ids + optional name match)
    """
    out: dict[str, float] = {}
    key = _norm(query)
    zs = load_zero_shot_gene_index()
    row = (zs.get("_by_alias") or {}).get(key)
    if row:
        for g in row.get("genes") or []:
            if isinstance(g, dict):
                gene = str(g.get("gene") or "").upper()
                score = float(g.get("score") or 0.7)
            else:
                gene = str(g).upper()
                score = 0.7
            if gene:
                out[gene] = max(out.get(gene, 0.0), score)

    ot = (kb_datasources or {}).get("opentargets") or {}
    for r in ot.get("disease_genes") or []:
        did = str(r.get("disease_id") or "")
        name = _norm(r.get("name") or "")
        did_n = _norm(did)
        hit = False
        if atlas_disease_id and did == atlas_disease_id:
            hit = True
        elif key and len(key) >= 5:
            if key == name or key == did_n:
                hit = True
            elif name and (key in name or name in key) and min(len(key), len(name)) >= 6:
                hit = True
        if not hit:
            continue
        for g in r.get("genes") or []:
            gene = str(g.get("gene") or "").upper()
            if not gene:
                continue
            out[gene] = max(out.get(gene, 0.0), float(g.get("score") or 0.55))
    return out


def lookup_zero_shot_record(query: str) -> dict[str, Any] | None:
    zs = load_zero_shot_gene_index()
    return (zs.get("_by_alias") or {}).get(_norm(query))
