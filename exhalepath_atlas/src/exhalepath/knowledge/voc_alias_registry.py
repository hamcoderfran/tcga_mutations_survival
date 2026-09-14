"""Unified VOC alias + PTR-MS m/z registry (Pillar 5).

Merges:
- ``ingest.real_breath_corpus.NAME_TO_VOC``
- ``datasources.ds01_metabolomics`` name map
- Magdeburg PTR-MS m/z → atlas VOC map
- ``voc_catalog.json`` aliases
- confounder-sensitive VOC seeds

Honesty: mapping aid for partner peak tables — not chemical identity proof.
"""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import DATA_DIR, PACKAGE_ROOT
from ..ingest.real_breath_corpus import NAME_TO_VOC


def _norm(name: str) -> str:
    n = re.sub(r"\([^)]*\)", "", str(name)).strip().lower()
    n = re.sub(r"[^a-z0-9]+", " ", n)
    return re.sub(r"\s+", " ", n).strip()


def _load_json(path: Path) -> Any | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _mz_map_path() -> Path | None:
    for p in (
        DATA_DIR / "real_breath" / "literature_panels" / "magdeburg_ptrms_mz_map.json",
        PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "magdeburg_ptrms_mz_map.json",
        PACKAGE_ROOT / "src" / "exhalepath" / "data" / "real_breath" / "literature_panels" / "magdeburg_ptrms_mz_map.json",
    ):
        if p.exists():
            return p
    return None


def _voc_catalog_path() -> Path | None:
    for p in (
        DATA_DIR / "knowledge" / "voc_catalog.json",
        PACKAGE_ROOT / "data" / "knowledge" / "voc_catalog.json",
        PACKAGE_ROOT / "src" / "exhalepath" / "data" / "knowledge" / "voc_catalog.json",
    ):
        if p.exists():
            return p
    return None


def _ds01_name_map() -> dict[str, str]:
    try:
        from ..datasources.ds01_metabolomics import _NAME_TO_VOC as ds01

        return {str(k).lower(): str(v) for k, v in ds01.items()}
    except Exception:
        return {}


# Confounder-sensitive breath VOCs commonly flagged in diligence (seed registry).
CONFOUNDER_SENSITIVE_SEEDS: dict[str, list[str]] = {
    "smoking": ["acetonitrile", "toluene", "benzene", "xylene", "styrene", "furan", "2_butanone"],
    "oral_microbiome": ["hydrogen_sulfide", "methyl_mercaptan", "indole", "dimethyl_disulfide", "dms"],
    "diet_gut": ["acetone", "isoprene", "ethanol", "methanol", "trimethylamine", "butyric_acid", "acetic_acid"],
    "ambient_exogenous": ["toluene", "benzene", "ethylbenzene", "limonene", "hexane"],
}


def build_voc_alias_registry() -> dict[str, Any]:
    """Build unified alias → voc_id registry with provenance."""
    aliases: dict[str, dict[str, Any]] = {}

    def _put(alias: str, voc_id: str, source: str) -> None:
        key = _norm(alias)
        if not key or not voc_id:
            return
        if key not in aliases:
            aliases[key] = {"voc_id": voc_id, "sources": [source], "raw_alias": alias}
        elif voc_id == aliases[key]["voc_id"]:
            if source not in aliases[key]["sources"]:
                aliases[key]["sources"].append(source)
        else:
            aliases[key].setdefault("conflicts", []).append(
                {"voc_id": voc_id, "source": source}
            )

    for alias, vid in NAME_TO_VOC.items():
        _put(alias, vid, "real_breath_corpus.NAME_TO_VOC")
    for alias, vid in _ds01_name_map().items():
        _put(alias, vid, "ds01_metabolomics._NAME_TO_VOC")

    catalog = _load_json(_voc_catalog_path() or Path())
    if isinstance(catalog, dict):
        for entry in catalog.get("vocs") or catalog.get("compounds") or []:
            if not isinstance(entry, dict):
                continue
            vid = entry.get("voc_id") or entry.get("id")
            if not vid:
                continue
            _put(str(vid).replace("_", " "), str(vid), "voc_catalog.voc_id")
            for a in entry.get("aliases") or []:
                _put(str(a), str(vid), "voc_catalog.aliases")

    mz_entries: list[dict[str, Any]] = []
    mz_doc = _load_json(_mz_map_path() or Path())
    if isinstance(mz_doc, dict):
        for row in mz_doc.get("mz_map") or []:
            mz = row.get("mz")
            vid = row.get("atlas_voc_id")
            mz_entries.append(dict(row))
            if mz is not None and vid:
                _put(f"m/z {mz}", str(vid), "magdeburg_ptrms_mz_map")
                _put(f"mz{mz}", str(vid), "magdeburg_ptrms_mz_map")
                _put(str(mz), str(vid), "magdeburg_ptrms_mz_map")

    return {
        "schema_version": "VocAliasRegistry-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Unified name/m/z → atlas voc_id map for partner ingest diligence. "
            "Conflicts and uncertain Magdeburg channels are retained; do not treat "
            "as chemical identity certification."
        ),
        "n_aliases": len(aliases),
        "aliases": aliases,
        "mz_map": mz_entries,
        "confounder_sensitive": CONFOUNDER_SENSITIVE_SEEDS,
        "clinical_claim": False,
    }


def map_name_to_voc(name: str, registry: dict[str, Any] | None = None) -> str | None:
    """Resolve a free-text / m/z label to atlas voc_id via the registry."""
    reg = registry or build_voc_alias_registry()
    key = _norm(name)
    if not key:
        return None
    hit = (reg.get("aliases") or {}).get(key)
    if hit and not hit.get("conflicts"):
        return str(hit["voc_id"])
    if hit and hit.get("voc_id"):
        return str(hit["voc_id"])
    # snake_case passthrough if already an atlas id key
    snake = key.replace(" ", "_")
    for alias, meta in (reg.get("aliases") or {}).items():
        if meta.get("voc_id") == snake or alias == snake:
            return snake
    return None


def audit_column_names(
    names: list[str],
    *,
    registry: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Report mapped vs unmapped rates for partner feature column names."""
    reg = registry or build_voc_alias_registry()
    skip = {
        "sample_id",
        "sample",
        "subject_id",
        "subject",
        "id",
        "label",
        "class",
        "group",
        "disease",
        "status",
    }
    cols = [str(n) for n in names if _norm(n) not in skip]
    mapped: list[dict[str, str]] = []
    unmapped: list[str] = []
    for c in cols:
        vid = map_name_to_voc(c, reg)
        if vid:
            mapped.append({"column": c, "voc_id": vid})
        else:
            unmapped.append(c)
    n = len(cols)
    rate = (len(unmapped) / n) if n else None
    return {
        "n_columns": n,
        "n_mapped": len(mapped),
        "n_unmapped": len(unmapped),
        "unmapped_rate": rate,
        "mapped": mapped,
        "unmapped": unmapped,
        "honesty": "Unmapped-rate is a diligence metric for partner tables, not data quality alone.",
    }


def write_voc_alias_registry(
    *,
    out_dir: Path | None = None,
) -> dict[str, Any]:
    report = build_voc_alias_registry()
    out = Path(out_dir) if out_dir else (DATA_DIR / "knowledge" / "voc_alias_registry")
    out.mkdir(parents=True, exist_ok=True)
    (out / "VOC_ALIAS_REGISTRY.json").write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# VOC alias / m/z registry",
        "",
        f"Generated: {report['generated_utc']}",
        "",
        f"> {report['honesty']}",
        "",
        f"- aliases: **{report['n_aliases']}**",
        f"- Magdeburg m/z rows: **{len(report.get('mz_map') or [])}**",
        "",
        "## Confounder-sensitive seeds",
        "",
    ]
    for k, vocs in (report.get("confounder_sensitive") or {}).items():
        lines.append(f"- `{k}`: {', '.join(f'`{v}`' for v in vocs)}")
    lines += [
        "",
        "## Regenerate",
        "",
        "```bash",
        "voc audit-voc-aliases",
        "```",
        "",
    ]
    (out / "VOC_ALIAS_REGISTRY.md").write_text("\n".join(lines))

    # confounder registry companion seed
    conf = {
        "schema_version": "ConfounderRegistry-1.0",
        "generated_utc": report["generated_utc"],
        "honesty": (
            "Seed registry of breath VOCs often sensitive to smoking, oral microbiome, "
            "diet/gut, or ambient exposure. Use with voc eval-confounder-ptr — not disease AUROC."
        ),
        "strata": report.get("confounder_sensitive"),
        "related_reports": [
            "knowledge/confounder_ptr/CONFOUNDER_PTR.json",
            "knowledge/voc_alias_registry/VOC_ALIAS_REGISTRY.json",
        ],
        "clinical_claim": False,
    }
    (out / "CONFOUNDER_REGISTRY.json").write_text(json.dumps(conf, indent=2) + "\n")

    pkg = PACKAGE_ROOT / "data" / "knowledge" / "voc_alias_registry"
    if out.resolve() != pkg.resolve() and (PACKAGE_ROOT / "data" / "knowledge").exists():
        pkg.mkdir(parents=True, exist_ok=True)
        for name in (
            "VOC_ALIAS_REGISTRY.json",
            "VOC_ALIAS_REGISTRY.md",
            "CONFOUNDER_REGISTRY.json",
        ):
            (pkg / name).write_text((out / name).read_text())
    return report


__all__ = [
    "CONFOUNDER_SENSITIVE_SEEDS",
    "audit_column_names",
    "build_voc_alias_registry",
    "map_name_to_voc",
    "write_voc_alias_registry",
]
