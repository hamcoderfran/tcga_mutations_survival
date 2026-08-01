"""Load curated great_stack prior packs (+ mirror under package data)."""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..config import KNOWLEDGE_DIR


def _candidate_dirs() -> list[Path]:
    roots = []
    # packaged alongside knowledge
    roots.append(Path(KNOWLEDGE_DIR).parent / "great_stack")
    # repo data/
    here = Path(__file__).resolve()
    roots.append(here.parents[3] / "data" / "great_stack")  # exhalepath_atlas/data/great_stack
    roots.append(here.parents[2] / "data" / "great_stack")
    # package-local copy
    roots.append(here.parent / "data")
    return roots


@lru_cache(maxsize=1)
def great_stack_dir() -> Path:
    for p in _candidate_dirs():
        if (p / "model_registry.json").exists():
            return p
    # prefer writable repo path
    return _candidate_dirs()[1]


def load_json(name: str) -> dict[str, Any]:
    path = great_stack_dir() / name
    if not path.exists():
        return {}
    return json.loads(path.read_text())


@lru_cache(maxsize=1)
def model_registry() -> dict[str, Any]:
    return load_json("model_registry.json")


@lru_cache(maxsize=1)
def opera_panel() -> dict[str, Any]:
    return (load_json("opera_adme_panel.json").get("panel") or {})


@lru_cache(maxsize=1)
def primekg_edges() -> list[dict[str, Any]]:
    return list(load_json("primekg_lite_graph.json").get("edges") or [])


@lru_cache(maxsize=1)
def omnipath_edges() -> list[dict[str, Any]]:
    return list(load_json("omnipath_signaling_lite.json").get("edges") or [])


@lru_cache(maxsize=1)
def agora_routes() -> dict[str, Any]:
    return load_json("agora_microbiome_routes.json")


@lru_cache(maxsize=1)
def flux_proxy() -> dict[str, Any]:
    return load_json("humangem_flux_proxy.json")
