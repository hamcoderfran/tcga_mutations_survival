from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any


class DataSource(ABC):
    """One priority datasource (1–12)."""

    priority: int
    key: str
    title: str
    description: str

    def __init__(self, root: Path):
        self.root = Path(root)
        self.out_dir = self.root / "data" / "datasources" / self.key
        self.out_dir.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        """Download / build local tables. Returns artifact paths."""

    @abstractmethod
    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        """Write/merge knowledge JSON fragments used at predict time."""

    def write_json(self, name: str, payload: Any) -> Path:
        path = self.out_dir / name
        path.write_text(json.dumps(payload, indent=2, default=str))
        return path

    def write_manifest(self, **kwargs: Any) -> Path:
        doc = {
            "priority": self.priority,
            "key": self.key,
            "title": self.title,
            **kwargs,
        }
        return self.write_json("manifest.json", doc)
