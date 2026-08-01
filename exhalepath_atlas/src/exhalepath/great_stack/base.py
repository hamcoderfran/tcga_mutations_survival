"""Model adapter protocol."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from .types import ModelOutput, StackQuery


class StackModel(ABC):
    model_id: str = "base"
    family: str = "base"
    aspect: str = "general"
    default_weight: float = 1.0

    def __init__(self, ctx: dict[str, Any]):
        self.ctx = ctx

    @abstractmethod
    def predict(self, query: StackQuery) -> ModelOutput:
        raise NotImplementedError

    def _ok(
        self,
        *,
        voc_signals=None,
        aspects=None,
        metadata=None,
        notes=None,
        weight: float | None = None,
        status: str = "ok",
    ) -> ModelOutput:
        return ModelOutput(
            model_id=self.model_id,
            family=self.family,
            aspect=self.aspect,
            weight=float(weight if weight is not None else self.default_weight),
            status=status,
            voc_signals=list(voc_signals or []),
            aspects=list(aspects or []),
            metadata=dict(metadata or {}),
            notes=list(notes or []),
        )

    def _skip(self, reason: str) -> ModelOutput:
        return self._ok(status="skipped", notes=[reason], weight=0.0)
