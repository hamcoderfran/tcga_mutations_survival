"""ExhalePath core heads: hybrid, physiology, calibrator, zero-shot."""

from __future__ import annotations

from typing import Any

from ..base import StackModel
from ..types import AspectHit, StackQuery, VOCSignal


def _run_biomarker(ctx: dict[str, Any], query: StackQuery, *, mode: str):
    from ...biomarker import ExhaleBiomarkerEngine

    engine = ctx.get("engine")
    if engine is None:
        engine = ExhaleBiomarkerEngine(use_opentargets=True, reload_knowledge=False)
        ctx["engine"] = engine
    return engine.predict(
        query.disease,
        location=query.location,
        top_n=max(query.top_n, 25),
        genes=query.genes or None,
        description=query.description,
        mode=mode,
        sex=query.sex,
        age_years=query.age,
        smoking_status=query.smoking,
        comorbidities=query.comorbidities or None,
        explain=True,
    )


def _signals_from_report(report) -> list[VOCSignal]:
    out = []
    for i, p in enumerate(report.top_vocs, 1):
        out.append(
            VOCSignal(
                voc_id=p.voc_id,
                log2fc=float(p.log2_fold_change),
                confidence=float(p.confidence),
                rank=i,
                delta_ppb=float(p.delta_ppb),
                evidence=[f"drivers={','.join(p.top_pathway_drivers[:3])}"],
            )
        )
    return out


def _aspects_from_report(report) -> list[AspectHit]:
    aspects: list[AspectHit] = []
    for ps in list(report.result.bundle.pathway_scores or [])[:12]:
        aspects.append(
            AspectHit(
                id=ps.pathway_id,
                name=ps.name,
                score=float(ps.score),
                kind="pathway",
                evidence=["exhalepath_pathway_score"],
            )
        )
    for cs in list(report.result.bundle.cell_states or [])[:10]:
        aspects.append(
            AspectHit(
                id=getattr(cs, "cell_state_id", None) or cs.name,
                name=cs.name,
                score=float(cs.effective_source),
                kind="cell_state",
                evidence=["exhalepath_cell_state"],
            )
        )
    return aspects


class HybridExhalePathModel(StackModel):
    model_id = "exhalepath_hybrid"
    family = "exhalepath"
    aspect = "voc_quantity"
    default_weight = 1.2

    def predict(self, query: StackQuery):
        try:
            report = _run_biomarker(self.ctx, query, mode="hybrid")
            self.ctx["biomarker_report"] = report
            self.ctx["disease_resolved"] = {
                "disease_id": report.disease_id,
                "disease_name": report.disease_name,
                "location": report.location,
            }
            return self._ok(
                voc_signals=_signals_from_report(report),
                aspects=_aspects_from_report(report),
                metadata={"mode": "hybrid", "model_version": report.model_version},
            )
        except Exception as exc:  # noqa: BLE001
            return self._ok(status="degraded", notes=[f"hybrid failed: {exc}"])


class PhysiologyExhalePathModel(StackModel):
    model_id = "exhalepath_physiology"
    family = "exhalepath"
    aspect = "voc_quantity"
    default_weight = 0.9

    def predict(self, query: StackQuery):
        try:
            report = _run_biomarker(self.ctx, query, mode="physiology")
            return self._ok(
                voc_signals=_signals_from_report(report),
                aspects=_aspects_from_report(report),
                metadata={"mode": "physiology"},
            )
        except Exception as exc:  # noqa: BLE001
            return self._ok(status="degraded", notes=[f"physiology failed: {exc}"])


class CalibratorModel(StackModel):
    model_id = "exhalepath_calibrator"
    family = "ml"
    aspect = "voc_quantity"
    default_weight = 0.85

    def predict(self, query: StackQuery):
        # Prefer residual between hybrid and physiology as calibrator contribution
        try:
            hybrid = self.ctx.get("biomarker_report")
            if hybrid is None:
                hybrid = _run_biomarker(self.ctx, query, mode="hybrid")
                self.ctx["biomarker_report"] = hybrid
            phys = _run_biomarker(self.ctx, query, mode="physiology")
            phys_map = {p.voc_id: p for p in phys.top_vocs}
            signals = []
            for p in hybrid.top_vocs:
                q = phys_map.get(p.voc_id)
                if q is None:
                    cal = float(p.log2_fold_change)
                else:
                    cal = float(p.log2_fold_change) - 0.5 * float(q.log2_fold_change)
                signals.append(
                    VOCSignal(
                        voc_id=p.voc_id,
                        log2fc=cal,
                        confidence=max(0.2, float(p.confidence) * 0.9),
                        delta_ppb=float(p.delta_ppb) - (float(q.delta_ppb) if q else 0.0) * 0.5,
                        evidence=["calibrator_residual_hybrid_minus_physio"],
                    )
                )
            return self._ok(voc_signals=signals, metadata={"method": "hybrid_physio_residual"})
        except Exception as exc:  # noqa: BLE001
            return self._ok(status="degraded", notes=[f"calibrator failed: {exc}"])


class ZeroShotModel(StackModel):
    model_id = "zero_shot_mechanism"
    family = "zero_shot"
    aspect = "mechanism"
    default_weight = 0.8

    def predict(self, query: StackQuery):
        from ...knowledge.loader import KnowledgeBase

        kb: KnowledgeBase = self.ctx["kb"]
        d = kb.resolve_disease(
            query.disease,
            location_hint=query.location,
            description=query.description,
            zero_shot=True,
        )
        meta = d.get("metadata") or d.get("zero_shot_meta") or {}
        is_zs = bool(d.get("custom") or str(d.get("disease_id", "")).startswith("custom"))
        # also treat unresolved atlas-adjacent via flags
        zs_flag = bool(meta.get("zero_shot") or is_zs)
        signals = []
        for vid, fc in (d.get("voc_log2fc_prior") or {}).items():
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=float(fc),
                    confidence=0.55 if zs_flag else 0.4,
                    evidence=["zero_shot_or_resolved_prior"],
                )
            )
        aspects = []
        for pw, w in (d.get("pathway_bias") or {}).items():
            aspects.append(
                AspectHit(id=pw, name=pw, score=float(w), kind="pathway", evidence=["zero_shot_bias"])
            )
        for g in (d.get("driver_genes") or query.genes or [])[:12]:
            aspects.append(
                AspectHit(id=g, name=g, score=1.0, kind="gene", evidence=["zero_shot_gene"])
            )
        return self._ok(
            voc_signals=signals,
            aspects=aspects,
            metadata={
                "zero_shot": zs_flag,
                "disease_id": d.get("disease_id"),
                "category": d.get("category"),
                "resolver": meta.get("resolver") or meta.get("zero_shot_mode"),
            },
            notes=["mechanism-first zero-shot head"] if zs_flag else ["atlas-resolved mechanism head"],
        )
