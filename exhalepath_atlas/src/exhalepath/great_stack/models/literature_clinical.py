"""Literature priors, clinical comorbidity, Census, ChEMBL, pathway enrichment."""

from __future__ import annotations

from ..base import StackModel
from ..types import AspectHit, StackQuery, VOCSignal


class LiteraturePriorModel(StackModel):
    model_id = "literature_voc_prior"
    family = "literature"
    aspect = "voc_prior"
    default_weight = 0.7

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location, description=query.description)
        signals = [
            VOCSignal(
                voc_id=vid,
                log2fc=float(fc),
                confidence=0.65,
                evidence=["disease_voc_priors"],
            )
            for vid, fc in (d.get("voc_log2fc_prior") or {}).items()
        ]
        return self._ok(
            voc_signals=signals,
            metadata={"disease_id": d.get("disease_id"), "n": len(signals)},
        )


class CellCensusModel(StackModel):
    model_id = "cell_census"
    family = "single_cell"
    aspect = "cell_state"
    default_weight = 0.75

    def predict(self, query: StackQuery):
        report = self.ctx.get("biomarker_report")
        aspects: list[AspectHit] = []
        signals: list[VOCSignal] = []
        if report is not None:
            cells = list(report.result.bundle.cell_states or [])
            for cs in cells[:12]:
                aspects.append(
                    AspectHit(
                        id=getattr(cs, "cell_state_id", None) or cs.name,
                        name=cs.name,
                        score=float(cs.effective_source),
                        kind="cell_state",
                        evidence=["census_or_atlas_density"],
                    )
                )
            # modulate top physiology VOCs by mean cell activity
            mean_src = sum(float(c.effective_source) for c in cells[:8]) / max(1, min(8, len(cells)))
            for p in report.top_vocs[:20]:
                signals.append(
                    VOCSignal(
                        voc_id=p.voc_id,
                        log2fc=float(p.log2_fold_change) * (0.7 + 0.6 * min(mean_src, 1.5)),
                        confidence=0.55,
                        evidence=[f"cell_activity_scale={mean_src:.3f}"],
                    )
                )
        # also expose atlas cell states for disease site
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        site = (query.location or d.get("default_site") or "").lower()
        for cs in kb.cell_states:
            tissues = [str(t).lower() for t in (cs.get("tissues") or cs.get("sites") or [])]
            if site and tissues and not any(site in t or t in site for t in tissues):
                continue
            aspects.append(
                AspectHit(
                    id=cs.get("cell_state_id") or cs.get("id") or cs.get("name"),
                    name=cs.get("name") or cs.get("cell_state_id") or "cell",
                    score=float(cs.get("baseline_density") or cs.get("density") or 0.2),
                    kind="cell_state",
                    evidence=["cell_state_atlas"],
                )
            )
        return self._ok(voc_signals=signals, aspects=aspects[:20], metadata={"site": site})


class ComorbidityModel(StackModel):
    model_id = "comorbidity_clinical"
    family = "clinical"
    aspect = "clinical"
    default_weight = 0.7

    def predict(self, query: StackQuery):
        if not query.comorbidities:
            return self._skip("no comorbidities provided")
        kb = self.ctx["kb"]
        signals: list[VOCSignal] = []
        aspects: list[AspectHit] = []
        for c in query.comorbidities:
            d = kb.resolve_disease(c, zero_shot=True)
            aspects.append(
                AspectHit(
                    id=d.get("disease_id") or c,
                    name=d.get("name") or c,
                    score=0.65,
                    kind="comorbidity",
                    evidence=["comorbidity_input"],
                )
            )
            for vid, fc in (d.get("voc_log2fc_prior") or {}).items():
                signals.append(
                    VOCSignal(
                        voc_id=vid,
                        log2fc=0.65 * float(fc),
                        confidence=0.5,
                        evidence=[f"comorbidity:{c}"],
                    )
                )
            for pw, w in (d.get("pathway_bias") or {}).items():
                aspects.append(
                    AspectHit(
                        id=pw,
                        name=pw,
                        score=0.65 * float(w),
                        kind="pathway",
                        evidence=[f"comorbidity_pathway:{c}"],
                    )
                )
        # collapse duplicate VOC signals
        merged: dict[str, VOCSignal] = {}
        for s in signals:
            if s.voc_id not in merged:
                merged[s.voc_id] = s
            else:
                merged[s.voc_id].log2fc += s.log2fc
                merged[s.voc_id].evidence.extend(s.evidence)
        return self._ok(
            voc_signals=list(merged.values()),
            aspects=aspects,
            metadata={"comorbidities": query.comorbidities},
        )


class ChemblPharmModel(StackModel):
    model_id = "chembl_pharm"
    family = "pharmacology"
    aspect = "pharmacology"
    default_weight = 0.55

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        priors = kb.chembl_priors or {}
        pw_priors = priors.get("pathways") or {}
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        aspects = []
        signals = []
        for pw, w in (d.get("pathway_bias") or {}).items():
            ch = pw_priors.get(pw) if isinstance(pw_priors, dict) else None
            score = float(w)
            if isinstance(ch, dict):
                score *= 1.0 + float(ch.get("weight") or ch.get("score") or 0.0)
            elif isinstance(ch, (int, float)):
                score *= 1.0 + float(ch)
            aspects.append(
                AspectHit(
                    id=pw,
                    name=pw,
                    score=score,
                    kind="pathway",
                    evidence=["chembl_pathway_prior"],
                )
            )
        voc_phys = priors.get("voc_physchem") or {}
        for vid, fc in list((d.get("voc_log2fc_prior") or {}).items())[:20]:
            boost = 0.0
            if vid in voc_phys:
                boost = 0.15
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=float(fc) * (1.0 + boost),
                    confidence=0.45 + boost,
                    evidence=["chembl_physchem_boost"] if boost else ["chembl_passthrough"],
                )
            )
        return self._ok(voc_signals=signals, aspects=aspects, metadata={"n_chembl_pathways": len(pw_priors)})


class PathwayEnrichmentModel(StackModel):
    model_id = "pathway_enrichment"
    family = "enrichment"
    aspect = "pathways"
    default_weight = 0.8

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location, description=query.description)
        genes = set(g.upper() for g in (query.genes or d.get("driver_genes") or []))
        # Open Targets genes
        ot = kb.datasources.get("opentargets") or {}
        for row in ot.get("disease_genes") or []:
            if row.get("disease_id") == d.get("disease_id"):
                for g in row.get("genes") or []:
                    sym = g.get("gene") if isinstance(g, dict) else g
                    if sym:
                        genes.add(str(sym).upper())
        aspects = []
        signals = []
        for pw in kb.pathways.values():
            seeds = {str(g).upper() for g in (pw.get("seed_genes") or [])}
            if not seeds:
                continue
            overlap = genes & seeds
            if not overlap and pw["pathway_id"] not in (d.get("pathway_bias") or {}):
                continue
            bias = float((d.get("pathway_bias") or {}).get(pw["pathway_id"], 0.0))
            jacc = len(overlap) / max(1, len(seeds))
            score = bias + 2.0 * jacc + 0.15 * len(overlap)
            if score <= 0:
                continue
            aspects.append(
                AspectHit(
                    id=pw["pathway_id"],
                    name=pw.get("name") or pw["pathway_id"],
                    score=score,
                    kind="pathway",
                    evidence=[f"overlap_genes={','.join(sorted(overlap)[:6])}"],
                )
            )
            effects = pw.get("voc_effects") or {}
            if isinstance(effects, dict):
                for vid, eff in effects.items():
                    val = float(eff) if isinstance(eff, (int, float)) else float(
                        (eff or {}).get("log2fc") or (eff or {}).get("weight") or 0.0
                    )
                    signals.append(
                        VOCSignal(
                            voc_id=vid,
                            log2fc=val * (0.5 + 0.5 * min(score, 2.0)),
                            confidence=min(0.85, 0.4 + 0.1 * len(overlap)),
                            evidence=[f"enriched:{pw['pathway_id']}"],
                        )
                    )
        # merge VOC
        merged: dict[str, VOCSignal] = {}
        for s in signals:
            if s.voc_id not in merged:
                merged[s.voc_id] = s
            else:
                merged[s.voc_id].log2fc += s.log2fc
                merged[s.voc_id].confidence = max(merged[s.voc_id].confidence, s.confidence)
        return self._ok(
            voc_signals=list(merged.values()),
            aspects=sorted(aspects, key=lambda a: a.score, reverse=True)[:15],
            metadata={"n_genes": len(genes)},
        )
