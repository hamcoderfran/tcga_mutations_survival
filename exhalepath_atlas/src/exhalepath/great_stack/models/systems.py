"""Human-GEM flux, OPERA physchem, AGORA microbiome, PBPK transport."""

from __future__ import annotations

from collections import defaultdict

from ..base import StackModel
from ..data import agora_routes, flux_proxy, opera_panel
from ..types import AspectHit, StackQuery, VOCSignal


class HumanGEMFluxModel(StackModel):
    model_id = "humangem_flux"
    family = "gsmm"
    aspect = "metabolic_flux"
    default_weight = 0.95

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location, description=query.description)
        flux = flux_proxy()
        cap = flux.get("pathway_capacity") or {}
        producers = flux.get("voc_producers") or {}
        genes = {str(g).upper() for g in (query.genes or d.get("driver_genes") or [])}
        bias = {k: float(v) for k, v in (d.get("pathway_bias") or {}).items()}

        pathway_flux: dict[str, float] = {}
        for pid, meta in cap.items():
            base = float(meta.get("capacity") or 0.5)
            gset = {str(x).upper() for x in (meta.get("genes") or [])}
            cover = len(genes & gset) / max(1, len(gset))
            pathway_flux[pid] = base * (0.4 + bias.get(pid, 0.0)) * (0.7 + 1.2 * cover)

        aspects = [
            AspectHit(
                id=pid,
                name=pid,
                score=score,
                kind="pathway",
                evidence=["humangem_flux_capacity"],
            )
            for pid, score in sorted(pathway_flux.items(), key=lambda x: -x[1])[:15]
            if score > 0
        ]

        voc_flux: dict[str, float] = defaultdict(float)
        for vid, routes in producers.items():
            for r in routes:
                pid = r.get("pathway_id")
                w = float(r.get("weight") or 1.0)
                enzymes = {str(x).upper() for x in (r.get("enzymes") or [])}
                enz_hit = 1.0 + 0.3 * len(genes & enzymes)
                pf = pathway_flux.get(pid, 0.35 + bias.get(pid or "", 0.0))
                voc_flux[vid] += w * pf * enz_hit

        # if sparse producers, fall back to pathway voc_effects * flux
        if len(voc_flux) < 5:
            for pid, pf in pathway_flux.items():
                effects = (kb.pathways.get(pid) or {}).get("voc_effects") or {}
                if isinstance(effects, dict):
                    for vid, eff in effects.items():
                        val = float(eff) if isinstance(eff, (int, float)) else float(
                            (eff or {}).get("log2fc") or 0.0
                        )
                        voc_flux[vid] += val * pf

        if voc_flux:
            mx = max(abs(v) for v in voc_flux.values()) or 1.0
            scale = 1.4 / mx
        else:
            scale = 1.0
        signals = [
            VOCSignal(
                voc_id=vid,
                log2fc=float(v * scale),
                confidence=0.6,
                evidence=["humangem_style_fba_proxy"],
            )
            for vid, v in voc_flux.items()
            if vid in kb.vocs
        ]
        return self._ok(
            voc_signals=signals,
            aspects=aspects,
            metadata={"method": "gene_constrained_pathway_flux_proxy", "n_voc": len(signals)},
            notes=[
                "Uses Human-GEM-style capacity proxy (COBRApy optional). "
                "Install cobrapy + Human-GEM for full FBA."
            ],
        )


class OperaPhyschemModel(StackModel):
    model_id = "opera_physchem"
    family = "qsar"
    aspect = "transport_adme"
    default_weight = 0.85

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        panel = opera_panel()
        report = self.ctx.get("biomarker_report")
        base_fc = dict(d.get("voc_log2fc_prior") or {})
        if report is not None:
            for p in report.top_vocs:
                base_fc[p.voc_id] = float(p.log2_fold_change)
        signals = []
        aspects = []
        for vid, fc in base_fc.items():
            ad = panel.get(vid) or {}
            eff = float(ad.get("exhalation_efficiency") or 0.5)
            fu = float(ad.get("fu_proxy") or 0.5)
            # high efficiency & unbound fraction amplify alveolar appearance
            transport = 0.5 + 0.9 * eff + 0.3 * fu
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=float(fc) * transport,
                    confidence=0.55,
                    evidence=[
                        f"lambda={ad.get('lambda_blood_air')}",
                        f"exhale_eff={eff:.3f}",
                    ],
                )
            )
            aspects.append(
                AspectHit(
                    id=vid,
                    name=ad.get("name") or vid,
                    score=transport,
                    kind="transport",
                    evidence=["opera_style_adme"],
                )
            )
        return self._ok(
            voc_signals=signals,
            aspects=sorted(aspects, key=lambda a: -a.score)[:15],
            metadata={"n_panel": len(panel)},
            notes=["OPERA-style QSPR ADME from partition/pubchem/chembl fused panel"],
        )


class AgoraMicrobiomeModel(StackModel):
    model_id = "agora_microbiome"
    family = "microbiome"
    aspect = "microbiome"
    default_weight = 0.65

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        pack = agora_routes()
        cat = (d.get("category") or "metabolic").lower()
        taxa = list((pack.get("category_taxa") or {}).get(cat) or [])
        # GI / metabolic locations amplify microbiome
        loc = (query.location or d.get("default_site") or "").lower()
        loc_boost = 1.35 if any(x in loc for x in ("gut", "colon", "intestin", "liver", "oral")) else 1.0
        if cat in {"gastrointestinal", "metabolic", "infectious"}:
            loc_boost *= 1.15
        routes = pack.get("routes") or []
        voc_scores: dict[str, float] = defaultdict(float)
        tax_scores: dict[str, float] = defaultdict(float)
        for r in routes:
            tax = r.get("taxon") or "microbe"
            vid = r.get("voc_id")
            if not vid or vid not in kb.vocs:
                continue
            w = float(r.get("weight") or 1.0)
            if taxa and not any(t.lower() in tax.lower() or tax.lower() in t.lower() for t in taxa):
                w *= 0.35
            else:
                w *= 1.0
            w *= loc_boost
            voc_scores[vid] += w
            tax_scores[tax] += w
        if voc_scores:
            mx = max(voc_scores.values()) or 1.0
            scale = 0.9 / mx
        else:
            scale = 1.0
        signals = [
            VOCSignal(
                voc_id=vid,
                log2fc=float(s * scale),
                confidence=0.45,
                evidence=["agora_mvoc_route"],
            )
            for vid, s in voc_scores.items()
        ]
        aspects = []
        # Prefer disease-category taxa as first-class aspects
        for t in taxa:
            aspects.append(
                AspectHit(
                    id=t,
                    name=t,
                    score=float(tax_scores.get(t) or loc_boost),
                    kind="microbe",
                    evidence=["agora_category_taxa"],
                )
            )
        for t, s in sorted(tax_scores.items(), key=lambda x: -x[1]):
            if t in taxa or t.lower() == "microbe":
                continue
            aspects.append(
                AspectHit(
                    id=t,
                    name=t,
                    score=float(s),
                    kind="microbe",
                    evidence=["agora_mvoc_route"],
                )
            )
        return self._ok(
            voc_signals=signals,
            aspects=aspects[:12],
            metadata={"category": cat, "taxa": taxa, "loc_boost": loc_boost},
            notes=["AGORA2-style microbial emission prior (mVOC + curated fallbacks)"],
        )


class PBPKTransportModel(StackModel):
    model_id = "pbpk_transport"
    family = "pbpk"
    aspect = "transport_adme"
    default_weight = 0.8

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        panel = opera_panel()
        phys = kb.physio_constants or {}
        # cardiac output / ventilation proxies by age & sex
        age = float(query.age or 40.0)
        sex = (query.sex or "female").lower()
        vent = 6.0 * (1.05 if sex.startswith("m") else 1.0) * (1.0 - max(0.0, age - 60) * 0.004)
        co = 5.5 * (1.1 if sex.startswith("m") else 1.0) * (1.0 - max(0.0, age - 60) * 0.003)
        # Farhi: C_alv ~ C_blood / (λ + V_A/Q)
        va_q = max(0.5, vent / max(co, 0.1))
        base = dict(d.get("voc_log2fc_prior") or {})
        report = self.ctx.get("biomarker_report")
        if report is not None:
            for p in report.top_vocs:
                base[p.voc_id] = float(p.log2_fold_change)
                if p.physiology and p.physiology.alveolar_fraction is not None:
                    # prefer physiology alveolar fraction when present
                    pass
        signals = []
        for vid, fc in base.items():
            ad = panel.get(vid) or {}
            lam = float(ad.get("lambda_blood_air") or 10.0)
            farhi = 1.0 / (lam + va_q)
            # smoking increases airway oxidative VOC appearance
            smoke = 1.15 if (query.smoking or "").lower() in {"current", "former"} else 1.0
            delivered = float(fc) * (0.35 + 8.0 * farhi) * smoke
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=delivered,
                    confidence=0.58,
                    evidence=[f"farhi_frac={farhi:.4f}", f"VA/Q={va_q:.3f}"],
                )
            )
        aspects = [
            AspectHit(
                id="ventilation_l_min",
                name="Alveolar ventilation proxy",
                score=vent,
                kind="physiology",
                evidence=["pbpk"],
            ),
            AspectHit(
                id="cardiac_output_l_min",
                name="Cardiac output proxy",
                score=co,
                kind="physiology",
                evidence=["pbpk"],
            ),
            AspectHit(
                id="va_q",
                name="V_A/Q",
                score=va_q,
                kind="physiology",
                evidence=["pbpk"],
            ),
        ]
        return self._ok(
            voc_signals=signals,
            aspects=aspects,
            metadata={"vent": vent, "co": co, "va_q": va_q, "phys_keys": list(phys.keys())[:6]},
            notes=["Steady-state Farhi/PBPK-lite alveolar delivery model"],
        )
