"""PrimeKG-lite graph walks, OmniPath signaling, Open Targets genetics."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

from ..base import StackModel
from ..data import omnipath_edges, primekg_edges
from ..types import AspectHit, StackQuery, VOCSignal


class OpenTargetsGenesModel(StackModel):
    model_id = "opentargets_genes"
    family = "genetics"
    aspect = "genetics"
    default_weight = 0.75

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        did = d.get("disease_id")
        aspects: list[AspectHit] = []
        gene_scores: dict[str, float] = {}
        for g in d.get("driver_genes") or []:
            gene_scores[str(g).upper()] = max(gene_scores.get(str(g).upper(), 0.0), 0.8)
        for row in (kb.datasources.get("opentargets") or {}).get("disease_genes") or []:
            if row.get("disease_id") != did:
                continue
            for g in row.get("genes") or []:
                if isinstance(g, dict):
                    sym = str(g.get("gene") or "").upper()
                    score = float(g.get("score") or 0.5)
                else:
                    sym, score = str(g).upper(), 0.6
                if sym:
                    gene_scores[sym] = max(gene_scores.get(sym, 0.0), score)
        for g in query.genes:
            gene_scores[g.upper()] = max(gene_scores.get(g.upper(), 0.0), 1.0)
        for sym, score in sorted(gene_scores.items(), key=lambda x: -x[1])[:20]:
            aspects.append(
                AspectHit(id=sym, name=sym, score=score, kind="gene", evidence=["opentargets"])
            )
        # project genes → pathway → VOC
        signals: list[VOCSignal] = []
        gene_set = set(gene_scores)
        for pw in kb.pathways.values():
            seeds = {str(x).upper() for x in (pw.get("seed_genes") or [])}
            hit = gene_set & seeds
            if not hit:
                continue
            strength = sum(gene_scores[g] for g in hit) / max(1, len(hit))
            effects = pw.get("voc_effects") or {}
            if isinstance(effects, dict):
                for vid, eff in effects.items():
                    val = float(eff) if isinstance(eff, (int, float)) else float(
                        (eff or {}).get("log2fc") or 0.0
                    )
                    signals.append(
                        VOCSignal(
                            voc_id=vid,
                            log2fc=val * strength,
                            confidence=min(0.9, 0.35 + 0.15 * len(hit)),
                            evidence=[f"genes={','.join(sorted(hit)[:5])}"],
                        )
                    )
        merged: dict[str, VOCSignal] = {}
        for s in signals:
            if s.voc_id not in merged:
                merged[s.voc_id] = s
            else:
                merged[s.voc_id].log2fc += 0.5 * s.log2fc
        return self._ok(
            voc_signals=list(merged.values()),
            aspects=aspects,
            metadata={"n_genes": len(gene_scores), "disease_id": did},
        )


class PrimeKGGraphModel(StackModel):
    model_id = "primekg_graph"
    family = "graph"
    aspect = "knowledge_graph"
    default_weight = 0.9

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location, description=query.description)
        did = d.get("disease_id")
        edges = primekg_edges()
        # adjacency
        fwd: dict[str, list[tuple[str, str, float, dict[str, Any]]]] = defaultdict(list)
        for e in edges:
            fwd[str(e["s"])].append((str(e["o"]), str(e.get("o_type")), float(e.get("w") or 1.0), e))

        # multi-hop: disease -> gene/pathway -> voc (max depth 3)
        voc_scores: dict[str, float] = defaultdict(float)
        gene_scores: dict[str, float] = defaultdict(float)
        pathway_scores: dict[str, float] = defaultdict(float)
        evidence: dict[str, list[str]] = defaultdict(list)

        seeds = [did] if did else []
        # also seed by name token matches on disease nodes via priors
        for e in edges:
            if e.get("s_type") == "disease" and e.get("s") == did and e.get("o_type") == "voc":
                direction = float(e.get("dir") or e.get("w") or 0.0)
                if e.get("dir") is None:
                    direction = float(e.get("w") or 0.0)
                voc_scores[str(e["o"])] += direction
                evidence[str(e["o"])].append("direct_disease_voc")

        # hop1
        hop1 = fwd.get(str(did), [])
        for node, ntype, w, e in hop1:
            if ntype == "gene":
                gene_scores[node] += w
            elif ntype == "pathway":
                pathway_scores[node] += w
            elif ntype == "voc":
                voc_scores[node] += float(e.get("dir") or w)

        # hop2 gene/pathway → pathway/voc
        for gene, gs in list(gene_scores.items()):
            for node, ntype, w, e in fwd.get(gene, []):
                if ntype == "pathway":
                    pathway_scores[node] += gs * w
                elif ntype == "voc":
                    voc_scores[node] += gs * float(e.get("dir") or w)
                    evidence[node].append(f"gene_path:{gene}")
        for pw, ps in list(pathway_scores.items()):
            for node, ntype, w, e in fwd.get(pw, []):
                if ntype == "voc":
                    voc_scores[node] += 0.85 * ps * float(e.get("dir") or w)
                    evidence[node].append(f"pathway_path:{pw}")

        # normalize VOC magnitudes
        if voc_scores:
            mx = max(abs(v) for v in voc_scores.values()) or 1.0
            scale = 1.25 / mx
        else:
            scale = 1.0
        signals = [
            VOCSignal(
                voc_id=vid,
                log2fc=float(score * scale),
                confidence=min(0.9, 0.4 + 0.05 * len(evidence.get(vid, []))),
                evidence=evidence.get(vid, ["primekg_walk"])[:4],
            )
            for vid, score in voc_scores.items()
            if vid in kb.vocs
        ]
        aspects = [
            AspectHit(id=g, name=g, score=float(s), kind="gene", evidence=["primekg"])
            for g, s in sorted(gene_scores.items(), key=lambda x: -x[1])[:15]
        ] + [
            AspectHit(id=p, name=p, score=float(s), kind="pathway", evidence=["primekg"])
            for p, s in sorted(pathway_scores.items(), key=lambda x: -x[1])[:15]
        ]
        status = "ok" if signals or aspects else "degraded"
        return self._ok(
            voc_signals=signals,
            aspects=aspects,
            status=status,
            metadata={"n_edges": len(edges), "disease_id": did, "n_voc_hits": len(signals)},
            notes=[] if signals else ["sparse graph hit; used pathway/gene aspects only"],
        )


class OmniPathSignalingModel(StackModel):
    model_id = "omnipath_signaling"
    family = "network"
    aspect = "signaling"
    default_weight = 0.7

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        d = kb.resolve_disease(query.disease, location_hint=query.location)
        bias = {k: float(v) for k, v in (d.get("pathway_bias") or {}).items()}
        # build crosstalk if file empty
        edges = list(omnipath_edges())
        if not edges:
            pw_genes = {
                pw["pathway_id"]: set(str(g).upper() for g in (pw.get("seed_genes") or []))
                for pw in kb.pathways.values()
            }
            pids = list(pw_genes)
            for i, a in enumerate(pids):
                ga = pw_genes[a]
                if len(ga) < 2:
                    continue
                for b in pids[i + 1 :]:
                    inter = ga & pw_genes[b]
                    if len(inter) >= 1:
                        w = len(inter) / ((len(ga) * max(1, len(pw_genes[b]))) ** 0.5)
                        if w >= 0.08:
                            edges.append({"s": a, "o": b, "w": w, "shared": sorted(inter)[:6]})

        diffused = dict(bias)
        for e in edges:
            a, b, w = e["s"], e["o"], float(e.get("w") or 0.0)
            if a in bias:
                diffused[b] = diffused.get(b, 0.0) + bias[a] * w * 0.5
            if b in bias:
                diffused[a] = diffused.get(a, 0.0) + bias[b] * w * 0.5

        aspects = [
            AspectHit(id=pid, name=pid, score=score, kind="pathway", evidence=["omnipath_diffusion"])
            for pid, score in sorted(diffused.items(), key=lambda x: -abs(x[1]))[:15]
        ]
        signals: list[VOCSignal] = []
        for pid, score in diffused.items():
            pw = kb.pathways.get(pid) or {}
            effects = pw.get("voc_effects") or {}
            if not isinstance(effects, dict):
                continue
            for vid, eff in effects.items():
                val = float(eff) if isinstance(eff, (int, float)) else float(
                    (eff or {}).get("log2fc") or 0.0
                )
                signals.append(
                    VOCSignal(
                        voc_id=vid,
                        log2fc=val * score,
                        confidence=0.5,
                        evidence=[f"signaling:{pid}"],
                    )
                )
        merged: dict[str, VOCSignal] = {}
        for s in signals:
            if s.voc_id not in merged:
                merged[s.voc_id] = s
            else:
                merged[s.voc_id].log2fc += 0.5 * s.log2fc
        return self._ok(
            voc_signals=list(merged.values()),
            aspects=aspects,
            metadata={"n_signaling_edges": len(edges)},
        )
