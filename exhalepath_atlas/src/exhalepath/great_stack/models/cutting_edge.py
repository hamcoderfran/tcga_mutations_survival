"""Cutting-edge stack heads: phenotype/MONDO, counterfactual null, meta ensemble."""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict

from ..base import StackModel
from ..types import AspectHit, StackQuery, VOCSignal


_PHENOTYPE_THEMES: list[tuple[re.Pattern[str], str, dict[str, float]]] = [
    (
        re.compile(r"sweet|maple|ketoacid|bcaa|branched.?chain", re.I),
        "bcaa_phenotype",
        {"acetone": 0.85, "2_butanone": 0.7, "ammonia": 0.4},
    ),
    (
        re.compile(r"musty|phenylketon|hypopigment|fair skin", re.I),
        "pku_phenotype",
        {"phenol": 0.5, "benzaldehyde": 0.4, "acetone": 0.35},
    ),
    (
        re.compile(r"hyperammon|lethargy|urea.?cycle|otc", re.I),
        "urea_phenotype",
        {"ammonia": 1.0, "dimethyl_amine": 0.45},
    ),
    (
        re.compile(r"tremor|rigidity|bradykinesi|parkinson", re.I),
        "parkinson_phenotype",
        {"dms": 0.4, "acetone": 0.25, "isoprene": -0.2},
    ),
    (
        re.compile(r"hallucin|delusion|psychosis|schizophren", re.I),
        "psychosis_phenotype",
        {"dms": 0.45, "acetone": 0.3, "isoprene": -0.22},
    ),
    (
        re.compile(r"depress|anhedon|suicid", re.I),
        "depression_phenotype",
        {"acetone": 0.28, "isoprene": -0.2, "pentane": 0.25},
    ),
    (
        re.compile(r"inflamm|fever|arthritis|colitis|crohn", re.I),
        "inflammation_phenotype",
        {"hydrogen_sulfide": 0.45, "carbon_disulfide": 0.35, "isoprene": -0.25},
    ),
    (
        re.compile(r"oxidative|lipid.?peroxid|smoke|copd|emphysema", re.I),
        "oxidative_phenotype",
        {"pentane": 0.7, "hexanal": 0.6, "ethane": 0.45},
    ),
    (
        re.compile(r"diabet|polyuria|hyperglyc|ketoacidosis", re.I),
        "diabetes_phenotype",
        {"acetone": 0.9, "isopropanol": 0.45, "isoprene": -0.3},
    ),
    (
        re.compile(r"jaundice|cirrhosis|hepatic|cholestas", re.I),
        "liver_phenotype",
        {"ammonia": 0.55, "dms": 0.5, "limonene": 0.35},
    ),
    (
        re.compile(r"carcinoid|flushing|serotonin", re.I),
        "carcinoid_phenotype",
        {"acetone": 0.3, "ammonia": 0.25},
    ),
    (
        re.compile(r"porphyria|photosensitiv|abdominal.?pain.?crisis", re.I),
        "porphyria_phenotype",
        {"ammonia": 0.35, "pentane": 0.3},
    ),
]


def _text_blob(query: StackQuery) -> str:
    parts = [
        query.disease or "",
        query.description or "",
        " ".join(query.comorbidities or []),
        " ".join(query.genes or []),
    ]
    return " | ".join(p for p in parts if p)


class PhenotypeMondoModel(StackModel):
    """Map free-text phenotypes / MONDO-like descriptions → VOC themes."""

    model_id = "phenotype_mondo"
    family = "phenotype"
    aspect = "phenotype"
    default_weight = 1.05

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        blob = _text_blob(query)
        if len(blob.strip()) < 3:
            return self._skip("no phenotype text")

        hits: list[AspectHit] = []
        voc_acc: dict[str, float] = defaultdict(float)
        for pat, theme, effects in _PHENOTYPE_THEMES:
            if not pat.search(blob):
                continue
            hits.append(
                AspectHit(
                    id=theme,
                    name=theme.replace("_", " "),
                    score=1.0,
                    kind="phenotype",
                    evidence=[f"regex:{pat.pattern[:40]}"],
                )
            )
            for vid, mag in effects.items():
                if vid in kb.vocs:
                    voc_acc[vid] += mag

        # MONDO id from zero-shot resolver if present
        d = kb.resolve_disease(
            query.disease,
            location_hint=query.location,
            description=query.description,
            zero_shot=True,
        )
        meta = d.get("metadata") or {}
        mondo = meta.get("mondo_id") or d.get("mondo_id")
        if mondo:
            hits.append(
                AspectHit(
                    id=str(mondo),
                    name=str(d.get("name") or mondo),
                    score=0.9,
                    kind="mondo",
                    evidence=["resolved_mondo"],
                )
            )

        # strengthen with pathway bias projection when phenotype sparse
        if len(voc_acc) < 3:
            for pid, w in (d.get("pathway_bias") or {}).items():
                effects = (kb.pathways.get(pid) or {}).get("voc_effects") or {}
                if isinstance(effects, dict):
                    for vid, eff in effects.items():
                        val = float(eff) if isinstance(eff, (int, float)) else float(
                            (eff or {}).get("log2fc") or 0.0
                        )
                        voc_acc[vid] += val * float(w) * 0.35

        signals = [
            VOCSignal(
                voc_id=vid,
                log2fc=max(-1.8, min(1.8, fc)),
                confidence=0.58,
                evidence=["phenotype_theme_projection"],
            )
            for vid, fc in voc_acc.items()
            if vid in kb.vocs and abs(fc) >= 0.05
        ]
        return self._ok(
            voc_signals=signals,
            aspects=hits[:12],
            status="ok" if signals or hits else "degraded",
            metadata={
                "n_phenotype_hits": len(hits),
                "mondo_id": mondo,
                "disease_id": d.get("disease_id"),
            },
            notes=["phenotype/MONDO free-text → VOC theme projection"],
        )


class CounterfactualNullModel(StackModel):
    """Conservative null: weak generic noise + pull-to-zero for UQ contrast."""

    model_id = "counterfactual_null"
    family = "null"
    aspect = "null"
    default_weight = 0.35

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        # deterministic pseudo-noise from disease string (reproducible)
        seed = hashlib.sha256((query.disease or "unknown").encode()).hexdigest()
        signals: list[VOCSignal] = []
        # emit near-zero for common breath VOCs — fusion uses this as shrinkage prior
        commons = [
            "acetone",
            "isoprene",
            "pentane",
            "ethane",
            "hexanal",
            "ammonia",
            "ethanol",
            "acetaldehyde",
        ]
        for i, vid in enumerate(commons):
            if vid not in kb.vocs:
                continue
            # tiny deterministic wiggle ±0.08
            nibble = (int(seed[i * 2 : i * 2 + 2], 16) / 255.0 - 0.5) * 0.16
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=float(nibble),
                    confidence=0.25,
                    evidence=["counterfactual_null_baseline"],
                )
            )
        return self._ok(
            voc_signals=signals,
            aspects=[
                AspectHit(
                    id="null_baseline",
                    name="counterfactual null",
                    score=0.1,
                    kind="null",
                    evidence=["sha256_disease_seed"],
                )
            ],
            metadata={"method": "deterministic_null_shrinkage"},
            notes=["null model for epistemic contrast / anti-overclaim"],
        )


class MetaEnsembleModel(StackModel):
    """Second-order meta head: reinforce hybrid VOCs when physio agrees."""

    model_id = "meta_ensemble"
    family = "meta"
    aspect = "voc_quantity"
    default_weight = 1.1

    def predict(self, query: StackQuery):
        report = self.ctx.get("biomarker_report")
        if report is None:
            return self._skip("hybrid report not available for meta ensemble")

        kb = self.ctx["kb"]
        # Build meta vote: hybrid signal × agreement with literature prior if any
        lit = {}
        d = kb.resolve_disease(
            query.disease, location_hint=query.location, description=query.description
        )
        for vid, fc in (d.get("voc_log2fc_prior") or {}).items():
            lit[vid] = float(fc)

        signals: list[VOCSignal] = []
        for p in report.top_vocs[:30]:
            fc = float(p.log2_fold_change)
            conf = float(p.confidence)
            boost = 1.0
            if p.voc_id in lit and lit[p.voc_id] * fc > 0:
                boost = 1.25  # literature agrees on sign
                conf = min(0.95, conf + 0.08)
            elif p.voc_id in lit and lit[p.voc_id] * fc < 0:
                boost = 0.7
                conf = max(0.15, conf - 0.1)
            signals.append(
                VOCSignal(
                    voc_id=p.voc_id,
                    log2fc=fc * boost,
                    confidence=conf,
                    delta_ppb=float(p.delta_ppb) * boost,
                    evidence=["meta_hybrid_literature_agreement"],
                )
            )

        return self._ok(
            voc_signals=signals,
            aspects=[
                AspectHit(
                    id="meta_agreement",
                    name="hybrid×literature meta",
                    score=1.0,
                    kind="meta",
                    evidence=["second_order_ensemble"],
                )
            ],
            metadata={"method": "hybrid_lit_sign_agreement", "n_voc": len(signals)},
            notes=["meta ensemble reinforces atlas-calibrated hybrid when lit agrees"],
        )
