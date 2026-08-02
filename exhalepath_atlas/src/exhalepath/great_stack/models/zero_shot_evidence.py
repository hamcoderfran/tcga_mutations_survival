"""Zero-shot VOC evidence from literature gene→VOC direction priors.

Uses data/knowledge/zero_shot_disease_genes.json expected_voc_direction
plus mechanism-theme VOC families so unseen diseases still emit real
VOC signals into the fusion layer (not empty priors).
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from ..base import StackModel
from ..types import AspectHit, StackQuery, VOCSignal


_THEME_VOC: dict[str, dict[str, float]] = {
    "bcaa": {
        "acetone": 0.9,
        "2_butanone": 0.75,
        "3_methylbutanal": 0.7,
        "ammonia": 0.45,
    },
    "phenylalanine": {
        "phenol": 0.75,
        "benzaldehyde": 0.55,
        "acetone": 0.35,
    },
    "urea_cycle": {"ammonia": 1.05, "dimethyl_amine": 0.55, "trimethylamine": 0.4},
    "galactose": {"acetaldehyde": 0.55, "ethanol": 0.4, "acetone": 0.3},
    "sphingolipid": {"hexanal": 0.7, "nonanal": 0.55, "pentane": 0.45},
    "neuro": {"dms": 0.45, "acetone": 0.35, "isoprene": -0.25},
    "oxidative": {
        "pentane": 0.85,
        "hexanal": 0.75,
        "ethane": 0.55,
        "nonanal": 0.5,
    },
    "inflammation": {
        "hydrogen_sulfide": 0.55,
        "carbon_disulfide": 0.4,
        "isoprene": -0.3,
    },
    "cns": {"dms": 0.4, "acetone": 0.3, "isoprene": -0.2},
    "metabolic": {"acetone": 0.55, "isoprene": -0.25, "ethanol": 0.35},
    "lipid": {"pentane": 0.7, "hexanal": 0.6, "nonanal": 0.45},
    "microbiome": {
        "indole": 0.55,
        "pyrrole": 0.4,
        "hydrogen_sulfide": 0.5,
        "ammonia": 0.35,
    },
    "heme": {"carbon_disulfide": 0.45, "ammonia": 0.35, "pentane": 0.3},
    "default": {"acetone": 0.4, "isoprene": -0.2, "pentane": 0.35},
}


@lru_cache(maxsize=1)
def _zero_shot_prior() -> dict[str, Any]:
    from ...config import KNOWLEDGE_DIR

    path = Path(KNOWLEDGE_DIR) / "zero_shot_disease_genes.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def _match_prior_entry(disease: str, genes: list[str], prior: dict[str, Any]) -> dict[str, Any] | None:
    d = disease.lower().strip()
    gset = {g.upper() for g in genes}
    best = None
    best_score = 0
    for entry in prior.get("diseases") or []:
        score = 0
        name = str(entry.get("name") or "").lower()
        did = str(entry.get("disease_id") or "").lower().replace("_", " ")
        aliases = [str(a).lower() for a in (entry.get("aliases") or [])]
        if d and (d == name or d == did or d in aliases or name in d or d in name):
            score += 3
        eg = {
            str(g.get("gene") if isinstance(g, dict) else g).upper()
            for g in (entry.get("genes") or [])
        }
        overlap = len(gset & eg)
        score += overlap
        if score > best_score:
            best_score = score
            best = entry
    return best if best_score > 0 else None


def _themes_for(disease: str, genes: list[str], entry: dict[str, Any] | None) -> list[str]:
    themes: list[str] = []
    d = disease.lower()
    gset = {g.upper() for g in genes}
    cat = str((entry or {}).get("category") or "").lower()

    if any(x in d for x in ("maple", "msud", "bcaa")) or gset & {"BCKDHA", "BCKDHB", "DBT", "DLD"}:
        themes.append("bcaa")
    if "phenyl" in d or "pku" in d or gset & {"PAH", "PTS", "QDPR"}:
        themes.append("phenylalanine")
    if any(x in d for x in ("urea", "otc", "citrullin", "argin")) or gset & {
        "OTC",
        "ASS1",
        "ASL",
        "ARG1",
        "CPS1",
    }:
        themes.append("urea_cycle")
    if "galactos" in d or gset & {"GALT", "GALK1", "GALE"}:
        themes.append("galactose")
    if any(x in d for x in ("niemann", "gaucher", "fabry", "tay-sachs", "sphing", "krabbe", "leukodystroph")) or gset & {
        "NPC1",
        "GBA",
        "GLA",
        "HEXA",
        "GALC",
        "ARSA",
    }:
        themes.append("sphingolipid")
    if any(
        x in d
        for x in (
            "parkinson",
            "alzheimer",
            "huntington",
            "als",
            "schizophren",
            "bipolar",
            "depress",
            "autism",
            "epilep",
            "rett",
        )
    ):
        themes.extend(["neuro", "cns"])
    if any(x in d for x in ("inflam", "crohn", "colitis", "rheumat", "lupus", "psoria", "behçet", "behcet", "granulomatosis")):
        themes.extend(["inflammation", "oxidative"])
    if any(x in d for x in ("diabet", "obes", "metabol", "nafld", "nash", "cushing", "addison")):
        themes.extend(["metabolic", "lipid"])
    if any(x in d for x in ("cancer", "carcinoma", "leukemia", "lymphoma", "tumor", "melanoma", "carcinoid")):
        themes.extend(["oxidative", "lipid"])
    if "heme" in d or "porphyr" in d or gset & {"HMBS", "UROD", "CPOX"}:
        themes.append("heme")
    if cat == "metabolic":
        themes.append("metabolic")
    if cat in {"neurological", "neuro"}:
        themes.extend(["neuro", "cns"])
    if not themes:
        themes = ["default", "oxidative"]

    seen: set[str] = set()
    out: list[str] = []
    for t in themes:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def _resolve_voc_id(kb: Any, key: str) -> str | None:
    k = key.lower().replace(" ", "_").replace("-", "_")
    if k in kb.vocs:
        return k
    for vid, meta in kb.vocs.items():
        name = str(meta.get("name") or "").lower().replace(" ", "_").replace("-", "_")
        if name == k or vid.lower() == k:
            return vid
        if k in vid.lower() or k in name:
            return vid
    return None


class ZeroShotEvidenceModel(StackModel):
    model_id = "zero_shot_evidence"
    family = "zero_shot"
    aspect = "voc_prior"
    default_weight = 1.35

    def predict(self, query: StackQuery):
        kb = self.ctx["kb"]
        disease = query.disease or ""
        genes = list(query.genes or [])
        if not disease and not genes and not query.description:
            return self._skip("no disease/genes/description for zero-shot evidence")

        prior = _zero_shot_prior()
        entry = _match_prior_entry(disease, genes, prior)
        themes = _themes_for(disease, genes, entry)

        raw: dict[str, float] = {}
        for t in themes:
            for voc, mag in _THEME_VOC.get(t, {}).items():
                raw[voc] = raw.get(voc, 0.0) + mag * 0.55

        if entry:
            for voc, direction in (entry.get("expected_voc_direction") or {}).items():
                sign = 1.0 if float(direction) >= 0 else -1.0
                base = abs(float(direction)) if abs(float(direction)) > 0.2 else 0.75
                raw[str(voc)] = raw.get(str(voc), 0.0) + sign * base * 0.95
            # pathway_bias → VOC via atlas pathway voc_effects
            for pid, w in (entry.get("pathway_bias") or {}).items():
                effects = (kb.pathways.get(pid) or {}).get("voc_effects") or {}
                if not isinstance(effects, dict):
                    continue
                for vid, eff in effects.items():
                    val = float(eff) if isinstance(eff, (int, float)) else float(
                        (eff or {}).get("log2fc") or (eff or {}).get("weight") or 0.0
                    )
                    raw[vid] = raw.get(vid, 0.0) + val * float(w) * 0.45

        signals: list[VOCSignal] = []
        for key, fc in raw.items():
            vid = _resolve_voc_id(kb, key)
            if vid is None:
                continue
            fc = max(-2.0, min(2.0, float(fc)))
            if abs(fc) < 0.05:
                continue
            signals.append(
                VOCSignal(
                    voc_id=vid,
                    log2fc=fc,
                    confidence=0.68 if entry else 0.52,
                    evidence=[
                        "zero_shot_expected_voc_direction" if entry else "theme_projection",
                        f"themes={','.join(themes[:4])}",
                    ],
                )
            )

        aspects = [
            AspectHit(id=t, name=t, score=1.0, kind="mechanism_theme", evidence=["zero_shot_theme"])
            for t in themes[:8]
        ]
        if entry and entry.get("mondo_id"):
            aspects.append(
                AspectHit(
                    id=str(entry["mondo_id"]),
                    name=str(entry.get("name") or entry["mondo_id"]),
                    score=float(entry.get("confidence") or 0.7),
                    kind="mondo",
                    evidence=["zero_shot_disease_genes"],
                )
            )

        # Evidence channel is always available; regime detection uses
        # zero_shot_mechanism only (see fusion.detect_zero_shot).
        zs_mode = "prior_match" if entry else "theme_transfer"
        status = "ok" if signals else "degraded"
        return self._ok(
            voc_signals=signals,
            aspects=aspects,
            status=status,
            metadata={
                "themes": themes,
                "matched_disease": (entry or {}).get("disease_id"),
                "n_voc": len(signals),
                "evidence_channel": zs_mode,
            },
            notes=["literature/theme VOC projection for zero-shot & novel diseases"],
        )
