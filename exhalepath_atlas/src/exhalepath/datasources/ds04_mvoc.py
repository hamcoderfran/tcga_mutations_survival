"""Priority 4 — Microbial VOC (mVOC) pathways for dysbiosis / gut→breath."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .base import DataSource

# Curated microbe → VOC emission map (mVOC / literature)
MVOC_EMITTERS = [
    {
        "microbe_group": "clostridia_putrefaction",
        "taxa_examples": ["Clostridium", "Clostridioides difficile"],
        "vocs": ["indole", "phenol", "hydrogen_sulfide", "methyl_mercaptan", "dimethyl_disulfide"],
        "pathways": ["microbial_proteolysis_putrefaction"],
        "diseases": ["gut_dysbiosis", "clostridioides_difficile_infection", "inflammatory_bowel_disease"],
    },
    {
        "microbe_group": "fermentative_enterobacteriaceae",
        "taxa_examples": ["Escherichia", "Klebsiella", "Enterobacter"],
        "vocs": ["ethanol", "propanol", "acetaldehyde", "ethyl_acetate", "hydrogen_sulfide"],
        "pathways": ["gut_microbiome_fermentation"],
        "diseases": ["sibo", "gut_dysbiosis", "inflammatory_bowel_disease"],
    },
    {
        "microbe_group": "sulfate_reducers",
        "taxa_examples": ["Desulfovibrio"],
        "vocs": ["hydrogen_sulfide", "methyl_mercaptan", "carbon_disulfide"],
        "pathways": ["microbial_proteolysis_putrefaction", "methionine_transsulfuration"],
        "diseases": ["inflammatory_bowel_disease", "gut_dysbiosis"],
    },
    {
        "microbe_group": "choline_tma_producers",
        "taxa_examples": ["Anaerococcus", "Clostridium"],
        "vocs": ["trimethylamine", "dimethyl_amine"],
        "pathways": ["microbial_proteolysis_putrefaction"],
        "diseases": ["chronic_kidney_disease", "gut_dysbiosis", "heart_disease"],
    },
    {
        "microbe_group": "helicobacter_urease",
        "taxa_examples": ["Helicobacter pylori"],
        "vocs": ["ammonia", "hydrogen_sulfide"],
        "pathways": ["urea_cycle", "gut_microbiome_fermentation"],
        "diseases": ["helicobacter_pylori_infection"],
    },
    {
        "microbe_group": "oral_anaerobes",
        "taxa_examples": ["Porphyromonas", "Fusobacterium"],
        "vocs": ["hydrogen_sulfide", "methyl_mercaptan", "indole"],
        "pathways": ["microbial_proteolysis_putrefaction"],
        "diseases": ["periodontitis"],
    },
]


class MVOCSource(DataSource):
    priority = 4
    key = "mvoc"
    title = "Microbial VOC (mVOC) emitter atlas"
    description = "Microbe groups → VOC/pathway links for dysbiosis breath signals"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        doc = {
            "version": "1.0.0",
            "description": self.description,
            "ref": "http://bioinformatics.charite.de/mvoc/",
            "emitters": MVOC_EMITTERS,
        }
        path = self.write_json("mvoc_emitters.json", doc)
        # VOC → microbe index
        voc_idx: dict[str, list[str]] = {}
        for e in MVOC_EMITTERS:
            for v in e["vocs"]:
                voc_idx.setdefault(v, []).append(e["microbe_group"])
        idx_path = self.write_json("voc_to_microbes.json", voc_idx)
        return {
            "emitters": path,
            "voc_index": idx_path,
            "manifest": self.write_manifest(n_emitter_groups=len(MVOC_EMITTERS)),
        }

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        doc = json.loads((self.out_dir / "mvoc_emitters.json").read_text())
        out = knowledge_dir / "datasource_mvoc.json"
        out.write_text(json.dumps(doc, indent=2))
        # Strengthen disease priors for microbiome VOCs when missing
        priors_path = knowledge_dir / "disease_voc_priors.json"
        n = 0
        if priors_path.exists():
            priors = json.loads(priors_path.read_text())
            by_disease_boost: dict[str, dict[str, float]] = {}
            for e in doc["emitters"]:
                for did in e["diseases"]:
                    for voc in e["vocs"]:
                        by_disease_boost.setdefault(did, {})[voc] = max(
                            by_disease_boost.get(did, {}).get(voc, 0.0), 0.55
                        )
            for d in priors.get("diseases") or []:
                boosts = by_disease_boost.get(d["disease_id"])
                if not boosts:
                    # alias soft match
                    aliases = {a.lower() for a in [d["disease_id"], d.get("name", "")] + list(d.get("aliases") or [])}
                    for did, b in by_disease_boost.items():
                        if did.replace("_", " ") in " ".join(aliases) or did in aliases:
                            boosts = b
                            break
                if not boosts:
                    continue
                vp = dict(d.get("voc_log2fc_prior") or {})
                for voc, val in boosts.items():
                    if voc not in vp:
                        vp[voc] = val
                        n += 1
                d["voc_log2fc_prior"] = vp
            priors_path.write_text(json.dumps(priors, indent=2))
        return {"path": str(out), "n_prior_fills": n}
