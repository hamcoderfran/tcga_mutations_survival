from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from ..config import MODELS_DIR
from ..ingest.opentargets import OpenTargetsClient
from ..knowledge.comorbidity import (
    merge_disease_with_comorbidities,
    resolve_comorbid_diseases,
)
from ..knowledge.disease_gene_index import lookup_disease_genes
from ..knowledge.loader import KnowledgeBase, default_knowledge
from ..pathways.score import score_pathways
from ..physio.census_fractions import census_cell_state_fractions_for_disease
from ..physio.engine import PhysiologyEngine
from ..schemas import DiseaseQuery, PathwayScore, PredictionBundle, VOCPrediction
from .features import build_feature_vector, mechanistic_log2fc


@dataclass
class PredictionResult:
    bundle: PredictionBundle

    def to_dataframe(self) -> pd.DataFrame:
        rows = []
        for p in self.bundle.predictions:
            d = p.model_dump()
            phys = d.pop("physiology", None) or {}
            d["tissue_production"] = phys.get("tissue_production")
            d["blood_delivery"] = phys.get("blood_delivery")
            d["alveolar_fraction"] = phys.get("alveolar_fraction")
            d["lambda_blood_air"] = phys.get("lambda_blood_air")
            d["cell_states"] = ",".join(phys.get("contributing_cell_states") or [])
            d["chains"] = ",".join(phys.get("contributing_chains") or [])
            rows.append(d)
        return pd.DataFrame(rows)

    def top(self, n: int = 10) -> pd.DataFrame:
        df = self.to_dataframe()
        return df.reindex(df["log2_fold_change"].abs().sort_values(ascending=False).index).head(n)

    def cell_states_frame(self) -> pd.DataFrame:
        return pd.DataFrame([c.model_dump() for c in self.bundle.cell_states])


class ExhalePathPredictor:
    """
    Hybrid mechanistic + physiology + calibrated ML predictor for exhaled VOC ppb.

    Physiology stack:
      pathway chains × affected cell states (density×activity)
      → tissue production → blood transfer → Farhi alveolar release → ppb
    """

    def __init__(
        self,
        knowledge: KnowledgeBase | None = None,
        model_path: Path | None = None,
        use_opentargets: bool = True,
    ):
        self.kb = knowledge or default_knowledge()
        self.use_opentargets = use_opentargets
        self.model_path = Path(model_path or MODELS_DIR / "voc_calibrator.joblib")
        self._bundle: dict[str, Any] | None = None
        if self.model_path.exists():
            self._bundle = joblib.load(self.model_path)
        self.physio = PhysiologyEngine(self.kb)

    def _associated_genes(self, disease_name: str) -> dict[str, float]:
        if not self.use_opentargets:
            return {}
        try:
            df = OpenTargetsClient().disease_targets(disease_name)
            if df.empty:
                return {}
            return {r.gene_symbol: float(r.score) for r in df.itertuples()}
        except Exception:  # noqa: BLE001
            return {}

    def _offline_associated_genes(self, disease: dict) -> dict[str, float]:
        """Use fused Open Targets / GDC fragments when live OT is disabled."""
        out: dict[str, float] = {}
        did = disease.get("disease_id")
        ot = (self.kb.datasources.get("opentargets") or {}) if self.kb.datasources else {}
        for row in ot.get("disease_genes") or []:
            if row.get("disease_id") != did:
                continue
            for g in row.get("genes") or []:
                gene = str(g.get("gene") or "").upper()
                if not gene:
                    continue
                out[gene] = max(float(out.get(gene, 0.0)), float(g.get("score") or 0.55))
        gdc = (self.kb.datasources.get("gdc") or {}) if self.kb.datasources else {}
        for row in gdc.get("projects") or []:
            if row.get("disease_id") != did:
                continue
            for gene in row.get("driver_genes") or []:
                gu = str(gene).upper()
                out[gu] = max(float(out.get(gu, 0.0)), 0.7)
        # Curated disease driver_genes as soft associations
        for gene in disease.get("driver_genes") or []:
            gu = str(gene).upper()
            out[gu] = max(float(out.get(gu, 0.0)), 0.65)
        return out

    def _predict_log2fc(
        self,
        *,
        disease: dict,
        pathway_scores: list[PathwayScore],
        query: DiseaseQuery,
        voc_id: str,
    ) -> tuple[float, list[str], float]:
        mech_fc, drivers = mechanistic_log2fc(
            kb=self.kb,
            disease=disease,
            pathway_scores=pathway_scores,
            voc_id=voc_id,
            tumor=query.tumor,
        )
        mech_conf = float(disease.get("_mechanism_confidence") or 0.0)
        if disease.get("_unresolved"):
            # Zero-shot: scale confidence by mechanism evidence strength
            confidence = 0.28 + 0.45 * min(max(mech_conf, 0.0), 1.0)
        else:
            confidence = 0.62

        if self._bundle and voc_id in self._bundle.get("models", {}):
            model = self._bundle["models"][voc_id]
        elif self._bundle and "__global__" in self._bundle.get("models", {}):
            model = self._bundle["models"]["__global__"]
        else:
            model = None

        if model is not None:
            feats = build_feature_vector(
                kb=self.kb,
                disease=disease,
                pathway_scores=pathway_scores,
                tumor=query.tumor,
                voc_id=voc_id,
                age_years=query.age_years,
                sex=query.sex,
                smoking_status=query.smoking_status,
            )
            cols = self._bundle["feature_columns"]
            x = pd.DataFrame([[feats.get(c, 0.0) for c in cols]], columns=cols)
            ml_fc = float(model.predict(x)[0])
            category = (disease.get("category") or "").lower()
            mech_heavy = category in {
                "microbiome",
                "neurological",
                "neurodegenerative",
                "inflammatory",
                "metabolic",
                "infectious",
                "pulmonary",
                "cardiovascular",
            } or bool(disease.get("_unresolved"))
            w_mech = 0.75 if mech_heavy else 0.35
            log2fc = w_mech * mech_fc + (1.0 - w_mech) * ml_fc
            confidence = min(0.92, confidence + (0.12 if mech_heavy else 0.25))
            return log2fc, drivers, confidence

        return mech_fc, drivers, confidence

    def predict(self, query: DiseaseQuery | dict[str, Any]) -> PredictionResult:
        if isinstance(query, dict):
            query = DiseaseQuery.model_validate(query)

        loc_hint = None
        if query.tumor and query.tumor.primary_site:
            loc_hint = query.tumor.primary_site
        disease = self.kb.resolve_disease(
            query.disease,
            location_hint=loc_hint,
            description=query.description,
            copy_voc_priors=bool(query.copy_voc_priors_from_neighbor),
        )
        if query.comorbidities:
            comorbid = resolve_comorbid_diseases(self.kb, list(query.comorbidities))
            disease = merge_disease_with_comorbidities(
                disease,
                comorbid,
                weight=float(query.comorbidity_weight),
            )
        assoc = self._associated_genes(disease.get("name") or query.disease)
        # Offline OT/GDC/driver priors always available (even when live OT is off)
        for g, s in self._offline_associated_genes(disease).items():
            assoc[g] = max(float(assoc.get(g, 0.0)), float(s))
        # Free-text / rare-disease gene index (beyond atlas disease_id keys)
        for g, s in lookup_disease_genes(
            disease.get("name") or query.disease,
            kb_datasources=self.kb.datasources,
            atlas_disease_id=None
            if disease.get("_unresolved")
            else disease.get("disease_id"),
        ).items():
            assoc[g] = max(float(assoc.get(g, 0.0)), float(s))
        for g, s in (disease.get("_associated_gene_scores") or {}).items():
            assoc[str(g).upper()] = max(
                float(assoc.get(str(g).upper(), 0.0)), float(s)
            )
        # Also pull OT associations for comorbidities (soft)
        if query.comorbidities and self.use_opentargets:
            for c in disease.get("_comorbidities") or []:
                cname = c.get("name")
                if not cname:
                    continue
                extra = self._associated_genes(cname)
                for g, s in extra.items():
                    assoc[g] = max(float(assoc.get(g, 0.0)), 0.55 * float(s))
        # Offline comorbidity gene priors
        if query.comorbidities:
            for c in disease.get("_comorbidities") or []:
                cid = c.get("disease_id")
                if not cid:
                    continue
                cdis = self.kb.diseases.get(cid) or {"disease_id": cid}
                for g, s in self._offline_associated_genes(cdis).items():
                    assoc[g] = max(float(assoc.get(g, 0.0)), 0.5 * float(s))
        pathway_scores = score_pathways(
            kb=self.kb,
            disease=disease,
            mutated_genes=query.mutated_genes,
            tumor=query.tumor,
            pathway_overrides=query.pathway_overrides,
            associated_genes=assoc,
        )

        # Lean into Census single-cell fractions when user did not override them
        if not query.cell_state_fractions:
            census_fracs = census_cell_state_fractions_for_disease(disease.get("disease_id", ""))
            if census_fracs:
                query = query.model_copy(update={"cell_state_fractions": census_fracs})

        physio_result = None
        if query.mode in {"physiology", "hybrid"}:
            physio_result = self.physio.run(
                disease=disease,
                query=query,
                pathway_scores=pathway_scores,
                associated_genes=assoc,
            )

        preds: list[VOCPrediction] = []
        for voc_id, voc in self.kb.vocs.items():
            log2fc, drivers, conf = self._predict_log2fc(
                disease=disease,
                pathway_scores=pathway_scores,
                query=query,
                voc_id=voc_id,
            )
            healthy = float(voc["healthy_ppb_median"])
            log2fc = float(np.clip(log2fc, -3.5, 4.5))
            legacy_ppb = healthy * (2**log2fc)

            trace = physio_result.traces.get(voc_id) if physio_result else None
            if query.mode == "legacy" or trace is None:
                pred_ppb = legacy_ppb
            elif query.mode == "physiology":
                pred_ppb = trace.predicted_ppb
                conf = min(0.9, conf + 0.08)
                drivers = list(dict.fromkeys(drivers + trace.contributing_chains[:3]))
            else:  # hybrid
                category = (disease.get("category") or "").lower()
                w_phys = 0.7 if category in {"metabolic", "microbiome", "neurological"} else 0.55
                pred_ppb = PhysiologyEngine.blend_with_legacy(
                    physio_ppb=trace.predicted_ppb,
                    legacy_ppb=legacy_ppb,
                    healthy_ppb=healthy,
                    weight_physio=w_phys,
                )
                # Do not let an under-calibrated physio branch cancel legacy direction
                if legacy_ppb >= healthy and trace.predicted_ppb < healthy:
                    pred_ppb = max(pred_ppb, legacy_ppb)
                if legacy_ppb < healthy and trace.predicted_ppb > healthy:
                    pred_ppb = min(pred_ppb, legacy_ppb)
                conf = min(0.93, conf + 0.1)
                drivers = list(dict.fromkeys(drivers + trace.contributing_chains[:3]))

            # Exogenous / demographic multipliers AFTER hybrid blend so smoking/age
            # effects are not erased by physio↔legacy saturation (lit-compare fix).
            exo_log2 = 0.0
            if query.smoking_status == "current" and voc_id in {
                "benzene",
                "toluene",
                "pentane",
                "ethylbenzene",
            }:
                exo_log2 += 0.45
            elif query.smoking_status == "former" and voc_id in {"benzene", "toluene"}:
                exo_log2 += 0.15
            if query.age_years is not None and voc_id in {"pentane", "ethane", "hexanal"}:
                exo_log2 += 0.004 * (float(query.age_years) - 50.0)
            if query.sex == "male" and voc_id == "isoprene":
                exo_log2 += 0.05
            if query.sex == "female" and voc_id in {"acetone", "isopropanol"}:
                exo_log2 += 0.04
            if abs(exo_log2) > 1e-12:
                pred_ppb = float(pred_ppb) * (2.0**exo_log2)
                if exo_log2 > 0 and query.smoking_status in {"current", "former"}:
                    drivers = list(dict.fromkeys([*drivers, "smoking_exposure"]))
                if query.age_years is not None and voc_id in {"pentane", "ethane", "hexanal"}:
                    drivers = list(dict.fromkeys([*drivers, "age_oxidative"]))

            low = float(voc.get("healthy_ppb_low", healthy * 0.1))
            high = float(voc.get("healthy_ppb_high", healthy * 10))
            # Headroom so exogenous multipliers are not clipped away
            upper = high * 8.0 * (2.0 ** max(0.0, exo_log2))
            pred_ppb = float(np.clip(pred_ppb, max(low * 0.2, 1e-9), upper))
            if healthy > 0:
                log2fc = float(np.log2(pred_ppb / healthy))
            fold_change = float(pred_ppb / healthy) if healthy > 0 else float("nan")

            sigma = 0.35 * (1.1 - conf)
            if query.include_uncertainty:
                ci_low = float(np.clip(healthy * (2 ** (log2fc - 1.96 * sigma)), max(low * 0.2, 1e-9), high * 8))
                ci_high = float(np.clip(healthy * (2 ** (log2fc + 1.96 * sigma)), max(low * 0.2, 1e-9), high * 8))
            else:
                ci_low = None
                ci_high = None

            preds.append(
                VOCPrediction(
                    voc_id=voc_id,
                    name=voc["name"],
                    cas=voc.get("cas"),
                    healthy_ppb=healthy,
                    predicted_ppb=pred_ppb,
                    delta_ppb=pred_ppb - healthy,
                    log2_fold_change=log2fc,
                    fold_change=fold_change,
                    ci_low_ppb=ci_low,
                    ci_high_ppb=ci_high,
                    top_pathway_drivers=drivers,
                    confidence=conf,
                    physiology=trace,
                )
            )

        preds.sort(key=lambda p: abs(p.log2_fold_change), reverse=True)
        notes = [
            "Predictions fuse pathway chains, affected cell-state density/activity, "
            "blood transfer, Farhi alveolar release, literature priors, and optional ML calibrators.",
            "ppb values are model estimates for research / hypothesis generation — "
            "not clinical diagnoses. Validate against breath GC-MS / PTR-MS cohorts.",
        ]
        if physio_result:
            notes.extend(physio_result.notes)
        if disease.get("_comorbidities"):
            names = [
                c.get("name") or c.get("disease_id")
                for c in disease["_comorbidities"]
            ]
            notes.append(
                "Comorbidities fused into pathway bias, VOC priors, and cell-state "
                f"modulation: {', '.join(str(n) for n in names)} "
                f"(weight={query.comorbidity_weight})."
            )
        if disease.get("_unresolved"):
            mode = disease.get("_zero_shot_mode") or "near_healthy_fallback"
            mech_c = float(disease.get("_mechanism_confidence") or 0.0)
            notes.append(
                f"Disease '{query.disease}' was not in the curated atlas "
                f"(zero-shot mode={mode}, mechanism_confidence={mech_c:.2f}). "
                "Used token/category cues, ontology nearest-neighbor pathway transfer "
                "(no VOC prior copy by default), and disease→gene index when available."
            )
            if mode == "near_healthy_fallback":
                notes.append(
                    "Uncertainty flag: no usable mechanism evidence — panel stays near-healthy. "
                    "Supply genes, pathway overrides, tissue/cell fractions, or a phenotype description."
                )
            elif mech_c < 0.55:
                notes.append(
                    "Uncertainty flag: weak zero-shot mechanism evidence — interpret directional "
                    "shifts cautiously."
                )
            if disease.get("_cell_state_donor_ids"):
                notes.append(
                    "Cell-state priors transferred from ontology neighbors: "
                    + ", ".join(disease["_cell_state_donor_ids"])
                )
        if query.mutated_genes or query.pathway_overrides or query.cell_state_fractions:
            notes.append(
                "User mechanism inputs applied "
                f"(genes={len(query.mutated_genes)}, "
                f"pathway_overrides={len(query.pathway_overrides)}, "
                f"cell_fractions={len(query.cell_state_fractions)})."
            )
        if self._bundle is None:
            notes.append("No trained calibrator found; using mechanistic/physiology models only.")

        bundle = PredictionBundle(
            disease_id=disease["disease_id"],
            disease_name=disease["name"],
            query=query,
            pathway_scores=pathway_scores,
            predictions=preds,
            cell_states=physio_result.cell_states if physio_result else [],
            model_version=(
                f"{self._bundle['version']}+physio-1.0"
                if self._bundle
                else "physio-1.0"
            ),
            notes=notes,
            metadata={
                "n_vocs": len(preds),
                "n_pathways_scored": len(pathway_scores),
                "n_associated_genes": len(assoc),
                "calibrator_loaded": self._bundle is not None,
                "mode": query.mode,
                "n_cell_states": len(physio_result.cell_states) if physio_result else 0,
                "n_chains": len(physio_result.chain_fluxes) if physio_result else 0,
                "comorbidities": disease.get("_comorbidities") or [],
                "comorbidity_ids": disease.get("_comorbid_ids") or [],
                "zero_shot": bool(disease.get("_zero_shot")),
                "zero_shot_mode": disease.get("_zero_shot_mode"),
                "mechanism_confidence": float(disease.get("_mechanism_confidence") or 0.0),
                "zero_shot_evidence": disease.get("_zero_shot_evidence") or [],
                "cell_state_donor_ids": disease.get("_cell_state_donor_ids") or [],
                "category": disease.get("category"),
                "default_site": disease.get("default_site"),
                "unresolved": bool(disease.get("_unresolved")),
            },
        )
        return PredictionResult(bundle=bundle)
