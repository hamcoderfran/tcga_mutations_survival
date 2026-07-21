from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd

from ..config import MODELS_DIR
from ..ingest.opentargets import OpenTargetsClient
from ..knowledge.loader import KnowledgeBase, default_knowledge
from ..pathways.score import score_pathways
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
        confidence = 0.45 if disease.get("_unresolved") else 0.62

        if self._bundle and voc_id in self._bundle.get("models", {}):
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
            ml_fc = float(self._bundle["models"][voc_id].predict(x)[0])
            category = (disease.get("category") or "").lower()
            mech_heavy = category in {
                "microbiome",
                "neurological",
                "neurodegenerative",
                "inflammatory",
                "metabolic",
            } or bool(disease.get("_unresolved"))
            w_mech = 0.75 if mech_heavy else 0.35
            log2fc = w_mech * mech_fc + (1.0 - w_mech) * ml_fc
            confidence = min(0.92, confidence + (0.12 if mech_heavy else 0.25))
            return log2fc, drivers, confidence

        return mech_fc, drivers, confidence

    def predict(self, query: DiseaseQuery | dict[str, Any]) -> PredictionResult:
        if isinstance(query, dict):
            query = DiseaseQuery.model_validate(query)

        disease = self.kb.resolve_disease(query.disease)
        assoc = self._associated_genes(disease.get("name") or query.disease)
        pathway_scores = score_pathways(
            kb=self.kb,
            disease=disease,
            mutated_genes=query.mutated_genes,
            tumor=query.tumor,
            pathway_overrides=query.pathway_overrides,
            associated_genes=assoc,
        )

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
            if query.smoking_status == "current" and voc_id in {"benzene", "toluene", "pentane"}:
                log2fc += 0.45
            if query.smoking_status == "former" and voc_id in {"benzene", "toluene"}:
                log2fc += 0.15

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
                # Do not let an under-calibrated physio branch cancel a positive legacy signal
                if legacy_ppb >= healthy and trace.predicted_ppb < healthy:
                    pred_ppb = max(pred_ppb, legacy_ppb)
                conf = min(0.93, conf + 0.1)
                drivers = list(dict.fromkeys(drivers + trace.contributing_chains[:3]))

            low = float(voc.get("healthy_ppb_low", healthy * 0.1))
            high = float(voc.get("healthy_ppb_high", healthy * 10))
            pred_ppb = float(np.clip(pred_ppb, max(low * 0.2, 1e-9), high * 8))
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
        if disease.get("_unresolved"):
            notes.append(
                f"Disease '{query.disease}' was not in the curated atlas; "
                "used Open Targets associations (if available) + generic pathway scoring."
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
            },
        )
        return PredictionResult(bundle=bundle)
