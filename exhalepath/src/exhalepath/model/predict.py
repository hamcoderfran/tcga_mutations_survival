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
from ..schemas import DiseaseQuery, PathwayScore, PredictionBundle, VOCPrediction
from .features import build_feature_vector, mechanistic_log2fc


@dataclass
class PredictionResult:
    bundle: PredictionBundle

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame([p.model_dump() for p in self.bundle.predictions])

    def top(self, n: int = 10) -> pd.DataFrame:
        df = self.to_dataframe()
        return df.reindex(df["log2_fold_change"].abs().sort_values(ascending=False).index).head(n)


class ExhalePathPredictor:
    """
    Hybrid mechanistic + calibrated ML predictor for exhaled VOC ppb profiles.

    Cold-start: pathway priors + disease literature offsets.
    Warm-start: per-VOC gradient boosting calibrators trained on the corpus.
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
            # Blend mechanistic prior with ML calibrator
            log2fc = 0.35 * mech_fc + 0.65 * ml_fc
            confidence = min(0.92, confidence + 0.25)
            return log2fc, drivers, confidence

        return mech_fc, drivers, confidence

    def predict(self, query: DiseaseQuery | dict[str, Any]) -> PredictionResult:
        if isinstance(query, dict):
            query = DiseaseQuery.model_validate(query)

        disease = self.kb.resolve_disease(query.disease)
        assoc = self._associated_genes(disease.get("name") or query.disease)
        # Prefer explicit mutated genes; enrich with high-confidence OT genes as soft evidence
        pathway_scores = score_pathways(
            kb=self.kb,
            disease=disease,
            mutated_genes=query.mutated_genes,
            tumor=query.tumor,
            pathway_overrides=query.pathway_overrides,
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
            # Smoking strongly elevates some hydrocarbons
            if query.smoking_status == "current" and voc_id in {"benzene", "toluene", "pentane"}:
                log2fc += 0.45
            if query.smoking_status == "former" and voc_id in {"benzene", "toluene"}:
                log2fc += 0.15

            healthy = float(voc["healthy_ppb_median"])
            # Bound extreme fold-changes for clinical plausibility
            log2fc = float(np.clip(log2fc, -3.5, 4.5))
            pred_ppb = healthy * (2**log2fc)
            # Soft clamp to literature-ish dynamic range, then re-derive fold metrics
            low = float(voc.get("healthy_ppb_low", healthy * 0.1))
            high = float(voc.get("healthy_ppb_high", healthy * 10))
            pred_ppb = float(np.clip(pred_ppb, max(low * 0.2, 1e-9), high * 8))
            if healthy > 0:
                log2fc = float(np.log2(pred_ppb / healthy))
            fold_change = float(pred_ppb / healthy) if healthy > 0 else float("nan")

            sigma = 0.35 * (1.1 - conf)  # log2 space uncertainty
            if query.include_uncertainty:
                ci_low = float(healthy * (2 ** (log2fc - 1.96 * sigma)))
                ci_high = float(healthy * (2 ** (log2fc + 1.96 * sigma)))
                ci_low = float(np.clip(ci_low, max(low * 0.2, 1e-9), high * 8))
                ci_high = float(np.clip(ci_high, max(low * 0.2, 1e-9), high * 8))
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
                )
            )

        preds.sort(key=lambda p: abs(p.log2_fold_change), reverse=True)
        notes = [
            "Predictions fuse literature disease priors, pathway emission coefficients, "
            "and optional ML calibrators trained on a TCGA-scaled corpus.",
            "ppb values are model estimates for research / hypothesis generation — "
            "not clinical diagnoses. Validate against breath GC-MS / PTR-MS cohorts.",
        ]
        if disease.get("_unresolved"):
            notes.append(
                f"Disease '{query.disease}' was not in the curated atlas; "
                "used Open Targets associations (if available) + generic pathway scoring."
            )
        if self._bundle is None:
            notes.append("No trained calibrator found; using mechanistic pathway model only.")

        bundle = PredictionBundle(
            disease_id=disease["disease_id"],
            disease_name=disease["name"],
            query=query,
            pathway_scores=pathway_scores,
            predictions=preds,
            model_version=self._bundle["version"] if self._bundle else "mechanistic-0.1",
            notes=notes,
            metadata={
                "n_vocs": len(preds),
                "n_pathways_scored": len(pathway_scores),
                "n_associated_genes": len(assoc),
                "calibrator_loaded": self._bundle is not None,
            },
        )
        return PredictionResult(bundle=bundle)
