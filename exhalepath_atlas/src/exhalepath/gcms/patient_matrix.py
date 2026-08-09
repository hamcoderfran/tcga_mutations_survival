"""Load patient × VOC intensity matrices from Metabolomics Workbench GC-MS studies."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..config import DATA_DIR
from ..ingest.real_breath_corpus import _factor_map, _map_voc, _ms_rows
from .schema import SampleRecord

MW_DIR = DATA_DIR / "datasources" / "metabolomics"

# Bundled studies with disease vs control factors + GC-MS/related intensities
BUNDLED_DIAGNOSTIC_STUDIES: dict[str, dict[str, Any]] = {
    "ST000883": {
        "disease_id": "malaria",
        "disease_name": "Malaria",
        "positive_token": "positive",
        "negative_token": "negative",
        "modality": "gcms",
        "unit": "gcms_intensity",
        "doi": "10.1093/infdis/jiy072",
        "title": "Breathprinting Reveals Malaria-Associated Biomarkers (ST000883)",
        "n_expected_subjects": 35,
    },
    "ST000587": {
        "disease_id": "heart_failure",
        "disease_name": "Heart failure",
        "positive_token": "failure",
        "negative_token": "control",
        "modality": "gcms_ebc",
        "unit": "uM_ebc",
        "doi": None,
        "title": "Heart failure exhaled breath condensate (ST000587)",
        "n_expected_subjects": None,
    },
}


@dataclass
class PatientVOCMatrix:
    study_id: str
    disease_id: str
    disease_name: str
    modality: str
    unit: str
    matrix: pd.DataFrame  # index=subject_id, columns=voc_id
    labels: pd.Series  # 1=disease, 0=control
    samples: list[SampleRecord]
    metadata: dict[str, Any]

    @property
    def subject_ids(self) -> list[str]:
        return list(self.matrix.index.astype(str))

    def log1p(self) -> "PatientVOCMatrix":
        m = self.matrix.copy()
        m = np.log1p(m.clip(lower=0))
        return PatientVOCMatrix(
            study_id=self.study_id,
            disease_id=self.disease_id,
            disease_name=self.disease_name,
            modality=self.modality,
            unit=f"log1p({self.unit})",
            matrix=m,
            labels=self.labels.copy(),
            samples=list(self.samples),
            metadata={**self.metadata, "transform": "log1p"},
        )


def _label_from_factors(
    factors: str, *, positive_token: str, negative_token: str
) -> Optional[int]:
    f = factors.lower()
    # Prefer explicit negative first (avoids "positive" substring traps)
    if negative_token in f and positive_token not in f:
        return 0
    if positive_token in f and negative_token not in f:
        return 1
    if negative_token in f:
        return 0
    if positive_token in f:
        return 1
    return None


def load_mw_patient_matrix(
    study_id: str,
    *,
    feature_map: str = "atlas",
) -> PatientVOCMatrix:
    """Build patient×VOC matrix for a bundled MW study.

    feature_map:
      - ``atlas`` — default ExhalePath catalog mapper
      - ``malaria_lit`` — atlas + Schaber/Berna terpene/thioether/alkane remaps
    """
    meta = BUNDLED_DIAGNOSTIC_STUDIES.get(study_id)
    if meta is None:
        raise ValueError(
            f"Study {study_id} not in bundled diagnostic catalog. "
            f"Known: {sorted(BUNDLED_DIAGNOSTIC_STUDIES)}"
        )
    path = MW_DIR / study_id / "mwtab.json"
    if not path.exists():
        raise FileNotFoundError(path)
    data = json.loads(path.read_text())
    rows = _ms_rows(data)
    fac = _factor_map(study_id)
    if not rows or not fac:
        raise RuntimeError(f"No MS rows/factors for {study_id}")

    df = pd.DataFrame(rows)
    if feature_map == "malaria_lit":
        from .malaria_remap import map_voc_malaria_expanded

        df["voc_id"] = df["metabolite"].map(
            lambda n: map_voc_malaria_expanded(n, _map_voc)
        )
    elif feature_map == "atlas":
        df["voc_id"] = df["metabolite"].map(_map_voc)
    else:
        raise ValueError("feature_map must be atlas|malaria_lit")
    n_before = df["metabolite"].nunique()
    n_mapped = df.dropna(subset=["voc_id"])["metabolite"].nunique()
    df = df.dropna(subset=["voc_id"])
    df["sample"] = df["sample"].astype(str)

    labels: dict[str, int] = {}
    samples: list[SampleRecord] = []
    for sid, factors in fac.items():
        lab = _label_from_factors(
            factors,
            positive_token=meta["positive_token"],
            negative_token=meta["negative_token"],
        )
        if lab is None:
            continue
        labels[sid] = lab
        samples.append(
            SampleRecord(
                sample_id=sid,
                subject_id=sid,
                study_id=study_id,
                label="disease" if lab == 1 else "control",
                disease_id=meta["disease_id"],
                modality=meta["modality"],
                factors_raw=factors,
            )
        )

    keep = [s for s in df["sample"].unique() if s in labels]
    df = df[df["sample"].isin(keep)]
    # median if duplicate metabolite→voc mappings
    piv = df.pivot_table(
        index="sample", columns="voc_id", values="intensity", aggfunc="median"
    )
    piv = piv.reindex(keep)
    # drop all-NaN columns; fill remaining NaN with column median (train-unsafe if used before split — document)
    piv = piv.dropna(axis=1, how="all")
    for col in piv.columns:
        med = float(piv[col].median(skipna=True))
        if np.isnan(med):
            med = 0.0
        piv[col] = piv[col].fillna(med)

    y = pd.Series({s: labels[s] for s in piv.index}, name="label")
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=meta["disease_id"],
        disease_name=meta["disease_name"],
        modality=meta["modality"],
        unit=meta["unit"],
        matrix=piv.astype(float),
        labels=y.loc[piv.index].astype(int),
        samples=[s for s in samples if s.subject_id in set(piv.index)],
        metadata={
            "doi": meta.get("doi"),
            "title": meta.get("title"),
            "n_voc_features": int(piv.shape[1]),
            "n_subjects": int(piv.shape[0]),
            "n_positive": int((y == 1).sum()),
            "n_negative": int((y == 0).sum()),
            "unmapped_metabolites_dropped": True,
            "feature_map": feature_map,
            "n_library_metabolites": int(n_before),
            "n_library_metabolites_mapped": int(n_mapped),
            "voc_ids": list(piv.columns.astype(str)),
        },
    )


def list_bundled_diagnostic_studies() -> list[dict[str, Any]]:
    out = []
    for sid, meta in BUNDLED_DIAGNOSTIC_STUDIES.items():
        path = MW_DIR / sid / "mwtab.json"
        out.append(
            {
                "study_id": sid,
                "available": path.exists(),
                **{k: meta[k] for k in ("disease_id", "disease_name", "modality", "title")},
            }
        )
    return out


def export_patient_matrix_csv(matrix: PatientVOCMatrix, out_dir: Path) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    mp = out_dir / f"{matrix.study_id}_patient_voc_matrix.csv"
    matrix.matrix.to_csv(mp)
    paths["matrix"] = mp
    lp = out_dir / f"{matrix.study_id}_patient_labels.csv"
    matrix.labels.to_csv(lp, header=["label"])
    paths["labels"] = lp
    sp = out_dir / f"{matrix.study_id}_samples.json"
    sp.write_text(
        json.dumps([s.model_dump() for s in matrix.samples], indent=2)
    )
    paths["samples"] = sp
    return paths
