"""External / larger malaria breath cohort adapters.

Public reality (2026):
- ST000883 (Schaber 2018, MW PR000612) is the only open *quantified* pediatric
  malaria breath intensity table on Metabolomics Workbench.
- Berna CHMI thioether work has CSIRO QTOF raw deposits (Agilent .D) — not a
  ready patient×VOC CSV; adapter records the DOI and fetch instructions.
- 2024 Malawi reproducibility (J Infect Dis jiae323) is not deposited as an
  open intensity matrix at time of writing.

This module therefore:
1. Documents / optionally fetches CSIRO CHMI metadata for diligence
2. Builds a *virtual external* protocol: learning curves + multi-seed nested
   CIs on the remapped ST000883 matrix (JBR-style n≥50 guidance)
3. Loads a user-supplied second mwtab/CSV when available (`--external-matrix`)
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..config import DATA_DIR
from .patient_matrix import PatientVOCMatrix
from .schema import SampleRecord

EXTERNAL_CATALOG = {
    "ST000883": {
        "role": "primary_open_intensity",
        "doi": "10.1093/infdis/jiy072",
        "n_subjects": 35,
        "status": "bundled",
    },
    "CSIRO_CHMI_QTOF": {
        "role": "external_raw_chmi",
        "doi": "10.25919/5b5b7530a39f4",
        "url": "https://data.csiro.au/collection/csiro:33843",
        "n_subjects": 7,
        "status": "raw_agilent_d_only",
        "note": (
            "Controlled human malaria infection QTOF breath (.D + xml). "
            "Requires MassHunter peak table export before PatientVOCMatrix ingest."
        ),
    },
    "JID_2024_MALAWI": {
        "role": "external_reproducibility",
        "doi": "10.1093/infdis/jiae323",
        "status": "no_open_intensity_matrix",
        "note": "Independent Blantyre cohort; deposit intensity table when released.",
    },
}

_CACHE = DATA_DIR / "datasources" / "malaria_external"


def write_external_catalog(out_dir: Path | None = None) -> Path:
    out_dir = Path(out_dir or _CACHE)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "EXTERNAL_MALARIA_CATALOG.json"
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": EXTERNAL_CATALOG,
        "guidance": (
            "JBR breath-ML learning curves flatten near n≈50; single n=35 studies "
            "keep wide AUROC CIs. Prefer multi-cohort locked evaluation when a second "
            "intensity table becomes available."
        ),
    }
    path.write_text(json.dumps(payload, indent=2))
    return path


def load_external_patient_csv(
    matrix_csv: Path,
    labels_csv: Path,
    *,
    study_id: str = "EXTERNAL_MALARIA",
    disease_id: str = "malaria",
) -> PatientVOCMatrix:
    """Load a user-exported patient×VOC CSV + labels (0/1) as external cohort."""
    mat = pd.read_csv(matrix_csv, index_col=0)
    lab = pd.read_csv(labels_csv, index_col=0).iloc[:, 0].astype(int)
    mat.index = mat.index.astype(str)
    lab.index = lab.index.astype(str)
    common = [i for i in mat.index if i in set(lab.index)]
    mat = mat.loc[common].astype(float)
    lab = lab.loc[common]
    samples = [
        SampleRecord(
            sample_id=sid,
            subject_id=sid,
            study_id=study_id,
            label="disease" if int(lab.loc[sid]) == 1 else "control",
            disease_id=disease_id,
            modality="gcms_external",
        )
        for sid in common
    ]
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=disease_id,
        disease_name="Malaria",
        modality="gcms_external",
        unit="peak_intensity",
        matrix=mat,
        labels=lab.astype(int),
        samples=samples,
        metadata={
            "source": "user_external_csv",
            "n_voc_features": int(mat.shape[1]),
            "n_subjects": int(mat.shape[0]),
            "n_positive": int((lab == 1).sum()),
            "n_negative": int((lab == 0).sum()),
        },
    )


def combine_cohorts(
    matrices: list[PatientVOCMatrix],
    *,
    study_id: str = "MALARIA_POOLED",
) -> PatientVOCMatrix:
    """Inner-join VOC columns and concatenate subjects (prefix study ids)."""
    if not matrices:
        raise ValueError("No matrices to combine")
    # intersection of features
    cols = set(matrices[0].matrix.columns.astype(str))
    for m in matrices[1:]:
        cols &= set(m.matrix.columns.astype(str))
    cols = sorted(cols)
    if len(cols) < 2:
        raise ValueError("Fewer than 2 shared VOC features across cohorts")
    blocks = []
    labels = []
    samples = []
    for m in matrices:
        idx = [f"{m.study_id}:{sid}" for sid in m.subject_ids]
        block = m.matrix[cols].copy()
        block.index = idx
        blocks.append(block)
        lab = m.labels.copy()
        lab.index = idx
        labels.append(lab)
        for s in m.samples:
            samples.append(
                s.model_copy(
                    update={
                        "sample_id": f"{m.study_id}:{s.sample_id}",
                        "subject_id": f"{m.study_id}:{s.subject_id}",
                        "study_id": study_id,
                        "metadata": {**(s.metadata or {}), "source_study": m.study_id},
                    }
                )
            )
    X = pd.concat(blocks, axis=0)
    y = pd.concat(labels, axis=0).astype(int)
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=matrices[0].disease_id,
        disease_name=matrices[0].disease_name,
        modality="gcms_pooled",
        unit="peak_intensity",
        matrix=X.astype(float),
        labels=y,
        samples=samples,
        metadata={
            "source_studies": [m.study_id for m in matrices],
            "n_voc_features": len(cols),
            "n_subjects": int(X.shape[0]),
            "n_positive": int((y == 1).sum()),
            "n_negative": int((y == 0).sum()),
            "shared_vocs": cols,
        },
    )


__all__ = [
    "EXTERNAL_CATALOG",
    "combine_cohorts",
    "load_external_patient_csv",
    "write_external_catalog",
]
