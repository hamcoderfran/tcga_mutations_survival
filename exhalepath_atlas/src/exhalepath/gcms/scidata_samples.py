"""Per-sample Sci Data 2024 peak-table adapters (not cohort means).

Peak CSVs are **feature × sample** (rows = compounds with pubchem_CID /
IUPAC Name; columns = sample ids). This module melts/transposes to
patient × VOC matrices and joins age/sex/spirometry from
CBD_metadata_for_ver3.xlsx.

Honest limitations
------------------
- No healthy control arm → AUROC is one-vs-rest across pulmonary cohorts.
- Metadata has age + sex + spirometry; **smoking is absent**.
- Blank columns are not present; blank-ratio filters apply detection-fraction
  only unless blanks are supplied separately.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..ingest.public_breath import _load_peak_table, map_compound_to_voc
from ..knowledge.loader import KnowledgeBase
from .patient_matrix import PatientVOCMatrix
from .schema import SampleRecord

_PKG = Path(__file__).resolve().parents[1]
_DATA = _PKG / "data" / "public_breath"

COHORT_FILES = {
    "asthma": "Asthma_peaktable_ver3.csv",
    "copd": "COPD_peaktable_ver3.csv",
    "bronchiectasis": "Bronchi_peaktable_ver3.csv",
}

METADATA_FILE = "CBD_metadata_for_ver3.xlsx"

_SHEET_TO_COHORT = {
    "Asthma": "asthma",
    "COPD": "copd",
    "Bronchiectasis": "bronchiectasis",
}

COHORT_TO_DISEASE = {
    "asthma": "asthma",
    "copd": "copd",
    # atlas has no dedicated bronchiectasis prior — COPD is the documented proxy
    "bronchiectasis": "copd",
}


def _norm_sample_id(x: Any) -> str:
    s = str(x).strip()
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    if s.isdigit():
        return str(int(s))
    return s


def _load_metadata_xlsx(path: Path) -> Dict[str, Dict[str, Any]]:
    """Map cohort:sample_id → {age, sex, fev1_pp, fvc_pp, ...}."""
    out: Dict[str, Dict[str, Any]] = {}
    if not path.is_file():
        return out
    try:
        xl = pd.ExcelFile(path)
    except Exception:  # noqa: BLE001
        return out

    for sheet, cohort in _SHEET_TO_COHORT.items():
        if sheet not in xl.sheet_names:
            continue
        df = xl.parse(sheet)
        cols = {str(c).strip().lower(): c for c in df.columns}
        id_col = cols.get("id") or cols.get("sample") or df.columns[0]
        age_col = cols.get("age")
        sex_col = cols.get("sex")
        fev1_col = cols.get("fev10 pp") or cols.get("fev1 pp") or cols.get("fev1")
        fvc_col = cols.get("fvc pp") or cols.get("fvc")
        bmi_col = cols.get("bmi")

        for _, row in df.iterrows():
            sid = _norm_sample_id(row[id_col])
            key = f"{cohort}:{sid}"
            meta: Dict[str, Any] = {"cohort": cohort, "sample_id": sid}
            if age_col is not None and pd.notna(row.get(age_col)):
                try:
                    meta["age"] = float(row[age_col])
                except Exception:  # noqa: BLE001
                    pass
            if sex_col is not None and pd.notna(row.get(sex_col)):
                sex = str(row[sex_col]).strip().lower()
                if sex.startswith("m"):
                    meta["sex"] = "male"
                elif sex.startswith("f"):
                    meta["sex"] = "female"
                else:
                    meta["sex"] = sex
            for name, col in (("fev1_pp", fev1_col), ("fvc_pp", fvc_col), ("bmi", bmi_col)):
                if col is not None and pd.notna(row.get(col)):
                    try:
                        meta[name] = float(row[col])
                    except Exception:  # noqa: BLE001
                        pass
            # smoking absent in Sci Data 2024 deposit
            out[key] = meta
    return out


def load_scidata_long(
    data_dir: Optional[Path] = None,
    cohorts: Optional[Sequence[str]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (long peak table, covariates DataFrame).

    Long columns: cohort, sample_id, uid, pubchem_cid, iupac_name, intensity, voc_id
    """
    root = Path(data_dir) if data_dir else _DATA
    use = list(cohorts) if cohorts else list(COHORT_FILES.keys())
    meta_map = _load_metadata_xlsx(root / METADATA_FILE)
    kb = KnowledgeBase()

    frames: List[pd.DataFrame] = []
    cov_rows: List[Dict[str, Any]] = []
    seen_uid: set[str] = set()

    for cohort in use:
        fname = COHORT_FILES.get(cohort)
        if not fname:
            raise ValueError(f"Unknown cohort {cohort}; known={list(COHORT_FILES)}")
        path = root / fname
        if not path.is_file():
            raise FileNotFoundError(path)
        long = _load_peak_table(path)
        long["sample_id"] = long["sample_id"].map(_norm_sample_id)
        long["cohort"] = cohort
        long["uid"] = long["cohort"] + "_" + long["sample_id"].astype(str)
        long["voc_id"] = [
            map_compound_to_voc(pubchem_cid=r.pubchem_cid, iupac_name=r.iupac_name, kb=kb)
            for r in long.itertuples()
        ]
        frames.append(long)

        for sid in sorted(long["sample_id"].unique()):
            uid = f"{cohort}_{sid}"
            if uid in seen_uid:
                continue
            seen_uid.add(uid)
            m = meta_map.get(f"{cohort}:{sid}", {})
            cov_rows.append(
                {
                    "sample_id": uid,
                    "cohort": cohort,
                    "raw_sample_id": sid,
                    "age": m.get("age"),
                    "sex": m.get("sex"),
                    "fev1_pp": m.get("fev1_pp"),
                    "fvc_pp": m.get("fvc_pp"),
                    "bmi": m.get("bmi"),
                    "smoking": m.get("smoking"),  # expected None for Sci Data 2024
                }
            )

    peaks = pd.concat(frames, ignore_index=True)
    covariates = pd.DataFrame(cov_rows)
    return peaks, covariates


def load_scidata_intensity_matrix(
    data_dir: Optional[Path] = None,
    cohorts: Optional[Sequence[str]] = None,
    *,
    mapped_vocs_only: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Sample×feature intensity matrix (+ covariates + cohort labels).

    Feature columns are atlas ``voc_id`` when mapped; otherwise
    ``cid:<pubchem>`` so blank/detection filters can still run on the full table.
    """
    peaks, covariates = load_scidata_long(data_dir=data_dir, cohorts=cohorts)
    if mapped_vocs_only:
        peaks = peaks.dropna(subset=["voc_id"]).copy()
        peaks["feature"] = peaks["voc_id"]
    else:
        peaks = peaks.copy()
        peaks["feature"] = peaks["voc_id"].where(
            peaks["voc_id"].notna(),
            "cid:" + peaks["pubchem_cid"].astype(str),
        )

    piv = peaks.pivot_table(
        index="uid", columns="feature", values="intensity", aggfunc="median"
    )
    piv = piv.fillna(0.0)
    piv.index = piv.index.astype(str)
    piv.index.name = "sample_id"
    cov = covariates.set_index("sample_id").reindex(piv.index)
    cov.index.name = "sample_id"
    labels = cov["cohort"].copy()
    labels.index = labels.index.astype(str)
    return piv.astype(float), cov.reset_index(), labels


def load_scidata_ovr_matrix(
    positive_cohort: str,
    *,
    data_dir: Optional[Path] = None,
    cohorts: Optional[Sequence[str]] = None,
    mapped_vocs_only: bool = True,
) -> PatientVOCMatrix:
    """One-vs-rest PatientVOCMatrix for a Sci Data pulmonary cohort."""
    if positive_cohort not in COHORT_FILES:
        raise ValueError(
            f"positive_cohort must be one of {list(COHORT_FILES)}; got {positive_cohort}"
        )
    intensity, cov, cohort_labels = load_scidata_intensity_matrix(
        data_dir=data_dir, cohorts=cohorts, mapped_vocs_only=mapped_vocs_only
    )
    y = (cohort_labels == positive_cohort).astype(int)
    disease_id = COHORT_TO_DISEASE[positive_cohort]
    samples: List[SampleRecord] = []
    for sid in intensity.index.astype(str):
        row = cov.loc[cov["sample_id"] == sid]
        age = sex = smoking = None
        meta_extra: Dict[str, Any] = {"cohort": str(cohort_labels.loc[sid])}
        if len(row):
            r = row.iloc[0]
            age = None if pd.isna(r.get("age")) else float(r["age"])
            sex = None if pd.isna(r.get("sex")) else str(r["sex"])
            smoking = None if pd.isna(r.get("smoking")) else str(r["smoking"])
            for k in ("fev1_pp", "fvc_pp", "bmi", "raw_sample_id"):
                if k in r and pd.notna(r[k]):
                    meta_extra[k] = r[k]
        samples.append(
            SampleRecord(
                sample_id=sid,
                subject_id=sid,
                study_id=f"SCIDATA2024_{positive_cohort}_ovr",
                label="disease" if int(y.loc[sid]) == 1 else "control",
                disease_id=disease_id,
                age=age,
                sex=sex,
                smoking_status=smoking,
                modality="gcms_scidata",
                metadata=meta_extra,
            )
        )

    return PatientVOCMatrix(
        study_id=f"SCIDATA2024_{positive_cohort}_ovr",
        disease_id=disease_id,
        disease_name=positive_cohort.replace("_", " ").title(),
        modality="gcms_scidata",
        unit="peak_area",
        matrix=intensity.astype(float),
        labels=y.astype(int),
        samples=samples,
        metadata={
            "doi": "10.1038/s41597-024-03216-0",
            "title": f"Sci Data 2024 per-sample peak table · {positive_cohort} one-vs-rest",
            "positive_cohort": positive_cohort,
            "n_voc_features": int(intensity.shape[1]),
            "n_subjects": int(intensity.shape[0]),
            "n_positive": int((y == 1).sum()),
            "n_negative": int((y == 0).sum()),
            "label_scheme": "one_vs_rest_pulmonary_cohorts",
            "smoking_available": False,
            "mapped_vocs_only": mapped_vocs_only,
            "caveat": (
                "No healthy controls — negatives are other pulmonary cohorts. "
                "Not disease-vs-healthy AUROC."
            ),
        },
    )


def list_scidata_cohorts(data_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    root = Path(data_dir) if data_dir else _DATA
    out = []
    for cohort, fname in COHORT_FILES.items():
        path = root / fname
        out.append(
            {
                "cohort": cohort,
                "file": fname,
                "available": path.is_file(),
                "disease_id_proxy": COHORT_TO_DISEASE[cohort],
            }
        )
    out.append(
        {
            "metadata": METADATA_FILE,
            "available": (root / METADATA_FILE).is_file(),
            "fields": ["age", "sex", "FEV1 PP", "FVC PP", "BMI"],
            "smoking": False,
        }
    )
    return out


__all__ = [
    "COHORT_FILES",
    "COHORT_TO_DISEASE",
    "list_scidata_cohorts",
    "load_scidata_intensity_matrix",
    "load_scidata_long",
    "load_scidata_ovr_matrix",
]
