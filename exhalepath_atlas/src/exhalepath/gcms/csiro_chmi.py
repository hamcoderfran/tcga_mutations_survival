"""CSIRO CHMI malaria breath (csiro:33843) — labels ready; intensity pending peak export.

Public deposit (DOI 10.25919/5b5b7530a39f4) ships Agilent ``.D`` + large ``mzdata.xml``,
not a patient×VOC intensity CSV. This module bundles the overview-derived
baseline-vs-peak-parasitemia labels so a MassHunter (or equivalent) peak table
can be joined for leave-one-study-out evaluation against ST000883.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from ..config import DATA_DIR
from .external_malaria import load_external_patient_csv
from .patient_matrix import PatientVOCMatrix

CSIRO_DIR = DATA_DIR / "datasources" / "malaria_external"
CSIRO_LABELS = CSIRO_DIR / "csiro_chmi_labels.csv"
CSIRO_SAMPLES = CSIRO_DIR / "csiro_chmi_samples.csv"
CSIRO_MANIFEST = CSIRO_DIR / "csiro_chmi_manifest.json"
CSIRO_OVERVIEW = CSIRO_DIR / "csiro_33843_overview.xlsx"

STUDY_ID = "CSIRO_CHMI_33843"


def load_csiro_manifest() -> dict[str, Any]:
    if not CSIRO_MANIFEST.exists():
        raise FileNotFoundError(CSIRO_MANIFEST)
    return json.loads(CSIRO_MANIFEST.read_text())


def load_csiro_labels() -> pd.DataFrame:
    """Baseline (0) vs peak-parasitemia (1) sample labels (n=14 for 7 subjects)."""
    if not CSIRO_LABELS.exists():
        raise FileNotFoundError(CSIRO_LABELS)
    df = pd.read_csv(CSIRO_LABELS)
    df["sample_id"] = df["sample_id"].astype(str)
    df["label"] = df["label"].astype(int)
    return df


def csiro_readiness() -> dict[str, Any]:
    """Status block for reports: labels bundled; intensity requires offline peak table."""
    meta = load_csiro_manifest() if CSIRO_MANIFEST.exists() else {}
    labels_ok = CSIRO_LABELS.exists()
    return {
        "study_id": STUDY_ID,
        "labels_bundled": labels_ok,
        "overview_bundled": CSIRO_OVERVIEW.exists(),
        "intensity_status": meta.get("intensity_status", "unknown"),
        "n_subjects": meta.get("n_subjects"),
        "n_labeled_baseline_peak": meta.get("n_labeled_baseline_peak"),
        "doi": meta.get("doi"),
        "url": meta.get("url"),
        "labels_path": str(CSIRO_LABELS) if labels_ok else None,
        "note": meta.get("intensity_note"),
        "loso_ready_when": (
            "Provide --external-matrix patient×VOC CSV indexed by Sample name "
            f"(see {CSIRO_LABELS.name}) with --loso"
        ),
    }


def load_csiro_patient_matrix(
    matrix_csv: Path,
    *,
    labels_csv: Path | None = None,
) -> PatientVOCMatrix:
    """Join a user-exported CSIRO peak table to bundled (or custom) CHMI labels."""
    labels = Path(labels_csv) if labels_csv else CSIRO_LABELS
    # load_external expects labels CSV with index col = sample id and a 0/1 column
    lab = pd.read_csv(labels)
    if "sample_id" in lab.columns:
        lab = lab.set_index("sample_id")
    else:
        lab = lab.set_index(lab.columns[0])
    if "label" not in lab.columns:
        lab = lab.iloc[:, [0]]
        lab.columns = ["label"]
    else:
        lab = lab[["label"]]
    tmp = Path(matrix_csv).with_suffix(".csiro_labels_tmp.csv")
    # write a minimal labels file next to matrix for the shared loader
    # Prefer in-memory path via temp in same dir
    import tempfile

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".csv", delete=False, dir=str(Path(matrix_csv).parent)
    ) as fh:
        lab.to_csv(fh.name)
        tmp = Path(fh.name)
    try:
        return load_external_patient_csv(
            Path(matrix_csv), tmp, study_id=STUDY_ID, disease_id="malaria"
        )
    finally:
        tmp.unlink(missing_ok=True)


def export_csiro_label_bundle(out_dir: Path) -> dict[str, str]:
    """Copy bundled CSIRO label artifacts into an output directory."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    for src in (CSIRO_LABELS, CSIRO_SAMPLES, CSIRO_MANIFEST):
        if src.exists():
            dest = out_dir / src.name
            dest.write_text(src.read_text())
            paths[src.stem] = str(dest)
    ready = csiro_readiness()
    (out_dir / "CSIRO_CHMI_READINESS.json").write_text(json.dumps(ready, indent=2))
    paths["readiness"] = str(out_dir / "CSIRO_CHMI_READINESS.json")
    return paths


__all__ = [
    "CSIRO_DIR",
    "CSIRO_LABELS",
    "STUDY_ID",
    "csiro_readiness",
    "export_csiro_label_bundle",
    "load_csiro_labels",
    "load_csiro_manifest",
    "load_csiro_patient_matrix",
]
