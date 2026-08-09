"""BreathVOC interchange — JSON Schema + matrix bundle export/validate.

Productizes ``SampleRecord`` / ``SplitManifest`` / patient×VOC matrices as a
versioned open format labs can adopt alongside MetaboLights / OMNI tables.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from .patient_matrix import PatientVOCMatrix
from .schema import SampleRecord, SplitManifest

BREATHVOC_SCHEMA_VERSION = "BreathVOC-1.1"

BREATHVOC_JSON_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "$id": "https://exhalepath.dev/schemas/breathvoc-1.1.json",
    "title": "BreathVOC patient matrix bundle",
    "type": "object",
    "required": [
        "schema_version",
        "study_id",
        "disease_id",
        "matrix",
        "labels",
        "samples",
    ],
    "properties": {
        "schema_version": {"const": BREATHVOC_SCHEMA_VERSION},
        "study_id": {"type": "string"},
        "disease_id": {"type": "string"},
        "disease_name": {"type": "string"},
        "modality": {"type": "string"},
        "unit": {"type": "string"},
        "generated_utc": {"type": "string"},
        "matrix": {
            "type": "object",
            "description": "subject_id → {voc_id → intensity}",
        },
        "labels": {
            "type": "object",
            "description": "subject_id → 0|1",
        },
        "samples": {"type": "array", "items": {"type": "object"}},
        "split_manifest": {"type": ["object", "null"]},
        "metadata": {"type": "object"},
    },
}


def matrix_to_bundle(
    matrix: PatientVOCMatrix,
    *,
    split: SplitManifest | None = None,
) -> dict[str, Any]:
    return {
        "schema_version": BREATHVOC_SCHEMA_VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": matrix.study_id,
        "disease_id": matrix.disease_id,
        "disease_name": matrix.disease_name,
        "modality": matrix.modality,
        "unit": matrix.unit,
        "matrix": {
            str(i): {str(c): float(matrix.matrix.loc[i, c]) for c in matrix.matrix.columns}
            for i in matrix.matrix.index
        },
        "labels": {str(i): int(matrix.labels.loc[i]) for i in matrix.labels.index},
        "samples": [s.model_dump() for s in matrix.samples],
        "split_manifest": split.model_dump() if split is not None else None,
        "metadata": matrix.metadata,
    }


def bundle_to_matrix(bundle: dict[str, Any]) -> PatientVOCMatrix:
    validate_breathvoc_bundle(bundle)
    mat = pd.DataFrame(bundle["matrix"]).T.astype(float)
    lab = pd.Series(bundle["labels"], dtype=int)
    lab = lab.loc[mat.index]
    samples = [SampleRecord.model_validate(s) for s in bundle.get("samples") or []]
    return PatientVOCMatrix(
        study_id=bundle["study_id"],
        disease_id=bundle["disease_id"],
        disease_name=bundle.get("disease_name") or bundle["disease_id"],
        modality=bundle.get("modality") or "gcms",
        unit=bundle.get("unit") or "intensity",
        matrix=mat,
        labels=lab,
        samples=samples,
        metadata=dict(bundle.get("metadata") or {}),
    )


def validate_breathvoc_bundle(bundle: dict[str, Any]) -> list[str]:
    """Lightweight validator (no jsonschema dependency). Returns error list."""
    errors: list[str] = []
    if bundle.get("schema_version") != BREATHVOC_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {BREATHVOC_SCHEMA_VERSION}, "
            f"got {bundle.get('schema_version')}"
        )
    for key in ("study_id", "disease_id", "matrix", "labels", "samples"):
        if key not in bundle:
            errors.append(f"missing required key: {key}")
    if errors:
        return errors
    mat = bundle["matrix"]
    lab = bundle["labels"]
    if not isinstance(mat, dict) or not mat:
        errors.append("matrix must be a non-empty object")
    if not isinstance(lab, dict) or not lab:
        errors.append("labels must be a non-empty object")
    if isinstance(mat, dict) and isinstance(lab, dict):
        missing = [k for k in mat if k not in lab]
        if missing[:5]:
            errors.append(f"labels missing for matrix subjects e.g. {missing[:3]}")
    for s in bundle.get("samples") or []:
        try:
            SampleRecord.model_validate(s)
        except Exception as exc:  # noqa: BLE001
            errors.append(f"sample invalid: {exc}")
            break
    if bundle.get("split_manifest"):
        try:
            SplitManifest.model_validate(bundle["split_manifest"])
        except Exception as exc:  # noqa: BLE001
            errors.append(f"split_manifest invalid: {exc}")
    return errors


def export_breathvoc(
    matrix: PatientVOCMatrix,
    out_path: Path,
    *,
    split: SplitManifest | None = None,
) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bundle = matrix_to_bundle(matrix, split=split)
    out_path.write_text(json.dumps(bundle, indent=2))
    schema_path = out_path.with_name("breathvoc.schema.json")
    schema_path.write_text(json.dumps(BREATHVOC_JSON_SCHEMA, indent=2))
    return out_path


def load_and_validate_breathvoc(path: Path) -> tuple[PatientVOCMatrix, list[str]]:
    bundle = json.loads(Path(path).read_text())
    errors = validate_breathvoc_bundle(bundle)
    if errors:
        raise ValueError("BreathVOC validation failed: " + "; ".join(errors))
    return bundle_to_matrix(bundle), errors


__all__ = [
    "BREATHVOC_JSON_SCHEMA",
    "BREATHVOC_SCHEMA_VERSION",
    "bundle_to_matrix",
    "export_breathvoc",
    "load_and_validate_breathvoc",
    "matrix_to_bundle",
    "validate_breathvoc_bundle",
]
