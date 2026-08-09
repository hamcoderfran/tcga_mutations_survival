"""OEM library surface for partner LIMS embeds (recurring ACV).

Minimal stable API:
- ``score_sample`` — disease signature match on an observed VOC dict
- ``verify_locked_split`` — preregistration hash check against a matrix
- ``import_feature_table`` — BreathVOC JSON / OMNI CSV → PatientVOCMatrix
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal, Optional

from ..gcms.claim_ledger import build_claim_ledger
from ..gcms.interchange import load_and_validate_breathvoc
from ..gcms.locked_split import load_split_manifest, verify_split_manifest
from ..gcms.omni_partner import load_omni_style_csv
from ..gcms.patient_matrix import PatientVOCMatrix, load_mw_patient_matrix
from ..gcms.score import disease_signature, score_observed_vector


def cite_vocs_for_disease(
    disease_id: str,
    voc_ids: list[str],
) -> dict[str, dict[str, Any]]:
    """Map voc_id → evidence grade + DOI from the claim ledger."""
    ledger = build_claim_ledger(disease_ids=[disease_id], include_atlas_priors=True)
    out: dict[str, dict[str, Any]] = {}
    for c in ledger.get("claims") or []:
        if c.get("disease_id") != disease_id:
            continue
        vid = c.get("voc_id")
        if vid not in voc_ids:
            continue
        # Prefer quantified over atlas_prior when both exist
        prev = out.get(vid)
        grade = str(c.get("evidence_grade") or "unknown")
        rank = {
            "quantified": 0,
            "quantified_null": 1,
            "directional_only": 2,
            "mixed": 3,
            "mechanism_narrative": 4,
            "atlas_prior": 5,
        }.get(grade, 9)
        if prev is None or rank < prev.get("_rank", 99):
            out[vid] = {
                "voc_id": vid,
                "evidence_grade": grade,
                "doi": c.get("doi"),
                "log2fc_claim": c.get("log2fc"),
                "direction": c.get("direction"),
                "source": c.get("source"),
                "circularity_risk": c.get("circularity_risk"),
                "notes": c.get("notes"),
                "_rank": rank,
            }
    for v in out.values():
        v.pop("_rank", None)
    for vid in voc_ids:
        out.setdefault(
            vid,
            {
                "voc_id": vid,
                "evidence_grade": "uncited",
                "doi": None,
                "log2fc_claim": None,
                "direction": None,
                "source": None,
                "circularity_risk": True,
                "notes": "No ledger claim — do not present as literature-backed",
            },
        )
    return out


def score_sample(
    disease_id: str,
    vocs: dict[str, float],
    *,
    signature_source: Literal["hybrid", "stack", "literature"] = "hybrid",
    method: str = "cosine",
    top_n: int = 40,
) -> dict[str, Any]:
    """Score an observed VOC vector; attach ledger citations for overlapping VOCs."""
    sig = disease_signature(disease_id, source=signature_source, top_n=top_n)
    scored = score_observed_vector(vocs, sig, method=method)  # type: ignore[arg-type]
    used = list(scored.get("vocs_used") or [])
    cites = cite_vocs_for_disease(disease_id, used or list(vocs.keys())[:20])
    # rank signature VOCs by |log2fc| with cites
    ranked = sorted(sig.items(), key=lambda kv: -abs(kv[1]))[:20]
    top_vocs = []
    for vid, fc in ranked:
        cite = cites.get(vid) or cite_vocs_for_disease(disease_id, [vid]).get(vid)
        top_vocs.append(
            {
                "voc_id": vid,
                "signature_log2fc": float(fc),
                "observed": vocs.get(vid),
                "evidence_grade": (cite or {}).get("evidence_grade"),
                "doi": (cite or {}).get("doi"),
                "circularity_risk": (cite or {}).get("circularity_risk"),
            }
        )
    return {
        "disease_id": disease_id,
        "signature_source": signature_source,
        "score": scored.get("score"),
        "n_overlap": scored.get("n_overlap"),
        "method": method,
        "top_vocs": top_vocs,
        "citations": cites,
        "positioning": (
            "Cut the cost of wrong VOC panels — ledger-cited hypotheses, "
            "not clinical diagnostic claims."
        ),
    }


def verify_locked_split(
    manifest_path: Path,
    *,
    study_id: str | None = None,
    matrix: PatientVOCMatrix | None = None,
) -> dict[str, Any]:
    """Verify a locked split manifest against a bundled or provided matrix."""
    manifest = load_split_manifest(Path(manifest_path))
    if matrix is None:
        sid = study_id or manifest.dataset_id
        matrix = load_mw_patient_matrix(sid)
    problems = verify_split_manifest(manifest, matrix)
    return {
        "ok": len(problems) == 0,
        "problems": problems,
        "dataset_id": manifest.dataset_id,
        "content_sha256": manifest.content_sha256,
        "n_splits": manifest.n_splits,
        "strategy": manifest.strategy,
    }


def import_feature_table(
    path: Path,
    *,
    fmt: Literal["auto", "breathvoc", "omni"] = "auto",
    disease_id: str = "malaria",
    study_id: str = "PARTNER_IMPORT",
) -> dict[str, Any]:
    """One-click import of partner BreathVOC JSON or OMNI-style CSV."""
    path = Path(path)
    kind = fmt
    if kind == "auto":
        if path.suffix.lower() == ".json" or "breathvoc" in path.name.lower():
            kind = "breathvoc"
        else:
            kind = "omni"
    if kind == "breathvoc":
        matrix, errors = load_and_validate_breathvoc(path)
        return {
            "format": "breathvoc",
            "study_id": matrix.study_id,
            "disease_id": matrix.disease_id,
            "n_subjects": int(matrix.matrix.shape[0]),
            "n_vocs": int(matrix.matrix.shape[1]),
            "validation_errors": errors,
            "matrix": matrix,
        }
    matrix = load_omni_style_csv(
        path, study_id=study_id, disease_id=disease_id
    )
    return {
        "format": "omni",
        "study_id": matrix.study_id,
        "disease_id": matrix.disease_id,
        "n_subjects": int(matrix.matrix.shape[0]),
        "n_vocs": int(matrix.matrix.shape[1]),
        "validation_errors": [],
        "matrix": matrix,
    }


def kit_info() -> dict[str, Any]:
    return {
        "name": "exhalepath.oem",
        "surface": ["score_sample", "verify_locked_split", "import_feature_table"],
        "positioning": "Cut the cost of wrong VOC panels — embeddable diligence, not IVD",
        "not_a_medical_device": True,
    }


__all__ = [
    "cite_vocs_for_disease",
    "import_feature_table",
    "kit_info",
    "score_sample",
    "verify_locked_split",
]
