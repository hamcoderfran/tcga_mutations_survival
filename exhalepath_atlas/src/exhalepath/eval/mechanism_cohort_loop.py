"""Mechanism → cohort closed loop: pack/prior → signature → nested patient AUROC.

The revolutionary claim is not a higher AUROC — it is a *feedback loop*:
mechanism hypotheses are scored on locked patient matrices with optimism gaps,
so packs can be ablated and improved against real intensities.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.metrics import roc_auc_score

from ..gcms.locked_split import lock_split
from ..gcms.patient_matrix import list_bundled_diagnostic_studies, load_mw_patient_matrix
from ..gcms.score import disease_signature, score_patients
from .patient_diagnostic import _fold_control_relative_scores


def _nested_signature_auroc(
    study_id: str,
    disease_id: str,
    *,
    source: str,
    n_splits: int = 5,
    seed: int = 42,
) -> dict[str, Any]:
    matrix = load_mw_patient_matrix(study_id).log1p()
    try:
        sig = disease_signature(
            disease_id,
            source=source,  # type: ignore[arg-type]
            top_n=40,
            voc_ids=list(matrix.matrix.columns.astype(str)),
        )
    except Exception as exc:  # noqa: BLE001
        return {"status": "error", "error": str(exc), "source": source}

    if len(sig) < 2:
        return {
            "status": "skipped",
            "reason": "signature_overlap_<2",
            "source": source,
            "n_sig": len(sig),
        }

    manifest = lock_split(
        matrix, strategy="stratified_kfold", n_splits=n_splits, seed=seed
    )
    fold_aucs = []
    for fold in manifest.folds:
        sc = _fold_control_relative_scores(
            matrix, sig, fold.train_subject_ids, fold.test_subject_ids, method="cosine"
        )
        y_te = matrix.labels.loc[fold.test_subject_ids].to_numpy(dtype=int)
        if len(np.unique(y_te)) < 2:
            continue
        fold_aucs.append(float(roc_auc_score(y_te, sc.to_numpy())))

    full = score_patients(matrix, sig, method="cosine", reference="control_mean")
    y = matrix.labels.loc[full.index].astype(int)
    non_nested = float(roc_auc_score(y, full)) if len(np.unique(y)) >= 2 else None
    nested = float(np.mean(fold_aucs)) if fold_aucs else None
    return {
        "status": "ok",
        "source": source,
        "n_signature_vocs": len(sig),
        "signature_vocs": sorted(sig.keys()),
        "nested_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": (non_nested - nested)
        if (non_nested is not None and nested is not None)
        else None,
        "n_folds_scored": len(fold_aucs),
        "split_sha256": manifest.content_sha256,
    }


def evaluate_mechanism_cohort_loop(
    *,
    out_dir: Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    """Run closed loop on all bundled diagnostic MW studies × signature sources."""
    studies = list_bundled_diagnostic_studies()
    rows = []
    for st in studies:
        if not st.get("available"):
            continue
        sid = st["study_id"]
        did = st["disease_id"]
        entry: dict[str, Any] = {
            "study_id": sid,
            "disease_id": did,
            "disease_name": st.get("disease_name"),
            "sources": {},
        }
        for source in ("hybrid", "stack", "literature"):
            try:
                entry["sources"][source] = _nested_signature_auroc(
                    sid, did, source=source, seed=seed
                )
            except Exception as exc:  # noqa: BLE001
                entry["sources"][source] = {"status": "error", "error": str(exc)}
        rows.append(entry)

    report = {
        "title": "Mechanism → cohort closed loop",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Nested AUROC tests whether mechanism/literature signatures transfer to "
            "patient intensities. Literature may be partly circular with atlas priors. "
            "Research enablement — not a diagnostic claim."
        ),
        "studies": rows,
        "how_to_improve_packs": [
            "Ablate VOCs with negative contribution on nested folds",
            "Prefer quantified literature claims over atlas priors (see claim ledger)",
            "Add external intensity cohort + LOSO before raising clinical language",
        ],
    }

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "MECHANISM_COHORT_LOOP.json").write_text(
            json.dumps(report, indent=2, default=str)
        )
        lines = [
            "# Mechanism → cohort closed loop",
            "",
            f"Generated: {report['generated_utc']}",
            "",
            f"> {report['honesty']}",
            "",
        ]
        for st in rows:
            lines.append(f"## {st['study_id']} ({st['disease_id']})")
            lines.append("")
            for src, block in (st.get("sources") or {}).items():
                lines.append(
                    f"- `{src}`: nested={_pct(block.get('nested_auroc'))} "
                    f"non-nested={_pct(block.get('non_nested_auroc'))} "
                    f"gap={_pct(block.get('optimism_gap'))} "
                    f"n_vocs={block.get('n_signature_vocs')}"
                )
            lines.append("")
        (out_dir / "MECHANISM_COHORT_LOOP.md").write_text("\n".join(lines))
    return report


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


__all__ = ["evaluate_mechanism_cohort_loop"]
