"""Partner intensity diligence pipeline (Pillar 2).

OMNI/BreathVOC CSV → import → lock-split → signature score summary.

Honesty: plumbing for partner diligence. Bundled fixtures may be derived from
public ST000883 splits — not multi-site prospective clinical validation.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import pandas as pd

from ..config import DATA_DIR, PACKAGE_ROOT
from ..gcms.interchange import export_breathvoc
from ..gcms.locked_split import lock_split
from ..gcms.score import disease_signature, score_observed_vector
from ..knowledge.voc_alias_registry import audit_column_names
from ..oem.kit import import_feature_table


def _default_partner_csv() -> Path:
    base = DATA_DIR / "datasources" / "partner_loso"
    path = base / "PARTNER_SITE_A_omni_feature_table.csv"
    if not path.exists():
        path = PACKAGE_ROOT / "data" / "datasources" / "partner_loso" / "PARTNER_SITE_A_omni_feature_table.csv"
    return path


def run_partner_diligence_pipeline(
    *,
    path: Path | None = None,
    fmt: str = "auto",
    disease_id: str = "malaria",
    out_dir: Path = Path("runs/partner_diligence"),
    seed: int = 42,
    n_splits: int = 5,
    run_diagnostic: bool = True,
) -> dict[str, Any]:
    """Import partner table, lock a patient split, optionally score signatures."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if path is None:
        path = _default_partner_csv()
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Partner feature table not found: {path}")

    packed = import_feature_table(path, fmt=fmt, disease_id=disease_id)  # type: ignore[arg-type]
    matrix = packed["matrix"]
    study_id = matrix.study_id or path.stem

    alias_audit: dict[str, Any] | None = None
    try:
        if path.suffix.lower() in {".csv", ".tsv"}:
            raw_cols = list(pd.read_csv(path, nrows=0).columns.astype(str))
            alias_audit = audit_column_names(raw_cols)
    except Exception as exc:  # noqa: BLE001
        alias_audit = {"ok": False, "error": str(exc)}

    export_breathvoc(matrix, out_dir / f"{study_id}.breathvoc.json")
    matrix.matrix.to_csv(out_dir / f"{study_id}_matrix.csv")
    matrix.labels.to_frame("label").to_csv(out_dir / f"{study_id}_labels.csv")

    split_path = out_dir / f"{study_id}_split_v1.json"
    manifest = lock_split(
        matrix,
        strategy="stratified_kfold",
        n_splits=n_splits,
        seed=seed,
        out_path=split_path,
    )

    diagnostic: dict[str, Any] | None = None
    if run_diagnostic:
        try:
            sig = disease_signature(disease_id, source="hybrid")
            scores: list[float] = []
            y: list[int] = []
            for sid, row in matrix.matrix.iterrows():
                vocs = {str(c): float(row[c]) for c in matrix.matrix.columns if row[c] == row[c]}
                sc = score_observed_vector(vocs, sig, method="cosine")
                if sc.get("score") is None:
                    continue
                scores.append(float(sc["score"]))
                y.append(int(matrix.labels.loc[sid]))
            pos = [s for s, yy in zip(scores, y) if yy == 1]
            neg = [s for s, yy in zip(scores, y) if yy == 0]
            diagnostic = {
                "mode": "partner_signature_score_summary",
                "n_scored": len(y),
                "n_positive": int(sum(y)),
                "n_negative": int(len(y) - sum(y)),
                "score_mean_pos": float(np.mean(pos)) if pos else None,
                "score_mean_neg": float(np.mean(neg)) if neg else None,
                "note": (
                    "Signature score separation summary for partner matrices. "
                    "Full nested AUROC uses voc eval-patient-diagnostic on registered "
                    "study adapters. Not a clinical claim."
                ),
            }
        except Exception as exc:  # noqa: BLE001
            diagnostic = {"ok": False, "error": str(exc)}

    meta = {k: v for k, v in packed.items() if k != "matrix"}
    unmapped_rate = (alias_audit or {}).get("unmapped_rate")
    report = {
        "schema_version": "PartnerDiligence-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Partner intensity diligence pipeline (import → lock-split → score summary). "
            "Bundled OMNI fixtures may be derived from public ST000883 — not multi-site "
            "prospective clinical validation."
        ),
        "input_path": str(path),
        "format": meta.get("format"),
        "disease_id": disease_id,
        "study_id": study_id,
        "n_subjects": meta.get("n_subjects") or int(matrix.matrix.shape[0]),
        "n_vocs": meta.get("n_vocs") or int(matrix.matrix.shape[1]),
        "unmapped_rate": unmapped_rate,
        "alias_audit": alias_audit,
        "split_manifest": str(split_path),
        "split_sha256": manifest.content_sha256,
        "import_meta": meta,
        "diagnostic": diagnostic,
        "clinical_claim": False,
    }
    (out_dir / "PARTNER_DILIGENCE.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    sha = report.get("split_sha256") or ""
    lines = [
        "# Partner intensity diligence",
        "",
        f"Generated: {report['generated_utc']}",
        "",
        f"> {report['honesty']}",
        "",
        f"- input: `{path}`",
        f"- subjects×VOCs: {report['n_subjects']}×{report['n_vocs']}",
        f"- split: `{split_path.name}` sha256=`{sha[:16]}…`",
        f"- unmapped_rate: {report.get('unmapped_rate')}",
        "",
        "Next: `voc eval-patient-diagnostic` on bundled public studies for nested AUROC,",
        "or `voc panel-decision --disease <id>` for a go/no-go memo.",
        "",
    ]
    (out_dir / "PARTNER_DILIGENCE.md").write_text("\n".join(lines))
    return report


__all__ = ["run_partner_diligence_pipeline"]
