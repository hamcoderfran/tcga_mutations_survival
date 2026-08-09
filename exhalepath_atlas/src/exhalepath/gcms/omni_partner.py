"""OMNI-style / partner feature-table ingest + diligence LOSO fixture.

Owlstone OMNI publishes example feature tables (sample × compound). Partners
drop the same shape; we map columns to atlas ``voc_id`` and run LOSO.

The bundled diligence fixture is a *stratified split of ST000883* into two
OMNI-style CSVs — format + LOSO plumbing for buyer diligence slides.
It is **not** an independent clinical external validation cohort.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..config import DATA_DIR
from ..ingest.real_breath_corpus import _map_voc
from .external_malaria import leave_one_study_out, load_external_patient_csv
from .patient_matrix import PatientVOCMatrix, load_mw_patient_matrix
from .schema import SampleRecord

PARTNER_DIR = DATA_DIR / "datasources" / "partner_loso"
FIXTURE_MANIFEST = PARTNER_DIR / "PARTNER_LOSO_MANIFEST.json"


def _map_column(name: str) -> Optional[str]:
    n = str(name).strip()
    if n.lower() in {"sample_id", "sample", "subject_id", "subject", "label", "class", "group"}:
        return None
    # already atlas id?
    if _map_voc(n):
        return _map_voc(n)
    # snake / spaced
    cand = n.lower().replace(" ", "_").replace("-", "_")
    if _map_voc(cand):
        return _map_voc(cand)
    # try raw as voc_id if looks like atlas
    if all(c.isalnum() or c == "_" for c in cand) and len(cand) > 2:
        return cand
    return _map_voc(n)


def load_omni_style_csv(
    path: Path,
    *,
    study_id: str = "PARTNER_OMNI",
    disease_id: str = "malaria",
    label_col: str | None = None,
    sample_col: str | None = None,
) -> PatientVOCMatrix:
    """Load partner OMNI-style table: rows=samples, columns=compounds (+ label)."""
    df = pd.read_csv(path)
    cols = [str(c) for c in df.columns]
    # detect sample id
    if sample_col and sample_col in df.columns:
        sid_col = sample_col
    else:
        sid_col = None
        for c in cols:
            if c.lower() in {"sample_id", "sample", "subject_id", "subject", "id"}:
                sid_col = c
                break
        if sid_col is None:
            df = df.copy()
            df.insert(0, "sample_id", [f"S{i:03d}" for i in range(len(df))])
            sid_col = "sample_id"
    # detect label
    if label_col and label_col in df.columns:
        ycol = label_col
    else:
        ycol = None
        for c in cols:
            if c.lower() in {"label", "class", "group", "disease", "status"}:
                ycol = c
                break
    if ycol is None:
        raise ValueError("OMNI-style CSV needs a label/class/group column (0/1 or case/control)")

    sample_ids = df[sid_col].astype(str)
    raw_y = df[ycol]
    labels = []
    for v in raw_y:
        s = str(v).strip().lower()
        if s in {"1", "true", "case", "disease", "positive", "pos", "malaria"}:
            labels.append(1)
        elif s in {"0", "false", "control", "negative", "neg", "healthy"}:
            labels.append(0)
        else:
            try:
                labels.append(1 if float(v) >= 0.5 else 0)
            except Exception as exc:  # noqa: BLE001
                raise ValueError(f"Cannot parse label {v!r}") from exc

    feat_cols = [c for c in df.columns if c not in {sid_col, ycol}]
    mapped: dict[str, pd.Series] = {}
    for c in feat_cols:
        vid = _map_column(str(c))
        if not vid:
            continue
        series = pd.to_numeric(df[c], errors="coerce")
        if vid in mapped:
            mapped[vid] = mapped[vid].combine(series, lambda a, b: np.nanmean([a, b]))
        else:
            mapped[vid] = series
    if len(mapped) < 2:
        raise ValueError("Fewer than 2 mappable VOC columns in partner table")
    mat = pd.DataFrame(mapped)
    mat.index = sample_ids.to_numpy()
    mat = mat.apply(pd.to_numeric, errors="coerce")
    for col in mat.columns:
        med = float(mat[col].median(skipna=True))
        if np.isnan(med):
            med = 0.0
        mat[col] = mat[col].fillna(med)
    y = pd.Series(labels, index=mat.index, dtype=int, name="label")
    samples = [
        SampleRecord(
            sample_id=sid,
            subject_id=sid,
            study_id=study_id,
            label="disease" if int(y.loc[sid]) == 1 else "control",
            disease_id=disease_id,
            modality="omni_partner",
        )
        for sid in mat.index.astype(str)
    ]
    return PatientVOCMatrix(
        study_id=study_id,
        disease_id=disease_id,
        disease_name=disease_id.replace("_", " ").title(),
        modality="omni_partner",
        unit="peak_intensity",
        matrix=mat.astype(float),
        labels=y.astype(int),
        samples=samples,
        metadata={
            "source": "omni_style_csv",
            "path": str(path),
            "n_voc_features": int(mat.shape[1]),
            "n_subjects": int(mat.shape[0]),
            "n_positive": int((y == 1).sum()),
            "n_negative": int((y == 0).sum()),
        },
    )


def build_partner_loso_fixture(
    out_dir: Path | None = None,
    *,
    seed: int = 42,
) -> dict[str, Any]:
    """Split ST000883 into two OMNI-style CSVs for runnable LOSO diligence."""
    out_dir = Path(out_dir or PARTNER_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    m = load_mw_patient_matrix("ST000883", feature_map="malaria_lit")
    pos = list(m.labels[m.labels == 1].index.astype(str))
    neg = list(m.labels[m.labels == 0].index.astype(str))
    rng = np.random.default_rng(seed)
    rng.shuffle(pos)
    rng.shuffle(neg)
    a_ids = pos[: len(pos) // 2] + neg[: len(neg) // 2]
    b_ids = pos[len(pos) // 2 :] + neg[len(neg) // 2 :]

    def _write(study_id: str, ids: list[str]) -> Path:
        block = m.matrix.loc[ids].copy()
        block.insert(0, "sample_id", ids)
        block["label"] = m.labels.loc[ids].astype(int).to_numpy()
        # OMNI-ish: put label last
        cols = ["sample_id"] + [c for c in block.columns if c not in {"sample_id", "label"}] + ["label"]
        path = out_dir / f"{study_id}_omni_feature_table.csv"
        block[cols].to_csv(path, index=False)
        return path

    path_a = _write("PARTNER_SITE_A", a_ids)
    path_b = _write("PARTNER_SITE_B", b_ids)
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Diligence fixture only: stratified split of public ST000883 into two "
            "OMNI-style CSVs to demonstrate partner ingest + LOSO. NOT independent "
            "external clinical validation. Replace with real partner intensity tables "
            "for buyer diligence."
        ),
        "source_study": "ST000883",
        "feature_map": "malaria_lit",
        "seed": seed,
        "site_a": {"path": str(path_a), "n": len(a_ids), "n_pos": int(m.labels.loc[a_ids].sum())},
        "site_b": {"path": str(path_b), "n": len(b_ids), "n_pos": int(m.labels.loc[b_ids].sum())},
    }
    (out_dir / "PARTNER_LOSO_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def ensure_partner_fixture() -> dict[str, Any]:
    if FIXTURE_MANIFEST.exists():
        return json.loads(FIXTURE_MANIFEST.read_text())
    return build_partner_loso_fixture()


def run_partner_loso_diligence(
    *,
    site_a: Path | None = None,
    site_b: Path | None = None,
    disease_id: str = "malaria",
    seed: int = 42,
    max_features: int = 8,
) -> dict[str, Any]:
    """Load two OMNI-style tables (or bundled fixture) and run LOSO."""
    if site_a is None or site_b is None:
        man = ensure_partner_fixture()
        site_a = Path(man["site_a"]["path"])
        site_b = Path(man["site_b"]["path"])
        fixture_note = man.get("honesty")
    else:
        fixture_note = "User-supplied partner OMNI-style tables"
        man = None

    a = load_omni_style_csv(Path(site_a), study_id="PARTNER_SITE_A", disease_id=disease_id)
    b = load_omni_style_csv(Path(site_b), study_id="PARTNER_SITE_B", disease_id=disease_id)
    from .score import disease_signature

    try:
        sig = disease_signature(
            disease_id,
            source="hybrid",
            top_n=40,
            voc_ids=sorted(set(a.matrix.columns) | set(b.matrix.columns)),
        )
    except Exception:  # noqa: BLE001
        sig = {c: 1.0 for c in list(a.matrix.columns)[:8]}

    loso = leave_one_study_out(
        [a.log1p(), b.log1p()], signature=sig, max_features=max_features, seed=seed
    )
    return {
        "title": "Partner LOSO diligence slide",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": fixture_note,
        "fixture_manifest": man,
        "site_a": {"path": str(site_a), **a.metadata},
        "site_b": {"path": str(site_b), **b.metadata},
        "loso": loso,
        "buyer_slide": {
            "headline": "Cut the cost of wrong VOC panels — multi-site LOSO before assay spend",
            "mean_auroc_signature": loso.get("mean_auroc_signature"),
            "mean_auroc_sparse": loso.get("mean_auroc_sparse"),
            "n_shared_vocs": len(loso.get("shared_vocs") or []),
            "do_not_claim": "Not clinical diagnostic SOTA; diligence on partner-format intensity tables",
        },
    }


__all__ = [
    "PARTNER_DIR",
    "build_partner_loso_fixture",
    "ensure_partner_fixture",
    "load_omni_style_csv",
    "run_partner_loso_diligence",
]
