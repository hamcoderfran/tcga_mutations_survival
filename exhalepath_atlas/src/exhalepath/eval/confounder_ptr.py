"""Non-malaria metadata-rich confounder pack: ST003200 PTR + Sci Data age/sex.

ST003200 (n=504 healthy PTR-TOF-MS) has Sex + Smoking Status + Age — the open
breath cohort where confounder effect sizes can be measured cleanly.
Sci Data 2024 pulmonary GC-MS has age/sex (no smoking) via one-vs-rest strata.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from ..config import DATA_DIR
from ..gcms.metrics import confounder_proxy_auroc, stratified_auroc
from ..gcms.patient_matrix import PatientVOCMatrix
from ..ingest.real_breath_corpus import _map_voc, _ms_rows
from ..gcms.schema import SampleRecord

MW_DIR = DATA_DIR / "datasources" / "metabolomics"


def load_st003200_ptr_matrix() -> tuple[PatientVOCMatrix, pd.DataFrame]:
    """Patient×VOC matrix + demographics table for ST003200."""
    path = MW_DIR / "ST003200" / "mwtab.json"
    if not path.exists():
        raise FileNotFoundError(path)
    data = json.loads(path.read_text())
    rows = _ms_rows(data)
    df = pd.DataFrame(rows)
    df["voc_id"] = df["metabolite"].map(_map_voc)
    df = df.dropna(subset=["voc_id"])
    df["sample"] = df["sample"].astype(str)
    piv = df.pivot_table(
        index="sample", columns="voc_id", values="intensity", aggfunc="median"
    )
    piv = piv.dropna(axis=1, how="all")
    for col in piv.columns:
        med = float(piv[col].median(skipna=True))
        if np.isnan(med):
            med = 0.0
        piv[col] = piv[col].fillna(med)

    demo_rows = []
    for s in data.get("SUBJECT_SAMPLE_FACTORS") or []:
        sid = str(s.get("Sample ID") or "")
        fac = s.get("Factors") or {}
        add = s.get("Additional sample data") or {}
        sex = str(fac.get("Sex") or "").strip()
        if sex.lower() == "female":
            sex = "Female"
        elif sex.lower() == "male":
            sex = "Male"
        smoke = str(fac.get("Smoking Status") or "").strip()
        age_raw = add.get("Age (years)") or add.get("Age")
        try:
            age = float(age_raw) if age_raw not in (None, "", "-") else None
        except (TypeError, ValueError):
            age = None
        demo_rows.append(
            {
                "sample_id": sid,
                "sex": sex if sex in ("Male", "Female") else None,
                "smoking_status": smoke if smoke and smoke != "-" else None,
                "age": age,
            }
        )
    demo = pd.DataFrame(demo_rows).drop_duplicates("sample_id").set_index("sample_id")
    common = [i for i in piv.index if i in demo.index]
    piv = piv.loc[common]
    demo = demo.loc[common]
    # dummy labels (healthy cohort) — required by PatientVOCMatrix
    labels = pd.Series(0, index=piv.index, dtype=int, name="label")
    samples = [
        SampleRecord(
            sample_id=sid,
            subject_id=sid,
            study_id="ST003200",
            label="control",
            disease_id="healthy",
            modality="ptr_ms",
            factors_raw=json.dumps(demo.loc[sid].to_dict()),
        )
        for sid in common
    ]
    matrix = PatientVOCMatrix(
        study_id="ST003200",
        disease_id="healthy",
        disease_name="Healthy breath PTR",
        modality="ptr_ms",
        unit="ppbv",
        matrix=piv.astype(float),
        labels=labels,
        samples=samples,
        metadata={
            "n_subjects": len(common),
            "n_voc_features": int(piv.shape[1]),
            "doi": "10.1007/s11306-024-02139-6",
            "role": "confounder_effect_size_reference",
        },
    )
    return matrix, demo


def evaluate_confounder_ptr(
    *,
    out_dir: Path | None = None,
    include_scidata: bool = True,
) -> dict[str, Any]:
    """VOC→confounder probes (smoking/sex) + age tertile summaries; Sci Data strata."""
    matrix, demo = load_st003200_ptr_matrix()
    X = np.log1p(matrix.matrix.clip(lower=0).to_numpy(dtype=float))

    # Smoking: current vs never (drop ex / missing)
    smoke = demo["smoking_status"].astype(str)
    smoke_mask = smoke.isin(["Current smoker", "Non smoker"])
    y_smoke = (smoke.loc[smoke_mask] == "Current smoker").astype(int).to_numpy()
    X_smoke = X[smoke_mask.to_numpy()]
    smoke_auc = confounder_proxy_auroc(X_smoke, y_smoke)

    # Sex
    sex = demo["sex"]
    sex_mask = sex.isin(["Male", "Female"])
    y_sex = (sex.loc[sex_mask] == "Male").astype(int).to_numpy()
    X_sex = X[sex_mask.to_numpy()]
    sex_auc = confounder_proxy_auroc(X_sex, y_sex)

    # Age tertiles — report mean VOC z by tertile (effect size, not disease AUROC)
    age = demo["age"]
    age_ok = age.dropna()
    age_tertile = pd.qcut(age_ok, 3, labels=["young", "mid", "older"])
    voc_by_age: dict[str, Any] = {}
    X_age = matrix.matrix.loc[age_ok.index]
    for level in ["young", "mid", "older"]:
        idx = age_tertile[age_tertile == level].index
        if len(idx) < 5:
            continue
        means = X_age.loc[idx].mean(axis=0)
        voc_by_age[level] = {
            "n": int(len(idx)),
            "top_mean_vocs": means.sort_values(ascending=False).head(5).round(3).to_dict(),
        }

    # Stratified "dummy" — show metadata coverage for smoking/sex/age
    coverage = {
        "n": int(len(demo)),
        "sex": demo["sex"].value_counts(dropna=False).to_dict(),
        "smoking_status": demo["smoking_status"].value_counts(dropna=False).to_dict(),
        "age_non_null": int(age.notna().sum()),
        "age_median": float(age.median()) if age.notna().any() else None,
    }

    scidata_block: Optional[dict[str, Any]] = None
    if include_scidata:
        try:
            from ..gcms.scidata_samples import list_scidata_cohorts, load_scidata_ovr_matrix
            from ..gcms.score import disease_signature, score_patients

            cohorts = [c for c in list_scidata_cohorts() if c.get("cohort")]
            strata_rows = []
            for c in cohorts[:3]:
                if not c.get("available"):
                    continue
                cid = str(c["cohort"])
                try:
                    m = load_scidata_ovr_matrix(cid)
                except Exception:  # noqa: BLE001
                    continue
                age_map = {s.subject_id: s.age for s in m.samples}
                sex_map = {s.subject_id: s.sex for s in m.samples}
                try:
                    sig = disease_signature(
                        m.disease_id or "bronchiectasis",
                        source="hybrid",
                        top_n=30,
                        voc_ids=list(m.matrix.columns.astype(str)),
                    )
                except Exception:  # noqa: BLE001
                    sig = {col: 1.0 for col in list(m.matrix.columns)[:5]}
                scores = score_patients(
                    m.log1p(), sig, method="cosine", reference="control_mean"
                )
                y = m.labels.loc[scores.index].astype(int)
                strata = {
                    "age": [age_map.get(sid) for sid in scores.index],
                    "sex": [sex_map.get(sid) for sid in scores.index],
                }
                st = stratified_auroc(y.tolist(), scores.tolist(), strata, min_n=8)
                strata_rows.append(
                    {
                        "cohort": m.study_id,
                        "n": int(m.metadata.get("n_subjects") or len(y)),
                        "overall_auroc": st.get("overall"),
                        "strata": st.get("strata"),
                        "note": (
                            "Sci Data OVR pulmonary — age/sex when present; "
                            "no smoking; no healthy arm"
                        ),
                    }
                )
            scidata_block = {
                "n_cohorts_reported": len(strata_rows),
                "cohorts": strata_rows,
                "smoking_available": False,
            }
        except Exception as exc:  # noqa: BLE001
            scidata_block = {"error": str(exc)}

    report: dict[str, Any] = {
        "title": "Confounder / demographics pack (ST003200 PTR + Sci Data GC-MS)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "ST003200 is healthy-only — VOC→smoking/sex AUROC measures confounder "
            "signal in the feature space, not disease discrimination. Sci Data OVR "
            "strata are pulmonary case-mix, not case–control vs healthy."
        ),
        "st003200": {
            "n_subjects": matrix.metadata.get("n_subjects"),
            "n_voc_features": matrix.metadata.get("n_voc_features"),
            "voc_ids": list(matrix.matrix.columns.astype(str)),
            "coverage": coverage,
            "confounder_proxy_auroc": {
                "smoking_current_vs_never": smoke_auc,
                "sex_male_vs_female": sex_auc,
                "n_smoking": int(smoke_mask.sum()),
                "n_sex": int(sex_mask.sum()),
            },
            "voc_means_by_age_tertile": voc_by_age,
        },
        "scidata_gcms_age_sex": scidata_block,
        "caveats": [
            "ST003200: no disease labels — use as confounder effect-size reference only.",
            "Sci Data: one-vs-rest pulmonary cohorts; smoking metadata absent.",
            "Do not treat confounder-proxy AUROC as diagnostic performance.",
        ],
    }

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "CONFOUNDER_PTR.json").write_text(
            json.dumps(report, indent=2, default=str)
        )
        (out_dir / "CONFOUNDER_PTR.md").write_text(_md(report))
        demo.to_csv(out_dir / "ST003200_demographics.csv")
        matrix.matrix.to_csv(out_dir / "ST003200_voc_matrix.csv")
        report["artifacts"] = {
            "json": str(out_dir / "CONFOUNDER_PTR.json"),
            "md": str(out_dir / "CONFOUNDER_PTR.md"),
        }
    return report


def _md(report: dict[str, Any]) -> str:
    s = report.get("st003200") or {}
    c = (s.get("confounder_proxy_auroc") or {})
    lines = [
        "# Confounder / demographics pack",
        "",
        f"Generated: {report.get('generated_utc')}",
        "",
        f"> {report.get('honesty')}",
        "",
        "## ST003200 healthy PTR (n=504)",
        "",
        f"- Features: {s.get('n_voc_features')} VOCs — {', '.join(s.get('voc_ids') or [])}",
        f"- Smoking current vs never — confounder-proxy AUROC: **{_pct(c.get('smoking_current_vs_never'))}** "
        f"(n={c.get('n_smoking')})",
        f"- Sex male vs female — confounder-proxy AUROC: **{_pct(c.get('sex_male_vs_female'))}** "
        f"(n={c.get('n_sex')})",
        "",
        "## Sci Data GC-MS age/sex strata",
        "",
    ]
    sc = report.get("scidata_gcms_age_sex") or {}
    if sc.get("error"):
        lines.append(f"- Unavailable: {sc['error']}")
    else:
        lines.append(f"- Cohorts reported: {sc.get('n_cohorts_reported')}")
        lines.append(f"- Smoking available: {sc.get('smoking_available')}")
        for row in sc.get("cohorts") or []:
            lines.append(
                f"- `{row.get('cohort')}` n={row.get('n')} overall AUROC={_pct(row.get('overall_auroc'))}"
            )
    lines += ["", "## Caveats", ""]
    for x in report.get("caveats") or []:
        lines.append(f"- {x}")
    lines.append("")
    return "\n".join(lines)


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


__all__ = ["evaluate_confounder_ptr", "load_st003200_ptr_matrix"]
