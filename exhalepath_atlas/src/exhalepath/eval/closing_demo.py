"""One demo that closes: paste note → top VOCs + optimism gap + paper zip.

Lead message: cut the cost of wrong VOC panels — never clinical AUROC theater.
Fixed diseases: malaria (MW ST000883) or asthma (Sci Data OVR + literature panel).
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal, Optional

from ..gcms.omni_partner import run_partner_loso_diligence
from ..gcms.paper_pack import export_paper_pack
from ..nl.patient_template import parse_patient_template
from ..oem.kit import cite_vocs_for_disease, score_sample


DEFAULT_NOTES = {
    "malaria": (
        "8yo child from endemic region with fever, chills, and headache; "
        "suspected Plasmodium falciparum malaria; non-smoker household"
    ),
    "asthma": (
        "34F with poorly controlled asthma, nocturnal wheeze, on inhaled steroids; "
        "never smoker; BMI 26"
    ),
}


def _run_optimism_gap(disease: str, out_dir: Path) -> dict[str, Any]:
    if disease == "malaria":
        from .malaria_diagnostic import evaluate_malaria_diagnostic

        report = evaluate_malaria_diagnostic(
            out_dir=out_dir / "diagnostic",
            feature_map="malaria_lit",
            max_features=6,
            seed=42,
            loso=False,
        )
        t = report.get("transferable_signature") or {}
        sk = report.get("fit_on_cohort_sparse_kbest") or {}
        return {
            "study_id": report.get("study_id"),
            "nested_auroc": t.get("nested_auroc"),
            "non_nested_auroc": t.get("non_nested_auroc"),
            "optimism_gap": t.get("optimism_gap"),
            "sparse_ceiling": sk.get("mean_test_auroc"),
            "path": str(out_dir / "diagnostic" / "MALARIA_DIAGNOSTIC_UPGRADE.md"),
            "split_sha256": (report.get("nested") or {}).get("content_sha256"),
        }
    # asthma — Sci Data OVR nested-ish via patient_diagnostic style on OVR matrix
    from ..gcms.locked_split import lock_split
    from ..gcms.scidata_samples import load_scidata_ovr_matrix
    from ..gcms.score import disease_signature, score_patients
    from sklearn.metrics import roc_auc_score
    import numpy as np

    matrix = load_scidata_ovr_matrix("asthma").log1p()
    sig = disease_signature(
        "asthma", source="hybrid", top_n=40, voc_ids=list(matrix.matrix.columns.astype(str))
    )
    manifest = lock_split(
        matrix,
        strategy="stratified_kfold",
        n_splits=5,
        seed=42,
        out_path=out_dir / "diagnostic" / "asthma_split_manifest.json",
    )
    (out_dir / "diagnostic").mkdir(parents=True, exist_ok=True)
    from .patient_diagnostic import _fold_control_relative_scores

    fold_aucs = []
    for fold in manifest.folds:
        sc = _fold_control_relative_scores(
            matrix, sig, fold.train_subject_ids, fold.test_subject_ids, method="cosine"
        )
        y_te = matrix.labels.loc[fold.test_subject_ids].to_numpy(dtype=int)
        if len(set(y_te.tolist())) >= 2:
            fold_aucs.append(float(roc_auc_score(y_te, sc.to_numpy())))
    full = score_patients(matrix, sig, method="cosine", reference="control_mean")
    y = matrix.labels.loc[full.index].astype(int)
    non_nested = float(roc_auc_score(y, full)) if y.nunique() >= 2 else None
    nested = float(np.mean(fold_aucs)) if fold_aucs else None
    # write stub report for paper pack
    stub = {
        "nested_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": (non_nested - nested)
        if (non_nested is not None and nested is not None)
        else None,
        "caveat": "Sci Data asthma OVR — negatives are other pulmonary cohorts, not healthy",
    }
    (out_dir / "diagnostic" / "ASTHMA_OVR_GAP.json").write_text(json.dumps(stub, indent=2))
    (out_dir / "diagnostic" / "METHODS.md").write_text(
        "# Methods stub — asthma OVR closing demo\n\n"
        "Nested signature AUROC on Sci Data asthma one-vs-rest. "
        "Not disease-vs-healthy clinical accuracy.\n"
    )
    return {
        "study_id": "scidata:asthma",
        "nested_auroc": nested,
        "non_nested_auroc": non_nested,
        "optimism_gap": stub["optimism_gap"],
        "sparse_ceiling": None,
        "path": str(out_dir / "diagnostic" / "ASTHMA_OVR_GAP.json"),
        "split_sha256": manifest.content_sha256,
        "caveat": stub["caveat"],
    }


def run_closing_demo(
    *,
    disease: Literal["malaria", "asthma"] = "malaria",
    note: str | None = None,
    out_dir: Path | None = None,
    include_partner_loso: bool = True,
    llm: str = "rules",
) -> dict[str, Any]:
    """End-to-end buyer demo ≤5 minutes of wall clock on bundled data."""
    out_dir = Path(out_dir or f"runs/closing_demo_{disease}")
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    text = (note or DEFAULT_NOTES[disease]).strip()
    tpl = parse_patient_template(text, llm=llm)
    if not tpl.primary_disease:
        tpl.primary_disease = disease
        tpl.warnings.append("primary_disease_defaulted_to_demo_disease")

    # Hypothesis VOCs via OEM score_sample (ledger-cited)
    # Use a small observed vector from signature directions as stand-in "panel read"
    from ..gcms.score import disease_signature

    sig = disease_signature(disease, source="hybrid", top_n=15)
    # simulate a weakly aligned observed vector for demo scoring
    observed = {k: float(v) * 0.5 for k, v in list(sig.items())[:12]}
    scored = score_sample(disease, observed, signature_source="hybrid")
    cites = cite_vocs_for_disease(
        disease, [r["voc_id"] for r in scored.get("top_vocs") or []]
    )

    gap = _run_optimism_gap(disease, out_dir)

    partner = None
    if include_partner_loso and disease == "malaria":
        try:
            partner = run_partner_loso_diligence(disease_id="malaria", seed=42)
            (out_dir / "PARTNER_LOSO.json").write_text(
                json.dumps(partner, indent=2, default=str)
            )
        except Exception as exc:  # noqa: BLE001
            partner = {"status": "error", "error": str(exc)}

    pack = export_paper_pack(
        out_dir / "diagnostic",
        out_zip=out_dir / f"{disease}_closing_paper_pack.zip",
        study_id=gap.get("study_id"),
        disease_id=disease,
        split_sha256=gap.get("split_sha256"),
        signature_source="hybrid",
        methods_notes=[
            "Closing demo pack — research enablement, not a diagnostic device claim.",
            "Lead value: cut the cost of wrong VOC panels via optimism-gap + paper zip.",
        ],
    )

    report = {
        "title": "Closing demo — paste note → VOCs + optimism gap + paper zip",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "headline": "Cut the cost of wrong VOC panels",
        "not_selling": [
            "Not clinical diagnostic SOTA",
            "Not FDA/CE cleared",
            "Do not lead with Sci Data 100% OVR or n≈35 AUROC as product proof",
        ],
        "disease": disease,
        "note": text,
        "patient_template": tpl.model_dump(),
        "top_vocs_cited": scored.get("top_vocs"),
        "sample_score": {
            "score": scored.get("score"),
            "n_overlap": scored.get("n_overlap"),
        },
        "citation_summary": {
            "n_vocs": len(cites),
            "by_grade": _count_grades(cites),
        },
        "optimism_gap": gap,
        "partner_loso": {
            "buyer_slide": (partner or {}).get("buyer_slide"),
            "honesty": (partner or {}).get("honesty"),
            "mean_auroc_signature": ((partner or {}).get("loso") or {}).get(
                "mean_auroc_signature"
            ),
            "mean_auroc_sparse": ((partner or {}).get("loso") or {}).get(
                "mean_auroc_sparse"
            ),
        }
        if partner
        else None,
        "paper_pack": pack,
        "artifacts": {
            "md": str(out_dir / "CLOSING_DEMO.md"),
            "json": str(out_dir / "CLOSING_DEMO.json"),
            "zip": pack.get("zip_path"),
        },
    }
    (out_dir / "CLOSING_DEMO.json").write_text(json.dumps(report, indent=2, default=str))
    (out_dir / "CLOSING_DEMO.md").write_text(_md(report))
    (out_dir / "PATIENT_TEMPLATE.json").write_text(tpl.model_dump_json(indent=2))
    return report


def _count_grades(cites: dict[str, Any]) -> dict[str, int]:
    out: dict[str, int] = {}
    for c in cites.values():
        g = str(c.get("evidence_grade") or "unknown")
        out[g] = out.get(g, 0) + 1
    return out


def _pct(x: Any) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.1f}%"
    except (TypeError, ValueError):
        return str(x)


def _md(report: dict[str, Any]) -> str:
    gap = report.get("optimism_gap") or {}
    partner = report.get("partner_loso") or {}
    cites = report.get("citation_summary") or {}
    lines = [
        "# Closing demo",
        "",
        f"**{report.get('headline')}**",
        "",
        f"Generated: {report.get('generated_utc')}",
        f"Disease focus: `{report.get('disease')}`",
        "",
        "> Research enablement — not a medical device. Do not sell Sci Data OVR "
        "or n≈35 AUROC as clinical proof.",
        "",
        "## 1) Pasted note → patient template",
        "",
        f"```\n{report.get('note')}\n```",
        "",
        f"Parsed primary_disease: `{((report.get('patient_template') or {}).get('primary_disease'))}`",
        "",
        "## 2) Top VOCs (ledger-cited)",
        "",
        f"Citation mix: {cites.get('by_grade')}",
        "",
        "| VOC | signature log2fc | evidence_grade | DOI |",
        "|---|---:|---|---|",
    ]
    for r in (report.get("top_vocs_cited") or [])[:12]:
        lines.append(
            f"| `{r.get('voc_id')}` | {r.get('signature_log2fc')} | "
            f"{r.get('evidence_grade')} | {r.get('doi') or '—'} |"
        )
    lines += [
        "",
        "## 3) Optimism gap (why buyers pay)",
        "",
        f"- Study: `{gap.get('study_id')}`",
        f"- Transferable nested AUROC: **{_pct(gap.get('nested_auroc'))}**",
        f"- Non-nested AUROC: **{_pct(gap.get('non_nested_auroc'))}**",
        f"- Optimism gap: **{_pct(gap.get('optimism_gap'))}**",
        f"- Fit-on-cohort sparse ceiling: **{_pct(gap.get('sparse_ceiling'))}** "
        "(ceiling, not the product claim)",
        f"- Locked split sha256: `{gap.get('split_sha256')}`",
        "",
        "## 4) Partner LOSO diligence",
        "",
    ]
    if partner:
        lines += [
            f"- Signature LOSO AUROC: **{_pct(partner.get('mean_auroc_signature'))}**",
            f"- Sparse LOSO AUROC: **{_pct(partner.get('mean_auroc_sparse'))}**",
            f"- Honesty: {partner.get('honesty')}",
            "",
        ]
    else:
        lines.append("- (skipped for this disease)\n")
    pack = report.get("paper_pack") or {}
    lines += [
        "## 5) Paper zip",
        "",
        f"- Zip: `{pack.get('zip_path')}`",
        f"- Zip sha256: `{pack.get('zip_sha256')}`",
        f"- Files: {pack.get('n_files')}",
        "",
        "## Pitch line",
        "",
        "We cut the cost of being wrong about VOCs — locked nested eval, "
        "ledger-cited hypotheses, and a deposit-ready zip — before you spend "
        "the next assay cohort.",
        "",
    ]
    return "\n".join(lines)


__all__ = ["DEFAULT_NOTES", "run_closing_demo"]
