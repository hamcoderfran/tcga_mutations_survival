"""Industry-standard smoke + benchmark pack with honest commercial framing.

Metrics follow breathomics / TRIPOD+AI / diagnostic-ML norms:
- nested vs non-nested AUROC + optimism gap (JBR 2024)
- AUPRC, Youden sens/spec, bootstrap CI95
- stratified AUROC when confounder metadata exist
- directional holdout (literature / public breath)
- coverage audit of open compound universe
- paper-pack / deposit scaffold completeness

Commercial framing is R&D enablement / de-risking for VOC programs —
**not** a clinical diagnostic valuation claim.
"""

from __future__ import annotations

import json
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


def _safe(name: str, fn: Callable[[], Any]) -> dict[str, Any]:
    t0 = time.perf_counter()
    try:
        payload = fn()
        return {
            "name": name,
            "ok": True,
            "seconds": round(time.perf_counter() - t0, 3),
            "payload": payload,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "name": name,
            "ok": False,
            "seconds": round(time.perf_counter() - t0, 3),
            "error": str(exc),
            "traceback": traceback.format_exc(limit=4),
        }


def _pct(x: Any, digits: int = 1) -> str:
    if x is None:
        return "—"
    try:
        return f"{100 * float(x):.{digits}f}%"
    except (TypeError, ValueError):
        return str(x)


def run_industry_pack(
    *,
    out_dir: Path = Path("runs/industry_pack"),
    quick: bool = False,
) -> dict[str, Any]:
    """Run comprehensive smoke + industry metrics; write reports under out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    steps: list[dict[str, Any]] = []

    # 1) Patient-level MW diagnostic (industry gold for breath ML)
    def _mw():
        from .patient_diagnostic import evaluate_patient_diagnostic

        return evaluate_patient_diagnostic(
            study_id="ST000883",
            signature_source="hybrid",
            out_dir=out_dir / "ST000883",
            n_splits=5,
            seed=42,
        )

    steps.append(_safe("patient_diagnostic_ST000883", _mw))

    # 2) Signature benchmark (hybrid/stack/literature × cosine/dot)
    def _bench():
        from .patient_diagnostic import run_signature_benchmark

        return run_signature_benchmark(study_id="ST000883", out_dir=out_dir / "benchmark")

    steps.append(_safe("signature_benchmark_ST000883", _bench))

    # 3) Sci Data per-sample (all cohorts unless quick)
    def _scidata():
        from .scidata_samples import evaluate_all_scidata_cohorts, evaluate_scidata_samples

        if quick:
            return {
                "mode": "quick_asthma",
                "asthma": evaluate_scidata_samples(
                    positive_cohort="asthma",
                    out_dir=out_dir / "scidata" / "asthma",
                    n_splits=3,
                    make_paper_pack=True,
                ),
            }
        return evaluate_all_scidata_cohorts(out_dir=out_dir / "scidata")

    steps.append(_safe("scidata_per_sample", _scidata))

    # 4) Stack holdout
    def _stack():
        from .stack_holdout import evaluate_stack_holdout

        return evaluate_stack_holdout(out_dir=out_dir / "stack_holdout")

    steps.append(_safe("stack_holdout", _stack))

    # 5) Public breath directional
    def _pub():
        from .public_breath import evaluate_public_breath

        return evaluate_public_breath(
            top_k=15, mode="hybrid", out_dir=out_dir / "public_breath"
        )

    steps.append(_safe("public_breath", _pub))

    # 6) Coverage audit (offline)
    def _cov():
        from .coverage_audit import run_coverage_audit

        return run_coverage_audit(expand=True, offline=True)

    steps.append(_safe("coverage_audit_offline", _cov))

    # 7) Lit compare (lighter defaults)
    if not quick:

        def _lit():
            from .lit_compare import run_lit_demo_disease100

            return run_lit_demo_disease100(
                out_dir=out_dir / "lit_compare",
                demo_max_patients=120,
                n_diseases=40,
            )

        steps.append(_safe("lit_compare", _lit))

    # 8) MetaboLights + paper pack smoke
    def _deposit():
        from ..gcms import (
            export_metabolights_bundle,
            export_paper_pack,
            load_mw_patient_matrix,
            load_scidata_ovr_matrix,
        )

        mw = load_mw_patient_matrix("ST000883")
        mtbls = export_metabolights_bundle(mw, out_dir / "metabolights" / mw.study_id)
        sci = load_scidata_ovr_matrix("asthma", mapped_vocs_only=True)
        mtbls2 = export_metabolights_bundle(sci, out_dir / "metabolights" / sci.study_id)
        pack = export_paper_pack(
            out_dir / "ST000883",
            study_id="ST000883",
            disease_id="malaria",
            split_sha256=(
                ((steps[0].get("payload") or {}).get("nested") or {}).get("content_sha256")
            ),
        )
        return {
            "metabolights_mw": {k: str(v) for k, v in mtbls.items()},
            "metabolights_scidata": {k: str(v) for k, v in mtbls2.items()},
            "paper_pack": pack,
        }

    steps.append(_safe("deposit_scaffolds", _deposit))

    # Build scorecard from successful steps
    by_name = {s["name"]: s for s in steps}
    rows = _scorecard_rows(by_name)
    commercial = _commercial_brief(rows, by_name)
    prune = _prune_memo()

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "product": "voc-breath / ExhalePath Atlas",
        "claim_boundary": (
            "Cutting-edge open research-enablement stack for exhaled VOC / GC-MS R&D. "
            "Not a clinical diagnostic SOTA claim; not FDA/CE clearance."
        ),
        "steps": [
            {k: v for k, v in s.items() if k != "payload"}
            | (
                {"summary": _summarize_payload(s["name"], s.get("payload"))}
                if s.get("ok")
                else {}
            )
            for s in steps
        ],
        "industry_scorecard": {"rows": rows},
        "commercial": commercial,
        "draft_branch_prune": prune,
        "smoke_pass": all(s["ok"] for s in steps),
        "n_steps_ok": sum(1 for s in steps if s["ok"]),
        "n_steps": len(steps),
    }

    # Attach full payloads separately (large)
    payloads = {s["name"]: s.get("payload") for s in steps if s.get("ok")}
    (out_dir / "INDUSTRY_PACK.json").write_text(json.dumps(report, indent=2, default=str))
    (out_dir / "INDUSTRY_PACK_PAYLOADS.json").write_text(
        json.dumps(payloads, indent=2, default=str)
    )
    (out_dir / "INDUSTRY_PACK.md").write_text(_md_pack(report))
    (out_dir / "BUYER_BRIEF.md").write_text(_md_buyer(report))
    (out_dir / "DRAFT_BRANCH_PRUNE.md").write_text(_md_prune(prune))
    return report


def _summarize_payload(name: str, payload: Any) -> dict[str, Any]:
    if not isinstance(payload, dict):
        return {"type": type(payload).__name__}
    if name.startswith("patient_diagnostic"):
        m = payload.get("metrics") or {}
        n = payload.get("nested") or {}
        return {
            "auroc": m.get("auroc"),
            "auprc": m.get("auprc"),
            "nested_auroc": n.get("mean_test_auroc"),
            "optimism_gap": n.get("optimism_gap"),
            "logistic_nested": (n.get("logistic_baseline") or {}).get("mean_test_auroc"),
            "split_sha256": (n.get("content_sha256") or "")[:16],
        }
    if name == "signature_benchmark_ST000883":
        return {
            "rows": [
                {
                    "signature": r.get("signature"),
                    "method": r.get("method"),
                    "auroc": r.get("auroc"),
                    "nested_auroc": r.get("nested_auroc"),
                }
                for r in (payload.get("rows") or [])
                if r.get("ok")
            ]
        }
    if name == "scidata_per_sample":
        if "cohorts" in payload:
            return {"cohorts": payload.get("cohorts")}
        a = (payload.get("asthma") or {}).get("metrics") or {}
        return {"asthma_auroc": a.get("auroc"), "n": a.get("n_subjects")}
    if name == "stack_holdout":
        pb = payload.get("public_breath") or {}
        lit = payload.get("literature") or {}
        p10 = payload.get("priority10") or {}
        brk = payload.get("patient_break") or {}
        return {
            "public_breath_stack": pb.get("stack_mean_directional"),
            "public_breath_hybrid": pb.get("hybrid_mean_directional"),
            "literature_stack": lit.get("stack_mean_directional"),
            "priority10_stack": p10.get("stack_mean_directional"),
            "patient_break_hard": brk.get("hard_pass"),
            "patient_break_soft": brk.get("soft_pass"),
            "sota": (payload.get("sota_assessment") or {}).get("is_cutting_edge_systems_stack"),
        }
    if name == "public_breath":
        o = payload.get("overall") or {}
        return {
            "directional": o.get("mean_directional_accuracy"),
            "recall_at_k": o.get("mean_elevated_recall_at_k"),
            "n_cases": o.get("n_cases"),
        }
    if name == "coverage_audit_offline":
        m = payload.get("metrics") or {}
        return {
            "open_compound_coverage_pct": m.get("open_compound_coverage_pct"),
            "meets_99pct": m.get("meets_99pct_open_compound_target"),
            "secured": m.get("open_compound_secured_n"),
            "universe": m.get("open_compound_universe_n"),
        }
    if name == "lit_compare":
        return {
            "lit_concordance_pct": (payload.get("literature") or {}).get(
                "mean_concordance_pct"
            ),
            "n_diseases": (payload.get("disease100") or {}).get("n_diseases"),
        }
    if name == "deposit_scaffolds":
        return {
            "paper_pack_files": (payload.get("paper_pack") or {}).get("n_files"),
            "zip": (payload.get("paper_pack") or {}).get("zip_path"),
        }
    return {"keys": list(payload.keys())[:12]}


def _scorecard_rows(by_name: dict[str, dict]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(metric: str, value: Any, standard: str, why: str, status: str = "measured"):
        rows.append(
            {
                "metric": metric,
                "value": value if isinstance(value, str) else (
                    value if value is not None else "FAIL/unavailable"
                ),
                "standard": standard,
                "why_it_matters": why,
                "status": status,
            }
        )

    mw = by_name.get("patient_diagnostic_ST000883") or {}
    if mw.get("ok"):
        p = mw["payload"]
        m = p.get("metrics") or {}
        n = p.get("nested") or {}
        add(
            "Nested patient AUROC (ST000883 hybrid)",
            _pct(n.get("mean_test_auroc")),
            "TRIPOD+AI / JBR nested CV",
            "Patient-level holdout — unit of clinical inference; beats peak-level CV optimism",
        )
        add(
            "Optimism gap (non-nested − nested)",
            _pct(n.get("optimism_gap")),
            "JBR 2024 breath-ML leakage audit",
            "Quantifies overfit risk buyers must price into VOC programs",
        )
        add(
            "Same-feature nested logistic ceiling",
            _pct((n.get("logistic_baseline") or {}).get("mean_test_auroc")),
            "Internal ML baseline",
            "Upper bound when labels are fit on this tiny cohort — mechanism signatures are transferable, not maxed here",
        )
        add(
            "AUPRC / Youden sens-spec",
            f"AUPRC={_pct(m.get('auprc'))}; sens={_pct(m.get('sensitivity'))}; spec={_pct(m.get('specificity'))}",
            "STARD / diagnostic operating point",
            "Discrimination + thresholded clinical operating characteristics",
        )
        add(
            "Locked split SHA256",
            (n.get("content_sha256") or "")[:24] + "…",
            "Preregistration / audit trail",
            "Citeable split integrity for partner diligence",
        )
    else:
        add("Patient diagnostic ST000883", "FAILED", "required", mw.get("error") or "", "fail")

    sci = by_name.get("scidata_per_sample") or {}
    if sci.get("ok"):
        p = sci["payload"]
        if "cohorts" in p:
            vals = [
                f"{c['cohort']}={_pct(c.get('auroc'))}"
                for c in p["cohorts"]
                if c.get("ok")
            ]
            add(
                "Sci Data per-sample OVR AUROC",
                "; ".join(vals) or "—",
                "Multi-cohort peak-table ML",
                "Per-sample (not cohort-mean) discrimination across pulmonary cohorts; NO healthy arm",
            )
        else:
            a = (p.get("asthma") or {}).get("metrics") or {}
            add(
                "Sci Data asthma OVR AUROC",
                _pct(a.get("auroc")),
                "Multi-cohort peak-table ML",
                "Per-sample adapter smoke; one-vs-rest caveat applies",
            )

    st = by_name.get("stack_holdout") or {}
    if st.get("ok"):
        p = st["payload"]
        pb = p.get("public_breath") or {}
        lit = p.get("literature") or {}
        p10 = p.get("priority10") or {}
        brk = p.get("patient_break") or {}
        add(
            "Stack vs hybrid directional (public breath)",
            f"stack={_pct(pb.get('stack_mean_directional'))} hybrid={_pct(pb.get('hybrid_mean_directional'))}",
            "Systems holdout",
            "Multi-head fusion vs hybrid baseline on public panels",
        )
        add(
            "Stack directional (literature / priority-10)",
            f"lit={_pct(lit.get('stack_mean_directional'))}; p10={_pct(p10.get('stack_mean_directional'))}",
            "Systems holdout",
            "Directional panels — treat as systems evidence; partial prior overlap possible",
        )
        add(
            "PatientTemplate adversarial break",
            f"hard={brk.get('hard_pass')}/{brk.get('hard_total')} soft={brk.get('soft_pass')}/{brk.get('soft_total')}",
            "Adversarial / phenotype stress",
            "Stack must not collapse under comorbidity/smoking template attacks",
        )

    pub = by_name.get("public_breath") or {}
    if pub.get("ok"):
        o = (pub["payload"] or {}).get("overall") or {}
        add(
            "Public breath directional accuracy",
            _pct(o.get("mean_directional_accuracy")),
            "Literature/public panel concordance",
            "Hypothesis-generation accuracy on curated elevated/suppressed panels",
        )

    cov = by_name.get("coverage_audit_offline") or {}
    if cov.get("ok"):
        m = (cov["payload"] or {}).get("metrics") or {}
        add(
            "Open-compound coverage vs VOLATILOME universe",
            f"{m.get('open_compound_coverage_pct')}% "
            f"({m.get('open_compound_secured_n')}/{m.get('open_compound_universe_n')})",
            "Open data diligence",
            "Honest ≥99% target on *open* compound universe — not all gated industry atlases",
        )

    lit = by_name.get("lit_compare") or {}
    if lit.get("ok"):
        add(
            "Literature concordance (disease suite)",
            f"{(lit['payload'].get('literature') or {}).get('mean_concordance_pct')}%",
            "Panel overlay",
            "Directional agreement vs published breath panels (circularity risk flagged)",
        )

    dep = by_name.get("deposit_scaffolds") or {}
    if dep.get("ok"):
        pack = (dep["payload"] or {}).get("paper_pack") or {}
        add(
            "Paper pack completeness",
            f"{pack.get('n_files')} files · zip sha {(pack.get('zip_sha256') or '')[:12]}…",
            "Publication / partner diligence pack",
            "One zip: ROC + Methods + overlay + split hash — hours→minutes for BD packs",
        )

    add(
        "Clinical diagnostic SOTA claim",
        "NONE — research enablement only",
        "Regulatory honesty",
        "Required for credible enterprise sales; overclaim destroys diligence",
        "boundary",
    )
    return rows


def _commercial_brief(rows: list[dict], by_name: dict) -> dict[str, Any]:
    """Honest commercial narrative — R&D OS, not clinical unicorn math."""
    return {
        "positioning": (
            "ExhalePath / voc-breath is an open **VOC R&D operating system**: "
            "mechanism-aware disease→VOC hypothesis generation + patient-level GC-MS "
            "eval harness + deposit/reporting scaffolds. Buyers are breath-biopsy platforms, "
            "pharma biomarker teams, and metabolomics CROs entering VOCs — not hospital IVD."
        ),
        "what_is_actually_revolutionary": [
            "Patient-level nested AUROC + locked split SHA256 on public GC-MS (rare in open stacks)",
            "Optimism-gap reporting that prevents $1–10M wasted follow-ups on peak-level CV mirages",
            "Sci Data per-sample adapters (not just cohort means) with age/sex strata + paper zip",
            "20-head mechanism stack with anti-dilution fusion + epistemic UQ for zero-shot diseases",
            "MetaboLights mzTab-M/ISA-Tab scaffolds + TRIPOD stub — BD diligence in one command",
            "Secure open-coverage audit with host allowlist (supply-chain hygiene for enterprise IT)",
        ],
        "what_is_NOT_revolutionary": [
            "Not multi-site prospective clinical AUROC SOTA",
            "Not FDA/CE cleared; not a medical device",
            "Sci Data OVR AUROCs can look 'perfect' because cohorts are chemically distinct pulmonary diseases without healthy controls — do not sell as clinical accuracy",
            "Literature directional % can partly overlap atlas priors (circularity)",
        ],
        "buyer_pain_and_savings": [
            {
                "pain": "Peak-level CV / leakage → false greenlight of VOC panels",
                "cost_range_usd": "0.5M–8M per failed biomarker program (assay, cohort, BD time)",
                "our_lever": "Nested patient AUROC + optimism gap + locked splits",
            },
            {
                "pain": "Months assembling Methods / TRIPOD / MetaboLights deposit packs",
                "cost_range_usd": "50k–400k FTE + delayed partnership diligence",
                "our_lever": "voc export-paper-pack + export-metabolights in minutes",
            },
            {
                "pain": "Zero-shot rare disease VOC hypotheses require PhD weeks",
                "cost_range_usd": "80k–250k per indication scout",
                "our_lever": "stack / zero-shot mechanism heads + literature overlay",
            },
            {
                "pain": "Gated Owlstone-class atlases block early scouting",
                "cost_range_usd": "license + delay; opportunity cost",
                "our_lever": "open HBDB/Sci Data/MW/HMDB alt-source matrix (ds15) for pre-license triage",
            },
        ],
        "credible_deal_math": {
            "note": (
                "Billions of enterprise value require either (a) platform adoption across "
                "many pharma programs or (b) clinical IVD clearance with reimbursable use. "
                "Today's evidence supports (a)-style **R&D software / partnership** economics, "
                "not (b). Do not pitch a $B valuation on ST000883 n≈35 AUROC."
            ),
            "near_term_acv_usd": {
                "single_pharma_biomarker_seat": "150k–600k / year",
                "breath_platform_oem_embed": "0.5M–3M / year + milestone",
                "CRO white-label eval suite": "100k–500k / year",
            },
            "path_to_large_exit": [
                "Become default pre-clinical VOC hypothesis + eval layer for 3–5 breath platforms",
                "Accumulate multi-site locked-split benchmarks partners cannot ignore",
                "Optional later: IVD partnership where *their* prospective data carry clinical claims",
            ],
            "tam_context": (
                "Breath biopsy / exhaled biomarker tooling sits inside a multi-billion "
                "non-invasive diagnostics + pharma biomarker adjacent market; "
                "software attach rates are a small slice — sell acceleration and de-risking, "
                "not 'we are the diagnostic'."
            ),
        },
        "pitch_one_liner": (
            "We cut the cost of being wrong about VOCs: locked patient-level eval, "
            "mechanism hypotheses for unseen diseases, and deposit-ready paper packs — "
            "so your VOC program spends money on real signal, not optimistic CV."
        ),
    }


def _prune_memo() -> dict[str, Any]:
    return {
        "cherry_picked": [
            "secure_fetch.py + http allowlist (#14)",
            "coverage_audit.py + COVERAGE_AUDIT artifacts (#14)",
            "lit_compare.py + lit_compare artifacts (#13) — NOT priors/predict wholesale",
            "post-blend smoking/age exo fix in predict.py (surgical from #13)",
            "ds15_alt_breath_sources + catalogs (#17)",
            "HBDB 60-disease JSON + extract_hbdb_zenodo.py (#16)",
            "HMDB Wishart mirror constant + breath metabolite extracts (#18)",
        ],
        "skipped_as_regressive": [
            "Wholesale cli.py / predict.py from draft tips (would strip zero-shot hardening)",
            "disease_voc_priors.json replacements (#13/#15/#16)",
            "Recalibrated voc_calibrator.joblib from draft tips",
            "model-improve-100disease full merge (#15)",
            "Large optional Tedlar xlsx / hbdb_eval zip raw dumps (catalogs kept)",
        ],
        "action": "Close or leave draft PRs #13–#19; do not merge wholesale into main",
    }


def _md_pack(report: dict[str, Any]) -> str:
    lines = [
        "# Industry benchmark pack",
        "",
        f"Generated: {report['generated_utc']}",
        "",
        f"> {report['claim_boundary']}",
        "",
        f"Smoke: **{report['n_steps_ok']}/{report['n_steps']}** steps OK "
        f"({'PASS' if report['smoke_pass'] else 'PARTIAL'})",
        "",
        "## Scorecard",
        "",
        "| Metric | Value | Standard | Why it matters |",
        "|---|---|---|---|",
    ]
    for r in (report.get("industry_scorecard") or {}).get("rows") or []:
        lines.append(
            f"| {r['metric']} | {r['value']} | {r['standard']} | {r['why_it_matters']} |"
        )
    lines += ["", "## Step timings", "", "| Step | OK | Seconds |", "|---|---|---:|"]
    for s in report.get("steps") or []:
        lines.append(f"| {s['name']} | {s['ok']} | {s.get('seconds')} |")
        if not s.get("ok"):
            lines.append(f"| ↳ error | {s.get('error')} | |")
    lines += [
        "",
        "## Draft branch prune",
        "",
        "See `DRAFT_BRANCH_PRUNE.md` and `BUYER_BRIEF.md`.",
        "",
    ]
    return "\n".join(lines)


def _md_buyer(report: dict[str, Any]) -> str:
    c = report.get("commercial") or {}
    lines = [
        "# Buyer brief — selling VOC R&D enablement (honest)",
        "",
        f"**One-liner:** {c.get('pitch_one_liner')}",
        "",
        f"**Positioning:** {c.get('positioning')}",
        "",
        "## What is actually revolutionary (evidence-backed)",
        "",
    ]
    for x in c.get("what_is_actually_revolutionary") or []:
        lines.append(f"- {x}")
    lines += ["", "## What is NOT revolutionary (do not overclaim)", ""]
    for x in c.get("what_is_NOT_revolutionary") or []:
        lines.append(f"- {x}")
    lines += ["", "## Pain → savings levers", "", "| Pain | Cost of status quo | Our lever |", "|---|---|---|"]
    for row in c.get("buyer_pain_and_savings") or []:
        lines.append(f"| {row['pain']} | {row['cost_range_usd']} | {row['our_lever']} |")
    deal = c.get("credible_deal_math") or {}
    lines += [
        "",
        "## Credible deal math (not fake $B from n=35 AUROC)",
        "",
        deal.get("note") or "",
        "",
        "### Near-term ACV",
        "",
    ]
    for k, v in (deal.get("near_term_acv_usd") or {}).items():
        lines.append(f"- **{k}**: {v}")
    lines += ["", "### Path to larger outcomes", ""]
    for x in deal.get("path_to_large_exit") or []:
        lines.append(f"- {x}")
    lines += [
        "",
        f"TAM context: {deal.get('tam_context')}",
        "",
        "## Live scorecard snapshot",
        "",
        "| Metric | Value |",
        "|---|---|",
    ]
    for r in (report.get("industry_scorecard") or {}).get("rows") or []:
        lines.append(f"| {r['metric']} | {r['value']} |")
    lines += [
        "",
        "## Diligence commands",
        "",
        "```bash",
        "voc eval-industry-pack",
        "voc eval-patient-diagnostic --study ST000883 --signature hybrid",
        "voc eval-scidata-samples --all",
        "voc eval-coverage --offline",
        "voc export-metabolights --study ST000883",
        "```",
        "",
    ]
    return "\n".join(lines)


def _md_prune(prune: dict[str, Any]) -> str:
    lines = [
        "# Draft branch prune memo",
        "",
        "## Cherry-picked (additive)",
        "",
    ]
    for x in prune.get("cherry_picked") or []:
        lines.append(f"- {x}")
    lines += ["", "## Skipped (regressive / conflicting)", ""]
    for x in prune.get("skipped_as_regressive") or []:
        lines.append(f"- {x}")
    lines += ["", f"**Action:** {prune.get('action')}", ""]
    return "\n".join(lines)


__all__ = ["run_industry_pack"]
