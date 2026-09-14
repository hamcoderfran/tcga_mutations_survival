"""Locked patient-level truth leaderboard (Pillar 1).

Aggregates existing GC-MS diagnostic artifacts into one reproducible
leaderboard keyed by study × signature × nested AUROC × split SHA256.

Honesty: research enablement only — not clinical validation or FDA claims.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import DATA_DIR, PACKAGE_ROOT


def _knowledge_roots() -> list[Path]:
    roots = []
    for base in (DATA_DIR, PACKAGE_ROOT / "data"):
        p = base / "knowledge" / "gcms_diagnostic"
        if p.exists():
            roots.append(p)
    # de-dupe while preserving order
    seen: set[Path] = set()
    out: list[Path] = []
    for r in roots:
        rp = r.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(r)
    return out or [DATA_DIR / "knowledge" / "gcms_diagnostic"]


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def _split_sha_from_report(report: dict[str, Any]) -> str | None:
    nested = report.get("nested") or {}
    if nested.get("content_sha256"):
        return str(nested["content_sha256"])
    for key in ("split_manifest", "locked_split", "split"):
        blob = report.get(key) or {}
        if isinstance(blob, dict) and blob.get("content_sha256"):
            return str(blob["content_sha256"])
    return None


def build_locked_truth_leaderboard(
    *,
    gcms_dir: Path | None = None,
) -> dict[str, Any]:
    """Build leaderboard rows from archived patient-diagnostic reports."""
    roots = [gcms_dir] if gcms_dir else _knowledge_roots()
    rows: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str, str]] = set()

    for root in roots:
        if not root or not root.exists():
            continue
        # Prefer explicit hybrid/literature subdirs + base study dirs
        report_paths = sorted(root.glob("**/PATIENT_DIAGNOSTIC_REPORT.json"))
        for path in report_paths:
            report = _load_json(path)
            if not report:
                continue
            study_id = str(report.get("study_id") or path.parent.name)
            signature = str(report.get("signature_source") or "unknown")
            method = str(report.get("score_method") or "cosine")
            key = (study_id, signature, method)
            if key in seen_keys:
                continue
            seen_keys.add(key)
            metrics = report.get("metrics") or {}
            nested = report.get("nested") or {}
            nested_auroc = nested.get("mean_test_auroc")
            if nested_auroc is None:
                nested_auroc = metrics.get("nested_auroc")
            non_nested = nested.get("non_nested_auroc")
            if non_nested is None:
                non_nested = metrics.get("auroc")
            optimism = nested.get("optimism_gap")
            if optimism is None and nested_auroc is not None and non_nested is not None:
                try:
                    optimism = float(non_nested) - float(nested_auroc)
                except (TypeError, ValueError):
                    optimism = None
            rows.append(
                {
                    "study_id": study_id,
                    "disease_id": report.get("disease_id"),
                    "signature": signature,
                    "method": method,
                    "n_subjects": metrics.get("n_subjects")
                    or (report.get("matrix_metadata") or {}).get("n_subjects"),
                    "nested_auroc": nested_auroc,
                    "non_nested_auroc": non_nested,
                    "optimism_gap": optimism,
                    "auprc": metrics.get("auprc"),
                    "sensitivity": metrics.get("sensitivity"),
                    "specificity": metrics.get("specificity"),
                    "split_sha256": _split_sha_from_report(report),
                    "source_path": str(path.relative_to(root)) if root in path.parents else str(path),
                    "clinical_claim": False,
                }
            )

        # Enrich from SIGNATURE_BENCHMARK if present and study missing rows
        bench = _load_json(root / "SIGNATURE_BENCHMARK.json")
        if bench:
            for r in bench.get("rows") or []:
                key = (
                    str(bench.get("study_id") or "ST000883"),
                    str(r.get("signature") or "unknown"),
                    str(r.get("method") or "cosine"),
                )
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                nested_auroc = r.get("nested_auroc")
                non_nested = r.get("auroc")
                gap = None
                if nested_auroc is not None and non_nested is not None:
                    gap = float(non_nested) - float(nested_auroc)
                rows.append(
                    {
                        "study_id": key[0],
                        "disease_id": None,
                        "signature": key[1],
                        "method": key[2],
                        "n_subjects": None,
                        "nested_auroc": nested_auroc,
                        "non_nested_auroc": non_nested,
                        "optimism_gap": gap,
                        "auprc": r.get("auprc"),
                        "sensitivity": None,
                        "specificity": None,
                        "split_sha256": None,
                        "source_path": "SIGNATURE_BENCHMARK.json",
                        "clinical_claim": False,
                        "logistic_nested_ceiling": r.get("logistic_nested"),
                    }
                )

    def _sort_key(row: dict[str, Any]) -> tuple:
        n = row.get("nested_auroc")
        return (
            0 if n is None else 1,
            float(n) if n is not None else -1.0,
            str(row.get("study_id") or ""),
            str(row.get("signature") or ""),
        )

    rows.sort(key=_sort_key, reverse=True)
    return {
        "schema_version": "LockedTruthLeaderboard-1.0",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "honesty": (
            "Patient-level nested AUROC leaderboard for research diligence only. "
            "Not clinical diagnostic SOTA. Prefer nested_auroc + optimism_gap + "
            "split_sha256 over directional literature percentages."
        ),
        "n_rows": len(rows),
        "rows": rows,
        "ranking_rule": "nested_auroc descending (missing nested last)",
    }


def render_leaderboard_md(report: dict[str, Any]) -> str:
    lines = [
        "# Locked patient-level truth leaderboard",
        "",
        f"Generated: {report.get('generated_utc')}",
        "",
        f"> {report.get('honesty')}",
        "",
        "| Study | Disease | Signature | Method | Nested AUROC | Optimism gap | Split SHA256 | n |",
        "|---|---|---|---|---:|---:|---|---:|",
    ]
    for r in report.get("rows") or []:
        sha = r.get("split_sha256") or "—"
        if isinstance(sha, str) and len(sha) > 12:
            sha = sha[:12] + "…"
        nested = r.get("nested_auroc")
        gap = r.get("optimism_gap")
        nested_s = "—" if nested is None else f"{100 * float(nested):.1f}%"
        gap_s = "—" if gap is None else f"{100 * float(gap):.1f}%"
        lines.append(
            f"| {r.get('study_id')} | {r.get('disease_id') or '—'} | "
            f"{r.get('signature')} | {r.get('method')} | {nested_s} | {gap_s} | "
            f"`{sha}` | {r.get('n_subjects') or '—'} |"
        )
    lines += [
        "",
        "## How to regenerate",
        "",
        "```bash",
        "voc eval-leaderboard",
        "voc eval-patient-diagnostic --all --signature hybrid",
        "```",
        "",
    ]
    return "\n".join(lines)


def write_locked_truth_leaderboard(
    *,
    out_dir: Path | None = None,
    gcms_dir: Path | None = None,
) -> dict[str, Any]:
    report = build_locked_truth_leaderboard(gcms_dir=gcms_dir)
    out = Path(out_dir) if out_dir else (DATA_DIR / "knowledge" / "leaderboard")
    out.mkdir(parents=True, exist_ok=True)
    (out / "LOCKED_TRUTH_LEADERBOARD.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "LOCKED_TRUTH_LEADERBOARD.md").write_text(render_leaderboard_md(report))
    # also mirror under package knowledge when writing to DATA_DIR
    pkg = PACKAGE_ROOT / "data" / "knowledge" / "leaderboard"
    if out.resolve() != pkg.resolve() and (PACKAGE_ROOT / "data" / "knowledge").exists():
        pkg.mkdir(parents=True, exist_ok=True)
        (pkg / "LOCKED_TRUTH_LEADERBOARD.json").write_text(json.dumps(report, indent=2) + "\n")
        (pkg / "LOCKED_TRUTH_LEADERBOARD.md").write_text(render_leaderboard_md(report))
    return report


__all__ = [
    "build_locked_truth_leaderboard",
    "write_locked_truth_leaderboard",
    "render_leaderboard_md",
]
