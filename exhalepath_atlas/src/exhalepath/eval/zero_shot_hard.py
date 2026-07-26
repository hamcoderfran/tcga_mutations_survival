"""Adversarial hard-break eval for zero-shot disease prediction."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from ..biomarker import ExhaleBiomarkerEngine
from ..config import KNOWLEDGE_DIR
from .zero_shot_reliability import _pathway_map, _pred_map


def _load_hard(path: Path | None = None) -> dict[str, Any]:
    candidates = [
        path,
        KNOWLEDGE_DIR / "zero_shot_hard_break.json",
        Path("data/knowledge/zero_shot_hard_break.json"),
    ]
    for p in candidates:
        if p and Path(p).exists():
            return json.loads(Path(p).read_text())
    raise FileNotFoundError("zero_shot_hard_break.json not found")


def evaluate_zero_shot_hard(
    *,
    out_dir: Path | None = None,
    profiles_path: Path | None = None,
    mode: str = "hybrid",
) -> dict[str, Any]:
    doc = _load_hard(profiles_path)
    profiles = list(doc.get("profiles") or [])
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    rows: list[dict[str, Any]] = []
    n_pass = 0
    trap_counts: Counter[str] = Counter()
    trap_fail: Counter[str] = Counter()

    for prof in profiles:
        trap = str(prof.get("trap") or "other")
        trap_counts[trap] += 1
        genes = [str(g).upper() for g in (prof.get("genes") or [])]
        report = engine.predict(
            prof["disease"],
            location=prof.get("location"),
            genes=genes or None,
            description=prof.get("description"),
            pathway_overrides=prof.get("pathway_overrides"),
            comorbidities=prof.get("comorbidities"),
            mode=mode,
            explain=False,
        )
        meta = report.result.bundle.metadata or {}
        pred = _pred_map(report)
        pw = _pathway_map(report)
        peak = max((abs(v) for v in pred.values()), default=0.0)
        checks: list[dict[str, Any]] = []

        def add(name: str, ok: bool, detail: Any = None) -> None:
            checks.append({"name": name, "pass": bool(ok), "detail": detail})

        unresolved = bool(
            meta.get("unresolved") or str(report.disease_id).startswith("custom::")
        )
        if prof.get("expect_unresolved") is not None:
            add(
                "unresolved",
                unresolved == bool(prof["expect_unresolved"]),
                report.disease_id,
            )

        mech_c = float(meta.get("mechanism_confidence") or 0.0)
        zs_mode = meta.get("zero_shot_mode")
        if prof.get("expect_mechanism") is True:
            add(
                "mechanism_mode",
                (zs_mode == "mechanism")
                or (mech_c >= float(prof.get("min_mechanism_confidence") or 0.35))
                or (peak >= float(prof.get("min_peak_abs_log2fc") or 0.2)),
                {"mode": zs_mode, "mechanism_confidence": mech_c, "peak": peak},
            )
        if prof.get("expect_mechanism") is False:
            add(
                "near_healthy_fallback",
                zs_mode in {None, "near_healthy_fallback"}
                or peak <= float(prof.get("max_abs_log2fc") or 0.5),
                {"mode": zs_mode, "peak": peak},
            )

        if prof.get("near_healthy"):
            thr = float(prof.get("max_abs_log2fc") or 0.45)
            add("near_healthy", peak <= thr, {"peak": peak, "threshold": thr})

        if prof.get("min_peak_abs_log2fc") is not None:
            add(
                "min_peak_abs_log2fc",
                peak >= float(prof["min_peak_abs_log2fc"]),
                peak,
            )

        for voc in prof.get("must_elevate") or []:
            add(f"elevate:{voc}", float(pred.get(voc, 0.0)) > 0.05, pred.get(voc))
        for voc in prof.get("must_suppress") or []:
            add(f"suppress:{voc}", float(pred.get(voc, 0.0)) < -0.05, pred.get(voc))

        if prof.get("pathway_any"):
            hit = any(float(pw.get(p, 0.0)) > 0.05 for p in prof["pathway_any"])
            add(
                "pathway_any",
                hit,
                {p: pw.get(p) for p in prof["pathway_any"]},
            )

        if unresolved:
            d = engine.kb.resolve_disease(
                prof["disease"],
                location_hint=prof.get("location"),
                description=prof.get("description"),
            )
            add(
                "no_voc_prior_leak",
                not bool(d.get("voc_log2fc_prior")),
                d.get("voc_log2fc_prior"),
            )
            # Prefer disease-text site over forced wrong location when available
            prefer = prof.get("site_should_prefer")
            if prefer:
                got_site = str(d.get("default_site") or "").lower()
                add(
                    "site_prefer",
                    prefer.lower() in got_site or prefer.lower() in str(report.location.get("tissue_id") or "").lower(),
                    {"default_site": d.get("default_site"), "location": report.location},
                )

        passed = all(c["pass"] for c in checks) if checks else False
        if passed:
            n_pass += 1
        else:
            trap_fail[trap] += 1

        top = sorted(pred.items(), key=lambda x: abs(x[1]), reverse=True)[:8]
        failed = [c["name"] for c in checks if not c["pass"]]
        rows.append(
            {
                "id": prof.get("id"),
                "trap": trap,
                "difficulty": prof.get("difficulty"),
                "disease": prof["disease"],
                "disease_id": report.disease_id,
                "location": report.resolved_location,
                "passed": passed,
                "failed_checks": failed,
                "peak_abs_log2fc": peak,
                "mechanism_confidence": mech_c,
                "zero_shot_mode": zs_mode,
                "category": meta.get("category"),
                "default_site": meta.get("default_site"),
                "top_vocs": [{"voc_id": v, "log2fc": fc} for v, fc in top],
                "checks": checks,
                "notes": prof.get("notes"),
            }
        )

    n = len(rows)
    by_trap = {
        t: {
            "n": trap_counts[t],
            "failed": trap_fail[t],
            "pass_rate": (trap_counts[t] - trap_fail[t]) / trap_counts[t]
            if trap_counts[t]
            else 0.0,
        }
        for t in sorted(trap_counts)
    }
    report = {
        "n_profiles": n,
        "n_passed": n_pass,
        "n_failed": n - n_pass,
        "pass_rate": n_pass / n if n else 0.0,
        # Hard suite is meant to expose breaks — "pass" if most cases still work
        "pass": (n_pass / n if n else 0.0) >= 0.55,
        "by_trap": by_trap,
        "failed_ids": [r["id"] for r in rows if not r["passed"]],
        "profiles": rows,
    }

    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "zero_shot_hard_break.json").write_text(json.dumps(report, indent=2))
        lines = [
            "# Zero-shot hard-break evaluation",
            "",
            f"**Pass rate:** {100 * report['pass_rate']:.0f}% "
            f"({report['n_passed']}/{report['n_profiles']})",
            f"**Suite gate (≥55%):** {report['pass']}",
            "",
            "## By trap",
        ]
        for t, row in by_trap.items():
            lines.append(
                f"- `{t}`: {row['n'] - row['failed']}/{row['n']} "
                f"({100 * row['pass_rate']:.0f}%)"
            )
        lines.extend(["", "## Failures"])
        fails = [r for r in rows if not r["passed"]]
        if not fails:
            lines.append("- none")
        for r in fails:
            lines.append(
                f"- ✗ **{r['id']}** `{r['disease']}` trap=`{r['trap']}` "
                f"peak={r['peak_abs_log2fc']:.3f} mech={r['mechanism_confidence']:.2f} "
                f"mode={r['zero_shot_mode']} failed={r['failed_checks']}"
            )
            lines.append(
                "  top: "
                + ", ".join(f"{x['voc_id']}({x['log2fc']:+.2f})" for x in r["top_vocs"][:5])
            )
        lines.extend(["", "## Passes (compact)"])
        for r in rows:
            if not r["passed"]:
                continue
            lines.append(
                f"- ✓ {r['id']} peak={r['peak_abs_log2fc']:.3f} "
                f"mech={r['mechanism_confidence']:.2f}"
            )
        (out / "ZERO_SHOT_HARD_BREAK.md").write_text("\n".join(lines) + "\n")
    return report
