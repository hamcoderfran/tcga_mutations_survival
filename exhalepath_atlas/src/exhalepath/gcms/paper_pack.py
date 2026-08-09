"""One-zip research paper pack: figures + Methods + overlay + split hash."""

from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional


_DEFAULT_INCLUDE = (
    "PATIENT_DIAGNOSTIC_REPORT.md",
    "PATIENT_DIAGNOSTIC_REPORT.json",
    "TRIPOD_AI_BREATHVOC_CHECKLIST.md",
    "roc_curve.png",
    "patient_scores.csv",
    "METHODS.md",
    "literature_overlay.csv",
    "literature_overlay.json",
    "NEXT_EXPERIMENTS.md",
    "graphpad_voc_long.csv",
    "filter_report.json",
    "stratified_auroc.json",
    "SCIDATA_SAMPLE_EVAL.md",
    "SCIDATA_SAMPLE_EVAL.json",
    "README_METABOLIGHTS.md",
)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_methods_stub(
    out_path: Path,
    *,
    study_id: str,
    disease_id: str,
    n_subjects: int | None,
    split_sha256: str | None,
    signature_source: str | None = None,
    extra_notes: Optional[list[str]] = None,
) -> Path:
    """Minimal Methods.md for paper packs (research enablement)."""
    out_path = Path(out_path)
    notes = extra_notes or []
    lines = [
        "# Methods (auto-draft)",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Study",
        "",
        f"- Study / deposit id: `{study_id}`",
        f"- Disease / positive class: `{disease_id}`",
        f"- n subjects (analyzed): {n_subjects if n_subjects is not None else '—'}",
        f"- Signature source (if used): `{signature_source or 'n/a'}`",
        "",
        "## Splits",
        "",
        f"- Locked split content SHA256: `{split_sha256 or 'not recorded'}`",
        "- Train-only control means / scalers; no test leakage by construction.",
        "",
        "## Metrics",
        "",
        "- Primary: AUROC with bootstrap CI95; AUPRC; Youden sens/spec.",
        "- Nested vs non-nested optimism gap reported when splits available.",
        "- Stratified AUROC by age tertile / sex / smoking when metadata exist.",
        "",
        "## Preprocessing",
        "",
        "- Blank-ratio / detection-fraction feature filters when blanks or sparse peaks present.",
        "- Atlas VOC name mapping; unmapped peaks may be retained as `cid:<pubchem>` for ML baselines.",
        "",
        "## Honesty",
        "",
        "- Research enablement / hypothesis generation — **not** a clinical validation claim.",
        "",
    ]
    if notes:
        lines += ["## Notes", ""]
        for n in notes:
            lines.append(f"- {n}")
        lines.append("")
    out_path.write_text("\n".join(lines))
    return out_path


def collect_paper_pack_files(
    run_dir: Path,
    *,
    extra_globs: Iterable[str] = (),
) -> list[Path]:
    """Gather known report artifacts under a run directory."""
    run_dir = Path(run_dir)
    found: list[Path] = []
    names = set(_DEFAULT_INCLUDE)
    for p in run_dir.rglob("*"):
        if not p.is_file():
            continue
        if p.name in names or p.name.endswith("_split_manifest.json"):
            found.append(p)
            continue
        if p.suffix.lower() in {".png", ".pdf"} and "roc" in p.name.lower():
            found.append(p)
    for pattern in extra_globs:
        found.extend(run_dir.glob(pattern))
        found.extend(run_dir.rglob(pattern))
    # de-dupe preserving order
    seen: set[Path] = set()
    out: list[Path] = []
    for p in found:
        rp = p.resolve()
        if rp in seen:
            continue
        seen.add(rp)
        out.append(p)
    return out


def export_paper_pack(
    run_dir: Path,
    out_zip: Path | None = None,
    *,
    study_id: str | None = None,
    ensure_methods: bool = True,
    split_sha256: str | None = None,
    disease_id: str | None = None,
    n_subjects: int | None = None,
    signature_source: str | None = None,
    methods_notes: Optional[list[str]] = None,
) -> dict[str, Any]:
    """Zip figures + Methods + overlay + split hash into one paper pack.

    Returns dict with zip path, sha256 of zip, included relative paths, split hash.
    """
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(run_dir)

    # pull split hash from manifest if present
    if split_sha256 is None:
        for man in run_dir.rglob("*_split_manifest.json"):
            try:
                payload = json.loads(man.read_text())
                split_sha256 = payload.get("content_sha256") or split_sha256
                study_id = study_id or payload.get("dataset_id")
            except Exception:  # noqa: BLE001
                pass

    if ensure_methods and not (run_dir / "METHODS.md").exists():
        write_methods_stub(
            run_dir / "METHODS.md",
            study_id=study_id or run_dir.name,
            disease_id=disease_id or "unknown",
            n_subjects=n_subjects,
            split_sha256=split_sha256,
            signature_source=signature_source,
            extra_notes=methods_notes,
        )

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "study_id": study_id or run_dir.name,
        "split_content_sha256": split_sha256,
        "disease_id": disease_id,
        "n_subjects": n_subjects,
        "files": [],
    }
    files = collect_paper_pack_files(run_dir)
    out_zip = Path(out_zip or (run_dir / f"{manifest['study_id']}_paper_pack.zip"))
    out_zip.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            rel = path.relative_to(run_dir).as_posix()
            zf.write(path, arcname=rel)
            manifest["files"].append(
                {"path": rel, "sha256": _sha256_file(path), "bytes": path.stat().st_size}
            )
        zf.writestr("PAPER_PACK_MANIFEST.json", json.dumps(manifest, indent=2))

    return {
        "zip_path": str(out_zip),
        "zip_sha256": _sha256_file(out_zip),
        "n_files": len(manifest["files"]),
        "split_content_sha256": split_sha256,
        "manifest": manifest,
    }


def copy_overlay_into_run(run_dir: Path, overlay_paths: dict[str, Path]) -> None:
    """Copy literature overlay / GraphPad artifacts into a diagnostic run dir."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    for name, src in overlay_paths.items():
        src = Path(src)
        if src.is_file():
            shutil.copy2(src, run_dir / (name if name.endswith((".csv", ".json", ".md")) else src.name))


__all__ = [
    "collect_paper_pack_files",
    "export_paper_pack",
    "write_methods_stub",
]
