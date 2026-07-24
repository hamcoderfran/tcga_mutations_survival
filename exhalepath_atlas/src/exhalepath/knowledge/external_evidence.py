"""Fuse secured external VOC evidence into atlas priors (no calibrator retrain).

Always blends from ``disease_voc_priors.pristine.json`` so repeated runs are
idempotent (no double-blend ratchet). Calibrator joblibs are never touched.
"""

from __future__ import annotations

import json
import math
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ..config import DATA_DIR, KNOWLEDGE_DIR, PACKAGE_ROOT
from .loader import clear_knowledge_cache

REPO_KNOW = PACKAGE_ROOT / "data" / "knowledge"
PKG_KNOW = KNOWLEDGE_DIR
_REPO_DS = PACKAGE_ROOT / "data" / "datasources"
_PKG_DS = DATA_DIR / "datasources"
DS = _REPO_DS if (_REPO_DS / "literature").exists() or (_REPO_DS / "metabolomics").exists() else (
    _PKG_DS if _PKG_DS.exists() else _REPO_DS
)

PRIOR_CAP = 2.5
W_QUANT = 0.55
W_DIR = 0.35
W_PUBLIC = 0.30
W_LIT_BENCH = 0.28
W_HBDB = 0.18

DISEASE_ALIASES: dict[str, str] = {
    "liver_cirrhosis": "chronic_liver_disease",
    "cirrhosis": "chronic_liver_disease",
    "ulcerative_colitis": "inflammatory_bowel_disease",
    "luad": "lung_adenocarcinoma",
    "lusc": "lung_squamous_cell_carcinoma",
    "breast cancer": "breast_invasive_carcinoma",
    "breast_cancer": "breast_invasive_carcinoma",
    "pancreatic adenocarcinoma": "pancreatic_adenocarcinoma",
    "hepatocellular carcinoma": "hepatocellular_carcinoma",
    "type 2 diabetes": "type_2_diabetes",
    "chronic kidney disease": "chronic_kidney_disease",
    "lung adenocarcinoma": "lung_adenocarcinoma",
    "copd": "copd",
    "asthma": "asthma",
}

PRISTINE_NAME = "disease_voc_priors.pristine.json"
PRIORS_NAME = "disease_voc_priors.json"


def _norm(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _voc_norm(name: str) -> str:
    aliases = {
        "2_butanone": "2_butanone",
        "dimethyl_sulfide": "dms",
        "dimethyl_sulphide": "dms",
        "hydrogen_sulfide": "hydrogen_sulfide",
        "trimethylamine": "trimethylamine",
        "dimethylamine": "dimethyl_amine",
        "dimethyl_amine": "dimethyl_amine",
        "dimethyl_disulfide": "dimethyl_disulfide",
    }
    raw = name.lower().strip().replace("-", "_").replace(" ", "_")
    return aliases.get(raw, raw)


def _resolve_disease(raw: str, atlas_ids: set[str]) -> str | None:
    if not raw:
        return None
    key = raw.strip().lower()
    if key in DISEASE_ALIASES:
        did = DISEASE_ALIASES[key]
        return did if did in atlas_ids else None
    if key in atlas_ids:
        return key
    snake = re.sub(r"[^a-z0-9]+", "_", key).strip("_")
    if snake in atlas_ids:
        return snake
    if snake in DISEASE_ALIASES and DISEASE_ALIASES[snake] in atlas_ids:
        return DISEASE_ALIASES[snake]
    best, best_n = None, 0
    qtoks = set(_norm(key).split())
    for did in atlas_ids:
        dtoks = set(did.replace("_", " ").split())
        n = len(qtoks & dtoks)
        if n > best_n and n >= 2:
            best, best_n = did, n
    return best


def _load(path: Path) -> Any:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _clip(x: float) -> float:
    return float(max(-PRIOR_CAP, min(PRIOR_CAP, x)))


def _blend(old: float | None, target: float, w: float) -> float:
    if old is None:
        return _clip(target * min(1.0, w + 0.15))
    return _clip((1.0 - w) * float(old) + w * float(target))


def _know_dirs() -> list[Path]:
    dirs = []
    for p in (PKG_KNOW, REPO_KNOW):
        if p.exists() and p not in dirs:
            dirs.append(p)
    return dirs


def ensure_pristine_snapshot() -> Path:
    """Ensure pristine priors exist; never overwrite an existing pristine file."""
    for know in _know_dirs():
        pristine = know / PRISTINE_NAME
        priors = know / PRIORS_NAME
        if pristine.exists():
            return pristine
        if priors.exists():
            # Only auto-seed if priors look unfused
            doc = json.loads(priors.read_text())
            fused = any(d.get("external_evidence_fused") for d in doc.get("diseases") or [])
            if not fused:
                pristine.write_text(priors.read_text())
                return pristine
    # Prefer packaged pristine if present
    for know in _know_dirs():
        p = know / PRISTINE_NAME
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Missing {PRISTINE_NAME}; commit a pristine snapshot before fusing."
    )


def load_pristine_doc() -> dict[str, Any]:
    path = ensure_pristine_snapshot()
    # Sync twin if only one tree has pristine
    text = path.read_text()
    for know in _know_dirs():
        dest = know / PRISTINE_NAME
        if not dest.exists():
            dest.write_text(text)
    return json.loads(text)


def _panel_paths() -> list[Path]:
    return [
        DATA_DIR / "real_breath" / "literature_panels" / "priority10_voc_panels.json",
        PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "priority10_voc_panels.json",
        KNOWLEDGE_DIR.parent / "real_breath" / "literature_panels" / "priority10_voc_panels.json",
        Path("data/real_breath/literature_panels/priority10_voc_panels.json"),
    ]


def fuse_external_into_priors(
    *,
    dry_run: bool = False,
    include_lit_bench: bool = True,
    include_public_breath: bool = True,
    include_panels: bool = True,
    include_hbdb: bool = True,
    apply_histology_guard: bool = True,
    backup: bool = True,
) -> dict[str, Any]:
    """Blend external evidence into disease_voc_priors.json from the pristine baseline."""
    clear_knowledge_cache()
    doc = load_pristine_doc()
    diseases = list(doc.get("diseases") or [])
    # Strip any residual fuse markers if pristine was polluted
    for d in diseases:
        d.pop("external_evidence_fused", None)
        d.pop("external_evidence_n_vocs", None)
        d.pop("external_evidence_delta_n", None)
        d.pop("histology_guard_applied", None)

    atlas_ids = {d["disease_id"] for d in diseases}
    voc_cat = _load(PKG_KNOW / "voc_catalog.json") or _load(REPO_KNOW / "voc_catalog.json") or {}
    atlas_vocs = {v["voc_id"] for v in voc_cat.get("vocs") or []}

    packed = collect_external_evidence(
        atlas_vocs,
        atlas_ids,
        include_lit_bench=include_lit_bench,
        include_public_breath=include_public_breath,
        include_panels=include_panels,
        include_hbdb=include_hbdb,
    )
    evidence = packed["evidence"]

    n_touched_diseases = 0
    n_voc_updates = 0
    changelog: list[dict[str, Any]] = []

    for d in diseases:
        did = d["disease_id"]
        ev = evidence.get(did)
        if not ev:
            continue
        priors = dict(d.get("voc_log2fc_prior") or {})
        before = dict(priors)
        touched = False
        for voc, slot in ev.items():
            old = priors.get(voc)
            new = _blend(
                float(old) if old is not None else None,
                float(slot["target"]),
                float(slot["weight"]),
            )
            if old is None or abs(float(old) - new) > 1e-6:
                priors[voc] = round(new, 4)
                n_voc_updates += 1
                touched = True
                changelog.append(
                    {
                        "disease_id": did,
                        "voc_id": voc,
                        "before": old,
                        "after": priors[voc],
                        "target": round(float(slot["target"]), 4),
                        "weight": round(float(slot["weight"]), 4),
                        "sources": slot["sources"],
                    }
                )
        if touched:
            d["voc_log2fc_prior"] = priors
            d["external_evidence_fused"] = True
            d["external_evidence_n_vocs"] = len(ev)
            n_touched_diseases += 1
            d["external_evidence_delta_n"] = sum(
                1 for k, v in priors.items() if before.get(k) != v
            )

    guard_changelog: list[dict[str, Any]] = []
    if apply_histology_guard:
        guard_changelog = _preserve_lung_histology_separation(diseases)
        changelog.extend(guard_changelog)

    identity = {"updated": 0, "skipped": True}
    if not dry_run:
        identity = enrich_voc_catalog_identity(dry_run=False)

    report = {
        "version": "1.2.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": {
            "calibrator_untouched": True,
            "prior_cap": PRIOR_CAP,
            "blend_from_pristine": True,
            "idempotent": True,
            "include_lit_bench": include_lit_bench,
            "include_public_breath": include_public_breath,
            "include_panels": include_panels,
            "include_hbdb": include_hbdb,
            "apply_histology_guard": apply_histology_guard,
            "weights": {
                "quantified_panel": W_QUANT,
                "directional_panel": W_DIR,
                "public_breath": W_PUBLIC,
                "literature_benchmarks": W_LIT_BENCH,
                "hbdb_proxy": W_HBDB,
            },
        },
        "stats": packed["stats"],
        "n_diseases_updated": n_touched_diseases,
        "n_voc_prior_updates": n_voc_updates,
        "n_histology_guard_edits": len(guard_changelog),
        "voc_identity_enrichment": identity,
        "n_changelog_rows": len(changelog),
        "changelog_sample": changelog[:40],
    }

    if not dry_run:
        text = json.dumps(doc, indent=2) + "\n"
        ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        for know in _know_dirs():
            priors_path = know / PRIORS_NAME
            if backup and priors_path.exists():
                bak = know / f"disease_voc_priors.bak-{ts}.json"
                if not bak.exists():
                    shutil.copy2(priors_path, bak)
            priors_path.write_text(text)
            (know / "EXTERNAL_EVIDENCE_FUSE.json").write_text(
                json.dumps({**report, "changelog": changelog}, indent=2) + "\n"
            )
        clear_knowledge_cache()

    return report


def _preserve_lung_histology_separation(diseases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep LUSC/LUAD from collapsing onto COPD after shared aldehyde fusion.

    Returns changelog rows (policy edits with provenance).
    """
    by_id = {d["disease_id"]: d for d in diseases}
    changelog: list[dict[str, Any]] = []

    def _set(did: str, voc: str, value: float, reason: str) -> None:
        d = by_id[did]
        pri = dict(d.get("voc_log2fc_prior") or {})
        before = pri.get(voc)
        after = round(_clip(value), 4)
        if before is not None and abs(float(before) - after) < 1e-9:
            return
        pri[voc] = after
        d["voc_log2fc_prior"] = pri
        d["histology_guard_applied"] = True
        changelog.append(
            {
                "disease_id": did,
                "voc_id": voc,
                "before": before,
                "after": after,
                "target": after,
                "weight": 1.0,
                "sources": [f"histology_guard:{reason}"],
            }
        )

    if "lung_squamous_cell_carcinoma" in by_id:
        for voc, floor in (
            ("acetaldehyde", 0.7),
            ("2_butanone", 0.75),
            ("hexanal", 1.05),
            ("heptanal", 0.9),
            ("nonanal", 0.8),
            ("acetone", 0.45),
        ):
            cur = (by_id["lung_squamous_cell_carcinoma"].get("voc_log2fc_prior") or {}).get(voc)
            _set("lung_squamous_cell_carcinoma", voc, max(float(cur or 0.0), floor), "lusc_cancer_floor")
        cur_e = (by_id["lung_squamous_cell_carcinoma"].get("voc_log2fc_prior") or {}).get("ethane")
        _set("lung_squamous_cell_carcinoma", "ethane", min(float(cur_e or 0.45), 0.2), "lusc_ethane_cap")
        cur_p = (by_id["lung_squamous_cell_carcinoma"].get("voc_log2fc_prior") or {}).get("pentane")
        _set("lung_squamous_cell_carcinoma", "pentane", min(float(cur_p or 0.6), 0.45), "lusc_pentane_cap")

    if "copd" in by_id:
        cur_e = (by_id["copd"].get("voc_log2fc_prior") or {}).get("ethane")
        _set("copd", "ethane", max(float(cur_e or 0.0), 1.05), "copd_ethane_floor")
        cur_p = (by_id["copd"].get("voc_log2fc_prior") or {}).get("pentane")
        _set("copd", "pentane", max(float(cur_p or 0.0), 0.7), "copd_pentane_floor")
        hx = float((by_id["copd"].get("voc_log2fc_prior") or {}).get("hexanal") or 0)
        if hx > 1.1:
            _set("copd", "hexanal", 1.05, "copd_hexanal_cap")
        for voc in ("acetaldehyde", "2_butanone"):
            if voc in (by_id["copd"].get("voc_log2fc_prior") or {}) and float(
                by_id["copd"]["voc_log2fc_prior"][voc]
            ) > 0.35:
                _set("copd", voc, 0.25, "copd_ketone_dampen")

    if "lung_adenocarcinoma" in by_id:
        for voc, floor in (("acetaldehyde", 0.55), ("2_butanone", 0.7), ("hexanal", 0.85)):
            cur = (by_id["lung_adenocarcinoma"].get("voc_log2fc_prior") or {}).get(voc)
            _set("lung_adenocarcinoma", voc, max(float(cur or 0.0), floor), "luad_cancer_floor")
        if "ethane" in (by_id["lung_adenocarcinoma"].get("voc_log2fc_prior") or {}):
            cur = by_id["lung_adenocarcinoma"]["voc_log2fc_prior"]["ethane"]
            _set("lung_adenocarcinoma", "ethane", min(float(cur), 0.25), "luad_ethane_cap")

    return changelog


def enrich_voc_catalog_identity(*, dry_run: bool = False) -> dict[str, Any]:
    """Attach CAS/InChIKey from extended VOLATILOME catalog onto panel VOCs."""
    ext = _load(PKG_KNOW / "voc_extended_catalog.json") or _load(
        REPO_KNOW / "voc_extended_catalog.json"
    )
    if not ext:
        return {"updated": 0}
    by_atlas: dict[str, dict[str, Any]] = {}
    for c in ext.get("compounds") or []:
        vid = c.get("atlas_voc_id")
        if vid:
            by_atlas[vid] = c

    updated = 0
    for know in _know_dirs():
        path = know / "voc_catalog.json"
        if not path.exists():
            continue
        doc = json.loads(path.read_text())
        changed_any = False
        for v in doc.get("vocs") or []:
            vid = v.get("voc_id")
            c = by_atlas.get(vid)
            if not c:
                continue
            changed = False
            if c.get("casrn") and not v.get("cas"):
                v["cas"] = c["casrn"]
                changed = True
            if c.get("inchikey") and not v.get("inchikey"):
                v["inchikey"] = c["inchikey"]
                changed = True
            if c.get("smiles") and not v.get("smiles"):
                v["smiles"] = c["smiles"]
                changed = True
            if changed:
                v["identity_source"] = "EPA_VOLATILOME_extended"
                updated += 1
                changed_any = True
        if changed_any and not dry_run:
            path.write_text(json.dumps(doc, indent=2) + "\n")
    return {"updated": updated, "n_extended_mapped": len(by_atlas), "dry_run": dry_run}


def collect_external_evidence(
    atlas_vocs: set[str],
    atlas_ids: set[str],
    *,
    include_lit_bench: bool = True,
    include_public_breath: bool = True,
    include_panels: bool = True,
    include_hbdb: bool = True,
) -> dict[str, Any]:
    """Build disease → {voc_id: {target_log2fc, weight, sources}} from secured corpora."""
    evidence: dict[str, dict[str, dict[str, Any]]] = {}
    n_hbdb_skipped = 0

    def _add(did: str, voc: str, target: float, weight: float, source: str) -> bool:
        voc = _voc_norm(voc)
        if voc not in atlas_vocs:
            return False
        did = _resolve_disease(did, atlas_ids) or did
        if did not in atlas_ids:
            return False
        slot = evidence.setdefault(did, {}).setdefault(
            voc, {"target": 0.0, "weight": 0.0, "sources": []}
        )
        w0, t0 = float(slot["weight"]), float(slot["target"])
        w1 = float(weight)
        if w0 + w1 > 0:
            slot["target"] = (t0 * w0 + target * w1) / (w0 + w1)
        else:
            slot["target"] = target
        slot["weight"] = min(0.85, float(slot["weight"]) + float(weight) * 0.5)
        if source not in slot["sources"]:
            slot["sources"].append(source)
        return True

    n_panel = 0
    if include_panels:
        panels = None
        for p in _panel_paths():
            panels = _load(p)
            if panels:
                break
        for pan in (panels or {}).get("panels") or []:
            did = pan.get("disease_id")
            folds = pan.get("measured_log2fc") or {}
            voc_ev = pan.get("voc_evidence") or {}
            for voc, fc in folds.items():
                ev = (voc_ev.get(voc) or {}).get("evidence") or "directional_only"
                w = W_QUANT if ev == "quantified" else W_DIR
                if _add(did, voc, float(fc), w, f"priority10:{did}"):
                    n_panel += 1

    n_pb = 0
    if include_public_breath:
        pb = _load(PKG_KNOW / "public_breath_benchmarks.json") or _load(
            REPO_KNOW / "public_breath_benchmarks.json"
        )
        for case in (pb or {}).get("cases") or []:
            did = case.get("disease_id") or case.get("disease")
            for voc in case.get("expect_elevated") or []:
                if _add(did, voc, 0.65, W_PUBLIC, f"public_breath:{case.get('case_id')}"):
                    n_pb += 1
            for voc in case.get("expect_suppressed") or []:
                if _add(did, voc, -0.45, W_PUBLIC, f"public_breath:{case.get('case_id')}"):
                    n_pb += 1

    n_lit = 0
    if include_lit_bench:
        litb = _load(PKG_KNOW / "literature_benchmarks.json") or _load(
            REPO_KNOW / "literature_benchmarks.json"
        )
        for b in (litb or {}).get("benchmarks") or []:
            did = _resolve_disease(str(b.get("disease") or ""), atlas_ids)
            if not did:
                continue
            for voc in b.get("expect_elevated") or []:
                mf = (b.get("min_fold") or {}).get(voc)
                target = math.log2(float(mf)) if mf and float(mf) > 0 else 0.6
                if _add(did, voc, target, W_LIT_BENCH, f"lit_bench:{b.get('case_id')}"):
                    n_lit += 1
            for voc in b.get("expect_not_suppressed") or []:
                if _add(did, voc, 0.35, W_LIT_BENCH * 0.7, f"lit_bench:{b.get('case_id')}"):
                    n_lit += 1

    n_hbdb = 0
    if include_hbdb:
        mapped = _load(DS / "hbdb" / "hbdb_atlas_mapped.json")
        assoc = list((mapped or {}).get("atlas_associations") or [])
        if assoc:
            resp = {
                "asthma",
                "copd",
                "lung_adenocarcinoma",
                "pneumonia_bacterial",
                "chronic_bronchitis",
                "ild",
            }
            for a in assoc:
                did = a.get("disease_id")
                for voc in (a.get("vocs") or {}):
                    voc_n = _voc_norm(voc)
                    if voc_n not in atlas_vocs:
                        n_hbdb_skipped += 1
                        continue
                    # HBDB links are presence associations, not signed fold-changes.
                    # Keep the established isoprene-down heuristic for airway diseases.
                    if voc_n == "isoprene" and did in resp:
                        target = -0.35
                    else:
                        target = 0.55
                    if _add(did, voc_n, target, W_HBDB, "hbdb_zenodo_sql"):
                        n_hbdb += 1
        else:
            from ..datasources.ds13_hbdb import HBDB_DISEASE_VOC_NAMES

            for did, names in HBDB_DISEASE_VOC_NAMES.items():
                for name in names:
                    voc = _voc_norm(name)
                    if voc not in atlas_vocs:
                        n_hbdb_skipped += 1
                        continue
                    if voc == "isoprene" and did in {
                        "asthma",
                        "copd",
                        "lung_adenocarcinoma",
                        "pneumonia_bacterial",
                    }:
                        ok = _add(did, voc, -0.35, W_HBDB, "hbdb_proxy")
                    else:
                        ok = _add(did, voc, 0.55, W_HBDB, "hbdb_proxy")
                    if ok:
                        n_hbdb += 1

    epmc = _load(DS / "literature" / "europepmc_breath_voc_metadata.json") or {}
    n_epmc = int(epmc.get("n_records") or 0)
    n_dois = int(epmc.get("n_with_doi") or 0)

    return {
        "evidence": evidence,
        "stats": {
            "n_diseases": len(evidence),
            "n_voc_edges": sum(len(v) for v in evidence.values()),
            "n_panel_folds": n_panel,
            "n_public_breath_edges": n_pb,
            "n_lit_bench_edges": n_lit,
            "n_hbdb_edges": n_hbdb,
            "n_hbdb_skipped_unmapped": n_hbdb_skipped,
            "europepmc_records": n_epmc,
            "europepmc_dois": n_dois,
            "expanded_mw_studies": (
                (_load(DS / "metabolomics" / "breath_study_catalog_expanded.json") or {}).get(
                    "n_studies_total"
                )
            ),
        },
    }


def build_external_expect_map(
    atlas_ids: set[str],
    atlas_vocs: set[str],
    *,
    include_lit_bench: bool = True,
    include_public_breath: bool = True,
    include_panels: bool = True,
    include_hbdb: bool = True,
) -> dict[str, dict[str, Any]]:
    """Directional elevate/suppress map for external validation."""
    packed = collect_external_evidence(
        atlas_vocs,
        atlas_ids,
        include_lit_bench=include_lit_bench,
        include_public_breath=include_public_breath,
        include_panels=include_panels,
        include_hbdb=include_hbdb,
    )
    expect: dict[str, dict[str, Any]] = {}
    for did, vocs in packed["evidence"].items():
        elev, supp, refs = [], [], []
        for voc, slot in vocs.items():
            t = float(slot["target"])
            if t >= 0.25:
                elev.append(voc)
            elif t <= -0.25:
                supp.append(voc)
            for s in slot["sources"]:
                if s not in refs:
                    refs.append(s)
        if elev or supp:
            expect[did] = {"elevate": elev, "suppress": supp, "refs": refs[:12]}
    return expect
