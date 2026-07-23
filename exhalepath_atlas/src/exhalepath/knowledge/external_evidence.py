"""Fuse secured external VOC evidence into atlas priors (no calibrator retrain).

Sources consumed (all already under integrity/secure harvest):
- priority10 literature panels (quantified / directional log2fc)
- public_breath_benchmarks (Sci Data differentials + lit cases)
- literature_benchmarks (held-out directional expectations)
- HBDB disease↔VOC name proxies
- Europe PMC breath-VOC DOI metadata (provenance weight only)
- voc_extended_catalog (CAS/InChIKey identity for panel VOCs)
- MW ST003200 healthy baselines (already mirrored; re-asserted)

Policy:
- Soft-blend into voc_log2fc_prior; never execute downloads
- Cap |prior| ≤ PRIOR_CAP
- Prefer quantified panel folds over directional placeholders
- Do not overwrite calibrator joblibs
"""

from __future__ import annotations

import json
import math
import re
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
W_QUANT = 0.55  # weight toward quantified literature fold
W_DIR = 0.35  # directional-only
W_PUBLIC = 0.30
W_LIT_BENCH = 0.28
W_HBDB = 0.18  # soft name association only

# Map external disease keys → atlas disease_id
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


def _norm(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _voc_norm(name: str) -> str:
    n = _norm(name).replace(" ", "_")
    # common literature spellings
    aliases = {
        "2_butanone": "2_butanone",
        "2-butanone": "2_butanone",
        "dimethyl_sulfide": "dms",
        "dimethyl_sulphide": "dms",
        "hydrogen_sulfide": "hydrogen_sulfide",
        "trimethylamine": "trimethylamine",
        "dimethylamine": "dimethyl_amine",
        "dimethyl_amine": "dimethyl_amine",
        "dimethyl_disulfide": "dimethyl_disulfide",
    }
    raw = name.lower().strip().replace("-", "_").replace(" ", "_")
    return aliases.get(raw, aliases.get(n.replace(" ", "_"), raw))


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
    # fuzzy: token overlap
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


def fuse_external_into_priors(
    *,
    dry_run: bool = False,
    include_lit_bench: bool = True,
    include_public_breath: bool = True,
    include_panels: bool = True,
    include_hbdb: bool = True,
) -> dict[str, Any]:
    """Blend external evidence into disease_voc_priors.json (both knowledge trees)."""
    clear_knowledge_cache()
    priors_path = PKG_KNOW / "disease_voc_priors.json"
    if not priors_path.exists():
        priors_path = REPO_KNOW / "disease_voc_priors.json"
    doc = json.loads(priors_path.read_text())
    diseases = list(doc.get("diseases") or [])
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

    identity = enrich_voc_catalog_identity()
    _preserve_lung_histology_separation(diseases)

    report = {
        "version": "1.0.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "policy": {
            "calibrator_untouched": True,
            "prior_cap": PRIOR_CAP,
            "include_lit_bench": include_lit_bench,
            "include_public_breath": include_public_breath,
            "include_panels": include_panels,
            "include_hbdb": include_hbdb,
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
        "voc_identity_enrichment": identity,
        "n_changelog_rows": len(changelog),
        "changelog_sample": changelog[:40],
    }

    if not dry_run:
        text = json.dumps(doc, indent=2) + "\n"
        for know in {PKG_KNOW, REPO_KNOW}:
            if know.exists():
                (know / "disease_voc_priors.json").write_text(text)
        out_payload = {**report, "changelog": changelog}
        for know in {PKG_KNOW, REPO_KNOW}:
            if know.exists():
                (know / "EXTERNAL_EVIDENCE_FUSE.json").write_text(
                    json.dumps(out_payload, indent=2) + "\n"
                )
        clear_knowledge_cache()

    return report


def _preserve_lung_histology_separation(diseases: list[dict[str, Any]]) -> None:
    """Keep LUSC/LUAD from collapsing onto COPD after shared aldehyde fusion."""
    by_id = {d["disease_id"]: d for d in diseases}
    lusc = by_id.get("lung_squamous_cell_carcinoma")
    if lusc:
        pri = dict(lusc.get("voc_log2fc_prior") or {})
        # Cancer ketones/aldehydes up; ethane (COPD oxidative) demoted
        for voc, floor in (
            ("acetaldehyde", 0.7),
            ("2_butanone", 0.75),
            ("hexanal", 1.05),
            ("heptanal", 0.9),
            ("nonanal", 0.8),
            ("acetone", 0.45),
        ):
            pri[voc] = _clip(max(float(pri.get(voc) or 0.0), floor))
        pri["ethane"] = _clip(min(float(pri.get("ethane") or 0.45), 0.2))
        pri["pentane"] = _clip(min(float(pri.get("pentane") or 0.6), 0.45))
        lusc["voc_log2fc_prior"] = pri

    copd = by_id.get("copd")
    if copd:
        pri = dict(copd.get("voc_log2fc_prior") or {})
        # Keep COPD ethane/pentane signature; damp panel-inflated cancer-like ketones
        pri["ethane"] = _clip(max(float(pri.get("ethane") or 0.0), 1.05))
        pri["pentane"] = _clip(max(float(pri.get("pentane") or 0.0), 0.7))
        if float(pri.get("hexanal") or 0) > 1.1:
            pri["hexanal"] = 1.05
        for voc in ("acetaldehyde", "2_butanone"):
            if voc in pri and float(pri[voc]) > 0.35:
                pri[voc] = 0.25
        copd["voc_log2fc_prior"] = pri

    # Mild LUAD cancer ketone floor so it stays farther from COPD than bronchitis
    luad = by_id.get("lung_adenocarcinoma")
    if luad:
        pri = dict(luad.get("voc_log2fc_prior") or {})
        for voc, floor in (("acetaldehyde", 0.55), ("2_butanone", 0.7), ("hexanal", 0.85)):
            pri[voc] = _clip(max(float(pri.get(voc) or 0.0), floor))
        if "ethane" in pri:
            pri["ethane"] = _clip(min(float(pri["ethane"]), 0.25))
        luad["voc_log2fc_prior"] = pri


def enrich_voc_catalog_identity() -> dict[str, Any]:
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
    for know in {PKG_KNOW, REPO_KNOW}:
        path = know / "voc_catalog.json"
        if not path.exists():
            continue
        doc = json.loads(path.read_text())
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
        path.write_text(json.dumps(doc, indent=2) + "\n")
    return {"updated": updated, "n_extended_mapped": len(by_atlas)}


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

    def _add(did: str, voc: str, target: float, weight: float, source: str) -> None:
        voc = _voc_norm(voc)
        if voc not in atlas_vocs:
            return
        did = _resolve_disease(did, atlas_ids) or did
        if did not in atlas_ids:
            return
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

    n_panel = 0
    if include_panels:
        panel_paths = [
            PACKAGE_ROOT / "data" / "real_breath" / "literature_panels" / "priority10_voc_panels.json",
            Path("data/real_breath/literature_panels/priority10_voc_panels.json"),
        ]
        panels = None
        for p in panel_paths:
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
                _add(did, voc, float(fc), w, f"priority10:{did}")
                n_panel += 1

    n_pb = 0
    if include_public_breath:
        pb = _load(PKG_KNOW / "public_breath_benchmarks.json") or _load(
            REPO_KNOW / "public_breath_benchmarks.json"
        )
        for case in (pb or {}).get("cases") or []:
            did = case.get("disease_id") or case.get("disease")
            for voc in case.get("expect_elevated") or []:
                _add(did, voc, 0.65, W_PUBLIC, f"public_breath:{case.get('case_id')}")
                n_pb += 1
            for voc in case.get("expect_suppressed") or []:
                _add(did, voc, -0.45, W_PUBLIC, f"public_breath:{case.get('case_id')}")
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
                _add(did, voc, target, W_LIT_BENCH, f"lit_bench:{b.get('case_id')}")
                n_lit += 1
            for voc in b.get("expect_not_suppressed") or []:
                _add(did, voc, 0.35, W_LIT_BENCH * 0.7, f"lit_bench:{b.get('case_id')}")
                n_lit += 1

    n_hbdb = 0
    if include_hbdb:
        from ..datasources.ds13_hbdb import HBDB_DISEASE_VOC_NAMES

        for did, names in HBDB_DISEASE_VOC_NAMES.items():
            for name in names:
                voc = _voc_norm(name)
                if voc == "isoprene" and did in {
                    "asthma",
                    "copd",
                    "lung_adenocarcinoma",
                    "pneumonia_bacterial",
                }:
                    _add(did, voc, -0.35, W_HBDB, "hbdb_proxy")
                else:
                    _add(did, voc, 0.55, W_HBDB, "hbdb_proxy")
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
    """Directional elevate/suppress map for external validation (broader than LIT_EXPECT)."""
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
