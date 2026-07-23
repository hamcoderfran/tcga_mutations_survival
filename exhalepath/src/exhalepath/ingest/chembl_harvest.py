"""
Harvest ChEMBL bioactivities for ExhalePath pathway seed genes + VOC molecules.

Why ChEMBL (∼24M activities) helps ExhalePath
---------------------------------------------
ChEMBL does **not** provide exhaled VOC ppb labels. It provides the largest
public chemogenomic graph for the enzymes/genes in VOC biosynthetic pathways.

Training uses:
  1) activity rows (molecule × target × pChEMBL) for pathway genes
  2) distilled pathway ligandability priors → VOC emission features
  3) VOC physicochemical properties (AlogP, MW, PSA) → blood–air / ADME priors

Scale knobs:
  --max-per-target / --max-rows cap API harvests
  --offline-demo synthesizes multi-million-scale rows for local training tests
  Full dump: point train-chembl at a ChEMBL SQLite extract of the same schema
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ..config import CHEMBL_DIR, KNOWLEDGE_DIR
from ..knowledge.loader import KnowledgeBase
from .chembl_client import ChEMBLClient


def _gene_to_pathways(kb: KnowledgeBase) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for pid, p in kb.pathways.items():
        for g in p.get("seed_genes") or []:
            out.setdefault(str(g).upper(), []).append(pid)
    return out


def _activity_row(a: dict[str, Any], *, gene: str, pathways: list[str], target_id: str) -> dict[str, Any]:
    smiles = a.get("canonical_smiles") or ""
    return {
        "activity_id": a.get("activity_id"),
        "molecule_chembl_id": a.get("molecule_chembl_id"),
        "target_chembl_id": target_id,
        "gene": gene,
        "pathways": "|".join(pathways),
        "pchembl_value": float(a["pchembl_value"]) if a.get("pchembl_value") not in (None, "") else np.nan,
        "standard_type": a.get("standard_type"),
        "standard_value": a.get("standard_value"),
        "standard_units": a.get("standard_units"),
        "assay_type": a.get("assay_type"),
        "canonical_smiles": smiles,
        "smiles_len": len(smiles),
        "document_chembl_id": a.get("document_chembl_id"),
    }


def harvest_chembl_for_pathways(
    *,
    out_dir: Path | None = None,
    max_per_target: int = 5_000,
    max_rows: int | None = 500_000,
    max_genes: int | None = None,
    offline_demo: bool = False,
    demo_rows: int = 1_000_000,
    include_voc_properties: bool = True,
    install_knowledge: bool = True,
) -> dict[str, Path]:
    """
    Build ChEMBL training tables + distilled priors for ExhalePath.

    Returns paths to activity CSV, gene/target map, VOC properties, priors JSON.
    """
    out_dir = Path(out_dir or CHEMBL_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    kb = KnowledgeBase()
    gene_pw = _gene_to_pathways(kb)
    genes = sorted(gene_pw)
    if max_genes is not None:
        genes = genes[: max(1, max_genes)]

    if offline_demo:
        activities = _synthesize_demo_activities(genes, gene_pw, n_rows=demo_rows)
        target_map = [
            {
                "gene": g,
                "target_chembl_id": f"DEMO_{g}",
                "pref_name": f"Demo target {g}",
                "organism": "Homo sapiens",
                "n_activities": int((activities["gene"] == g).sum()),
            }
            for g in genes
        ]
        voc_props = _synthesize_voc_properties(kb)
    else:
        client = ChEMBLClient()
        target_map = []
        rows: list[dict[str, Any]] = []
        # Prefer well-annotated metabolic / oncogene targets first (denser ChEMBL coverage)
        priority = {
            "KRAS",
            "CYP2E1",
            "CYP1A1",
            "CYP3A4",
            "HMGCR",
            "IDO1",
            "ALOX15",
            "LDHA",
            "HK2",
            "GPX4",
            "PTGS2",
            "HIF1A",
            "MTOR",
            "PIK3CA",
        }
        genes = sorted(genes, key=lambda g: (0 if g in priority else 1, g))
        for gene in genes:
            if max_rows is not None and len(rows) >= max_rows:
                break
            targets = client.find_targets_for_gene(gene)
            if not targets:
                target_map.append(
                    {
                        "gene": gene,
                        "target_chembl_id": None,
                        "pref_name": None,
                        "organism": None,
                        "n_activities": 0,
                    }
                )
                continue
            n_here = 0
            chosen = None
            # Try up to 3 ranked targets until we collect activities
            for t in targets[:3]:
                if max_rows is not None and len(rows) >= max_rows:
                    break
                if n_here >= max_per_target:
                    break
                tid = t["target_chembl_id"]
                budget = max_per_target - n_here
                if max_rows is not None:
                    budget = min(budget, max_rows - len(rows))
                got = 0
                for a in client.iter_activities(target_chembl_id=tid, max_rows=budget):
                    rows.append(
                        _activity_row(
                            a, gene=gene, pathways=gene_pw.get(gene, []), target_id=tid
                        )
                    )
                    n_here += 1
                    got += 1
                if got and chosen is None:
                    chosen = t
            t = chosen or targets[0]
            target_map.append(
                {
                    "gene": gene,
                    "target_chembl_id": t.get("target_chembl_id"),
                    "pref_name": t.get("pref_name"),
                    "organism": t.get("organism"),
                    "n_activities": n_here,
                }
            )
        activities = pd.DataFrame(rows)
        voc_props = (
            harvest_voc_molecule_properties(kb, client=client)
            if include_voc_properties
            else pd.DataFrame()
        )

    act_path = out_dir / "chembl_pathway_activities.csv"
    map_path = out_dir / "chembl_gene_target_map.csv"
    voc_path = out_dir / "chembl_voc_properties.csv"
    prior_path = out_dir / "chembl_pathway_priors.json"
    manifest_path = out_dir / "chembl_manifest.json"

    activities.to_csv(act_path, index=False)
    pd.DataFrame(target_map).to_csv(map_path, index=False)
    if len(voc_props):
        voc_props.to_csv(voc_path, index=False)
    priors = distill_chembl_priors(activities, kb=kb, voc_props=voc_props)
    prior_path.write_text(json.dumps(priors, indent=2))
    know_prior = KNOWLEDGE_DIR / "chembl_pathway_priors.json"
    if install_knowledge:
        # Runtime predict features read this knowledge file
        know_prior.write_text(json.dumps(priors, indent=2))

    manifest = {
        "version": "1.0.0",
        "source": "offline_demo" if offline_demo else "chembl_api",
        "n_activity_rows": int(len(activities)),
        "n_genes_queried": len(genes),
        "n_genes_with_targets": int(sum(1 for t in target_map if t.get("target_chembl_id"))),
        "n_vocs_with_properties": int(len(voc_props)),
        "max_per_target": max_per_target,
        "max_rows": max_rows,
        "note": (
            "Rows supervise chemogenomic pathway modulation / VOC physchem — "
            "NOT exhaled breath ppb. Distill into priors, then blend with physiology."
        ),
        "paths": {
            "activities": str(act_path),
            "gene_target_map": str(map_path),
            "voc_properties": str(voc_path),
            "priors": str(prior_path),
            "knowledge_priors": str(know_prior) if install_knowledge else None,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return {
        "activities": act_path,
        "gene_target_map": map_path,
        "voc_properties": voc_path,
        "priors": prior_path,
        "manifest": manifest_path,
    }


def harvest_voc_molecule_properties(
    kb: KnowledgeBase | None = None,
    *,
    client: ChEMBLClient | None = None,
) -> pd.DataFrame:
    """Resolve VOC catalog names to ChEMBL molecules and pull physicochemical props."""
    kb = kb or KnowledgeBase()
    client = client or ChEMBLClient()
    rows = []
    for voc_id, voc in kb.vocs.items():
        queries = [voc.get("name"), voc_id.replace("_", " "), voc_id]
        mol = None
        hit_q = None
        for q in queries:
            if not q:
                continue
            # exact pref_name first
            try:
                data = client.search_molecule(str(q), limit=5)
            except Exception:  # noqa: BLE001
                data = []
            for m in data:
                pref = (m.get("pref_name") or "").lower()
                if pref == str(q).lower() or voc_id.replace("_", " ") in pref:
                    mol = m
                    hit_q = q
                    break
            if mol is None and data:
                mol = data[0]
                hit_q = q
            if mol:
                break
        if not mol:
            rows.append({"voc_id": voc_id, "molecule_chembl_id": None, "matched": False})
            continue
        props = mol.get("molecule_properties") or {}
        rows.append(
            {
                "voc_id": voc_id,
                "name": voc.get("name"),
                "molecule_chembl_id": mol.get("molecule_chembl_id"),
                "matched": True,
                "query": hit_q,
                "full_mwt": _f(props.get("full_mwt")),
                "alogp": _f(props.get("alogp")),
                "psa": _f(props.get("psa")),
                "hba": _f(props.get("hba")),
                "hbd": _f(props.get("hbd")),
                "rtb": _f(props.get("rtb")),
                "aromatic_rings": _f(props.get("aromatic_rings")),
                "cx_logp": _f(props.get("cx_logp")),
                "cx_logd": _f(props.get("cx_logd")),
            }
        )
    return pd.DataFrame(rows)


def distill_chembl_priors(
    activities: pd.DataFrame,
    *,
    kb: KnowledgeBase | None = None,
    voc_props: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """
    Collapse millions of activity rows into pathway/VOC priors usable at predict time.

    pathway_ligandability ≈ robust mean pChEMBL × log10(n_actives)
    (high = chemically tractable / heavily modulated enzyme set)
    """
    kb = kb or KnowledgeBase()
    pathway_stats: dict[str, Any] = {}
    if len(activities) and "pathways" in activities.columns:
        explode = activities.dropna(subset=["pchembl_value"]).copy()
        explode["pathway_id"] = explode["pathways"].astype(str).str.split("|")
        explode = explode.explode("pathway_id")
        explode = explode[explode["pathway_id"].astype(str).str.len() > 0]
        for pid, g in explode.groupby("pathway_id"):
            vals = g["pchembl_value"].astype(float)
            n = int(len(vals))
            mean_p = float(vals.mean()) if n else 0.0
            # ligandability score in ~[0, 1]
            lig = float(np.clip((mean_p - 4.0) / 4.0, 0, 1) * np.log10(n + 1) / 4.0)
            pathway_stats[str(pid)] = {
                "n_activities": n,
                "n_molecules": int(g["molecule_chembl_id"].nunique())
                if "molecule_chembl_id" in g
                else n,
                "mean_pchembl": mean_p,
                "median_pchembl": float(vals.median()) if n else 0.0,
                "ligandability": lig,
            }

    voc_phys = {}
    if voc_props is not None and len(voc_props):
        for r in voc_props.to_dict(orient="records"):
            if not r.get("matched"):
                continue
            alogp = r.get("alogp")
            # Heuristic: higher AlogP → lower blood:air λ tendency for hydrocarbons
            lambda_hint = None
            if alogp is not None and not (isinstance(alogp, float) and np.isnan(alogp)):
                # map alogp [-1,5] → lambda hint multiplier around physio defaults
                lambda_hint = float(np.clip(10 ** (1.5 - 0.35 * float(alogp)), 0.2, 2000))
            voc_phys[r["voc_id"]] = {
                "molecule_chembl_id": r.get("molecule_chembl_id"),
                "full_mwt": r.get("full_mwt"),
                "alogp": r.get("alogp"),
                "psa": r.get("psa"),
                "lambda_blood_air_hint": lambda_hint,
            }

    return {
        "version": "1.0.0",
        "description": (
            "Distilled ChEMBL chemogenomic priors for ExhalePath. "
            "ligandability scales pathway→VOC emission confidence; "
            "VOC AlogP hints refine blood–air partition."
        ),
        "n_source_activities": int(len(activities)),
        "pathways": pathway_stats,
        "voc_physchem": voc_phys,
    }


def _synthesize_demo_activities(
    genes: list[str], gene_pw: dict[str, list[str]], *, n_rows: int
) -> pd.DataFrame:
    """Multi-million-scale synthetic ChEMBL-like table for offline training drills.

    Includes a learnable signal: pChEMBL rises with smiles_len and pathway prior
    so `train-chembl` can show non-trivial R² on demo data.
    """
    rng = np.random.default_rng(7)
    n = int(max(1_000, n_rows))
    gene_arr = rng.choice(np.array(genes, dtype=object), size=n)
    smiles_len = rng.integers(10, 120, size=n)
    std_types = rng.choice(np.array(["IC50", "Ki", "Kd", "EC50"], dtype=object), size=n)
    # pathway-specific base potency
    pw_base = {
        "lipid_peroxidation": 0.4,
        "glycolysis_warburg": 0.2,
        "ketone_body_metabolism": 0.5,
        "cytochrome_p450_detox": 0.7,
        "Kras_mapk_proliferation": 0.6,
        "fatty_acid_oxidation": 0.35,
    }
    base = np.array(
        [
            max((pw_base.get(p, 0.1) for p in gene_pw.get(str(g), [])), default=0.1)
            for g in gene_arr
        ]
    )
    std_boost = np.array([{ "IC50": 0.0, "Ki": 0.15, "Kd": 0.1, "EC50": -0.1}[s] for s in std_types])
    pchembl = (
        5.2
        + base
        + 0.015 * smiles_len.astype(float)
        + std_boost
        + rng.normal(0, 0.35, size=n)
    ).clip(3.0, 10.0)
    rows = {
        "activity_id": np.arange(n),
        "molecule_chembl_id": [f"CHEMBL_DEMO_{i%50000}" for i in range(n)],
        "target_chembl_id": [f"DEMO_{g}" for g in gene_arr],
        "gene": gene_arr,
        "pathways": ["|".join(gene_pw.get(str(g), ["unknown"])) for g in gene_arr],
        "pchembl_value": pchembl,
        "standard_type": std_types,
        "standard_value": 10 ** (9 - pchembl),
        "standard_units": "nM",
        "assay_type": "B",
        "canonical_smiles": ["C" * int(l) for l in smiles_len],
        "smiles_len": smiles_len,
        "document_chembl_id": "DEMO_DOC",
    }
    return pd.DataFrame(rows)


def _synthesize_voc_properties(kb: KnowledgeBase) -> pd.DataFrame:
    rng = np.random.default_rng(3)
    rows = []
    for i, (voc_id, voc) in enumerate(kb.vocs.items()):
        # crude MW from formula if present
        formula = voc.get("formula") or ""
        mw = 50 + 10 * len(formula)
        alogp = float(rng.normal(1.2, 1.0))
        rows.append(
            {
                "voc_id": voc_id,
                "name": voc.get("name"),
                "molecule_chembl_id": f"CHEMBL_VOC_{i}",
                "matched": True,
                "query": voc.get("name"),
                "full_mwt": mw,
                "alogp": alogp,
                "psa": float(max(0, rng.normal(30, 20))),
                "hba": float(rng.integers(0, 4)),
                "hbd": float(rng.integers(0, 3)),
                "rtb": float(rng.integers(0, 6)),
                "aromatic_rings": float(rng.integers(0, 2)),
                "cx_logp": alogp,
                "cx_logd": alogp - 0.2,
            }
        )
    return pd.DataFrame(rows)


def _f(x: Any) -> float | None:
    if x is None or x == "":
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def apply_chembl_lambda_hints_to_physio(
    priors_path: Path | None = None,
    physio_path: Path | None = None,
) -> Path:
    """Merge ChEMBL AlogP-derived λ hints into physio_constants (setdefault only)."""
    priors_path = Path(priors_path or KNOWLEDGE_DIR / "chembl_pathway_priors.json")
    physio_path = Path(physio_path or KNOWLEDGE_DIR / "physio_constants.json")
    if not priors_path.exists():
        raise FileNotFoundError(priors_path)
    priors = json.loads(priors_path.read_text())
    phys = json.loads(physio_path.read_text())
    lam = phys.setdefault("blood_air_partition_lambda", {})
    updated = 0
    for voc_id, props in (priors.get("voc_physchem") or {}).items():
        hint = props.get("lambda_blood_air_hint")
        if hint is None:
            continue
        if voc_id not in lam:
            lam[voc_id] = float(hint)
            updated += 1
    phys["chembl_lambda_hints_applied"] = updated
    physio_path.write_text(json.dumps(phys, indent=2))
    return physio_path
