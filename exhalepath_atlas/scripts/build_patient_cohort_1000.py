#!/usr/bin/env python3
"""Build a stratified 1000-patient diversity cohort for ExhalePath stress/eval."""

from __future__ import annotations

import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTS = [
    ROOT / "data" / "knowledge" / "patient_cohort_1000.json",
    ROOT / "src" / "exhalepath" / "data" / "knowledge" / "patient_cohort_1000.json",
]

# Oversampled focus triad + related pulmonary / oncology / systemic diversity.
FOCUS_TRIAD = ("copd", "chronic_bronchitis", "lung_adenocarcinoma")

STAGES_CANCER = ["I", "II", "III", "IV", None]
STAGES_OTHER = [None, "II", "III"]
SMOKING = ["never", "former", "current"]
SEXES = ["female", "male", "other"]

COMORBID_POOL = [
    "obesity",
    "type_2_diabetes",
    "hypertension",
    "heart_disease",
    "depression",
    "nafld",
    "asthma",
    "copd",
    "chronic_kidney_disease",
]


def _load_diseases() -> list[dict]:
    path = ROOT / "data" / "knowledge" / "disease_voc_priors.json"
    return list(json.loads(path.read_text())["diseases"])


def _load_tissues() -> list[str]:
    path = ROOT / "data" / "knowledge" / "whole_body_tissues.json"
    return [t["tissue_id"] for t in json.loads(path.read_text())["tissues"]]


def _load_gene_map(diseases: list[dict]) -> dict[str, list[str]]:
    """Merge disease driver_genes with offline Open Targets + GDC fragments."""
    out: dict[str, list[str]] = {}
    for d in diseases:
        genes = [str(g).upper() for g in (d.get("driver_genes") or [])]
        out[d["disease_id"]] = genes

    ot = ROOT / "data" / "knowledge" / "datasource_opentargets.json"
    if ot.exists():
        for row in json.loads(ot.read_text()).get("disease_genes") or []:
            did = row.get("disease_id")
            if not did:
                continue
            genes = [str(g.get("gene") or "").upper() for g in (row.get("genes") or [])]
            out.setdefault(did, [])
            out[did] = list(dict.fromkeys([*genes, *out[did]]))[:16]

    gdc = ROOT / "data" / "knowledge" / "datasource_gdc.json"
    if gdc.exists():
        for row in json.loads(gdc.read_text()).get("projects") or []:
            did = row.get("disease_id")
            if not did:
                continue
            genes = [str(g).upper() for g in (row.get("driver_genes") or [])]
            out.setdefault(did, [])
            out[did] = list(dict.fromkeys([*genes, *out[did]]))[:16]
    return out


def _pick_genes(rng: random.Random, pool: list[str], *, force: bool = False) -> list[str]:
    if not pool:
        return []
    if not force and rng.random() < 0.25:
        return []
    k = rng.randint(1, min(5, len(pool)))
    return rng.sample(pool, k)


def _pick_comorbidities(rng: random.Random, primary: str) -> list[str]:
    if rng.random() < 0.45:
        return []
    k = rng.randint(1, 3)
    pool = [c for c in COMORBID_POOL if c != primary]
    return rng.sample(pool, min(k, len(pool)))


def _patient(
    rng: random.Random,
    *,
    pid: str,
    disease_id: str,
    disease_name: str,
    category: str,
    default_site: str | None,
    tissues: list[str],
    gene_pool: list[str],
    stratum: str,
) -> dict:
    is_cancer = category == "cancer" or "cancer" in disease_id or disease_id.startswith("lung_")
    # Location: prefer default site, else random tissue; ~15% mismatched site (adversarial).
    if rng.random() < 0.15 and tissues:
        location = rng.choice(tissues)
        site_match = "mismatched"
    else:
        location = default_site or rng.choice(tissues)
        site_match = "matched"

    smoking = rng.choice(SMOKING)
    # Pulmonary + lung cancer: enrich current/former smokers
    if disease_id in {*FOCUS_TRIAD, "lung_squamous_cell_carcinoma", "asthma", "ild"}:
        smoking = rng.choices(SMOKING, weights=[0.2, 0.35, 0.45], k=1)[0]

    stage = None
    metastatic = False
    if is_cancer:
        stage = rng.choice(STAGES_CANCER)
        metastatic = bool(stage == "IV" or (stage in {"III", "IV"} and rng.random() < 0.4))
    elif rng.random() < 0.1:
        stage = rng.choice(STAGES_OTHER)

    age = int(rng.gauss(58 if is_cancer or disease_id == "copd" else 48, 14))
    age = max(18, min(95, age))
    sex = rng.choice(SEXES)
    genes = _pick_genes(rng, gene_pool, force=disease_id in FOCUS_TRIAD or is_cancer)
    comorbidities = _pick_comorbidities(rng, disease_id)
    weight = round(rng.uniform(0.35, 1.0), 2) if comorbidities else 0.65

    return {
        "patient_id": pid,
        "stratum": stratum,
        "site_match": site_match,
        "query": {
            "disease": disease_id,
            "location": location,
            "genes": genes,
            "comorbidities": comorbidities,
            "comorbidity_weight": weight,
            "age_years": age,
            "sex": sex,
            "smoking_status": smoking,
            "stage": stage,
            "metastatic": metastatic,
        },
        "meta": {
            "disease_id": disease_id,
            "disease_name": disease_name,
            "category": category,
            "n_genes": len(genes),
            "n_comorbidities": len(comorbidities),
        },
    }


def build(n: int = 1000, seed: int = 42) -> dict:
    rng = random.Random(seed)
    diseases = _load_diseases()
    by_id = {d["disease_id"]: d for d in diseases}
    tissues = _load_tissues()
    gene_map = _load_gene_map(diseases)

    # Stratified quotas (sum → n)
    quotas: list[tuple[str, str, int]] = [
        ("focus", "copd", 120),
        ("focus", "chronic_bronchitis", 100),
        ("focus", "lung_adenocarcinoma", 120),
        ("focus", "lung_squamous_cell_carcinoma", 60),
        ("pulmonary", "asthma", 40),
        ("pulmonary", "ild", 25),
        ("pulmonary", "cystic_fibrosis", 20),
        ("pulmonary", "ards", 20),
        ("pulmonary", "pneumonia_bacterial", 25),
        ("pulmonary", "tuberculosis", 20),
        ("pulmonary", "mesothelioma", 15),
        ("pulmonary", "sleep_apnea", 15),
    ]
    used = sum(q for _, _, q in quotas)
    # Fill remainder across all other atlas diseases
    other_ids = [
        d["disease_id"]
        for d in diseases
        if d["disease_id"] not in {did for _, did, _ in quotas}
    ]
    remain = n - used
    per = max(1, remain // max(1, len(other_ids)))
    for did in other_ids:
        quotas.append(("diverse", did, per))
    # Adjust last quota to hit exactly n
    total = sum(q for _, _, q in quotas)
    if total < n:
        quotas.append(("diverse", "type_2_diabetes", n - total))
    elif total > n:
        # trim from diverse tail
        overflow = total - n
        new_q = []
        for stratum, did, q in quotas:
            if overflow > 0 and stratum == "diverse" and q > 1:
                cut = min(overflow, q - 1)
                q -= cut
                overflow -= cut
            if q > 0:
                new_q.append((stratum, did, q))
        quotas = new_q

    patients: list[dict] = []
    i = 0
    for stratum, did, q in quotas:
        d = by_id.get(did)
        if d is None:
            continue
        for _ in range(q):
            if len(patients) >= n:
                break
            i += 1
            patients.append(
                _patient(
                    rng,
                    pid=f"P{i:04d}",
                    disease_id=did,
                    disease_name=d.get("name") or did,
                    category=d.get("category") or "unspecified",
                    default_site=d.get("default_site"),
                    tissues=tissues,
                    gene_pool=gene_map.get(did) or list(d.get("driver_genes") or []),
                    stratum=stratum,
                )
            )
        if len(patients) >= n:
            break

    # If short (missing disease ids), pad with random atlas diseases
    while len(patients) < n:
        d = rng.choice(diseases)
        i += 1
        patients.append(
            _patient(
                rng,
                pid=f"P{i:04d}",
                disease_id=d["disease_id"],
                disease_name=d.get("name") or d["disease_id"],
                category=d.get("category") or "unspecified",
                default_site=d.get("default_site"),
                tissues=tissues,
                gene_pool=gene_map.get(d["disease_id"]) or [],
                stratum="pad",
            )
        )

    patients = patients[:n]
    strata: dict[str, int] = {}
    disease_counts: dict[str, int] = {}
    for p in patients:
        strata[p["stratum"]] = strata.get(p["stratum"], 0) + 1
        did = p["meta"]["disease_id"]
        disease_counts[did] = disease_counts.get(did, 0) + 1

    return {
        "version": "1.0.0",
        "n_patients": len(patients),
        "seed": seed,
        "description": (
            "1000-patient diversity cohort spanning atlas diseases, anatomic sites, "
            "cancer stages, smoking backgrounds, comorbidities, and driver genes "
            "(Open Targets + GDC + curated). Oversamples COPD, chronic bronchitis, "
            "and lung cancer for VOC / genetic-shift embedding alignment."
        ),
        "focus_triad": list(FOCUS_TRIAD) + ["lung_squamous_cell_carcinoma"],
        "stratum_counts": strata,
        "disease_counts": dict(sorted(disease_counts.items(), key=lambda kv: -kv[1])),
        "resources_used": [
            "disease_voc_priors",
            "voc_catalog",
            "pathway_voc_map",
            "whole_body_tissues",
            "datasource_opentargets",
            "datasource_gdc",
            "cell_state_atlas",
            "voc_pathway_chains",
            "chembl_pathway_priors",
            "physio_constants",
        ],
        "patients": patients,
    }


def main() -> None:
    doc = build(1000, seed=42)
    text = json.dumps(doc, indent=2) + "\n"
    for p in OUTS:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        print(f"Wrote {p} ({doc['n_patients']} patients; focus={ {k: doc['disease_counts'].get(k) for k in doc['focus_triad']} })")


if __name__ == "__main__":
    main()
