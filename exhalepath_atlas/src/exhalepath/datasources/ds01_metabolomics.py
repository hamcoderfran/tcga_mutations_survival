"""Priority 1 — Metabolomics Workbench / MetaboLights quantified breath studies."""

from __future__ import annotations

import json
import re
import statistics
from pathlib import Path
from typing import Any

import requests

from .base import DataSource

UA = {"User-Agent": "ExhalePathAtlas/1.0"}

# Curated public breath / VOC study accessions (always retained)
CURATED_STUDIES = [
    {
        "repo": "metabolomics_workbench",
        "accession": "ST003200",
        "title": "9 common breath VOCs in 504 healthy subjects (ppbv)",
        "disease_hints": ["healthy_baseline", "pulmonary"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST003200",
    },
    {
        "repo": "metabolomics_workbench",
        "accession": "ST002449",
        "title": "Valley fever / pneumonia breath VOCs",
        "disease_hints": ["pneumonia_bacterial", "infectious"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST002449",
    },
    {
        "repo": "metabolomics_workbench",
        "accession": "ST000883",
        "title": "Breathprinting malaria-associated biomarkers",
        "disease_hints": ["infectious", "malaria"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST000883",
    },
    {
        "repo": "metabolomics_workbench",
        "accession": "ST000587",
        "title": "Exhaled breath condensate in decompensated heart failure",
        "disease_hints": ["heart_failure", "cardiovascular"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST000587",
    },
    {
        "repo": "metabolomics_workbench",
        "accession": "ST001164",
        "title": "Cystic fibrosis acute pulmonary exacerbations (breath)",
        "disease_hints": ["cystic_fibrosis", "pulmonary"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST001164",
    },
    {
        "repo": "metabolights",
        "accession": "MTBLS2400",
        "title": "Epilepsy breath metabolomics (SESI-HRMS)",
        "disease_hints": ["epilepsy"],
        "url": "https://www.ebi.ac.uk/metabolights/MTBLS2400",
    },
    {
        "repo": "metabolights",
        "accession": "MTBLS2436",
        "title": "Breath / VOC MetaboLights study (GC-MS)",
        "disease_hints": ["pulmonary"],
        "url": "https://www.ebi.ac.uk/metabolights/MTBLS2436",
    },
    {
        "repo": "massive",
        "accession": "MSV000095340",
        "title": "Breath metabolites healthy vs asthmatic children",
        "disease_hints": ["asthma", "copd"],
        "url": "https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000095340",
    },
    {
        "repo": "figshare_scientific_data",
        "accession": "23522490",
        "title": "Clinical breathomics asthma/COPD/bronchiectasis",
        "disease_hints": ["asthma", "copd"],
        "url": "https://doi.org/10.6084/m9.figshare.23522490.v6",
    },
]

# Name aliases → ExhalePath voc_id
_NAME_TO_VOC = {
    "acetone": "acetone",
    "isoprene": "isoprene",
    "pentane": "pentane",
    "n-pentane": "pentane",
    "ethane": "ethane",
    "hexanal": "hexanal",
    "heptanal": "heptanal",
    "nonanal": "nonanal",
    "decanal": "decanal",
    "ethanol": "ethanol",
    "2-butanone": "2_butanone",
    "butanone": "2_butanone",
    "methyl ethyl ketone": "2_butanone",
    "dimethyl sulfide": "dms",
    "dimethylsulphide": "dms",
    "ammonia": "ammonia",
    "methanol": "methanol",
    "toluene": "toluene",
    "benzene": "benzene",
    "limonene": "limonene",
    "hydrogen sulfide": "hydrogen_sulfide",
    "acetaldehyde": "acetaldehyde",
    "formaldehyde": "formaldehyde",
    "indole": "indole",
    "phenol": "phenol",
    "dimethyl disulfide": "dimethyl_disulfide",
    "1-propanol": "propanol",
    "propanol": "propanol",
    "2-propanol": "isopropanol",
    "isopropanol": "isopropanol",
    "2-pentanone": "2_pentanone",
    "benzaldehyde": "benzaldehyde",
    "trimethylamine": "trimethylamine",
    "propanal": "propionaldehyde",
    "propionaldehyde": "propionaldehyde",
}


def _norm(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _map_voc(name: str) -> str | None:
    n = _norm(name)
    if not n:
        return None
    if n in _NAME_TO_VOC:
        return _NAME_TO_VOC[n]
    # Strip common adduct / ion suffixes then retry exact alias match only
    for suf in (" h2", " h", " m", " 1", " 2"):
        if n.endswith(suf):
            base = n[: -len(suf)].strip()
            if base in _NAME_TO_VOC:
                # adduct channels (e.g. Acetone_H2) are not the parent VOC ppb
                return None
    return None


def _summarize_values(vals: list[float]) -> dict[str, float]:
    vals = sorted(vals)
    n = len(vals)
    if n == 0:
        return {}
    def pct(p: float) -> float:
        if n == 1:
            return vals[0]
        i = min(n - 1, max(0, int(round(p * (n - 1)))))
        return vals[i]
    return {
        "n": float(n),
        "median": float(statistics.median(vals)),
        "mean": float(statistics.fmean(vals)),
        "q25": float(pct(0.25)),
        "q75": float(pct(0.75)),
        "min": float(vals[0]),
        "max": float(vals[-1]),
    }


class MetabolomicsReposSource(DataSource):
    priority = 1
    key = "metabolomics"
    title = "Metabolomics Workbench / MetaboLights / MassIVE"
    description = "Quantified public breath/VOC studies for validation and healthy baselines"

    def _mw_get_json(self, url: str) -> Any | None:
        r = requests.get(url, timeout=90, headers=UA)
        if not r.ok or not r.content or r.content.strip() in (b"[]", b""):
            return None
        # Prefer JSON object; /json suffix sometimes returns TSV
        try:
            return r.json()
        except Exception:  # noqa: BLE001
            return None

    def _harvest_mw_study(self, accession: str) -> dict[str, Any]:
        base = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{accession}"
        meta = self._mw_get_json(f"{base}/metabolites") or {}
        data = self._mw_get_json(f"{base}/data") or {}
        metabolites: list[dict[str, Any]] = []
        if isinstance(meta, dict):
            for row in meta.values():
                if not isinstance(row, dict):
                    continue
                name = row.get("metabolite_name") or row.get("refmet_name") or ""
                metabolites.append(
                    {
                        "metabolite_name": name,
                        "refmet_name": row.get("refmet_name"),
                        "metabolite_id": row.get("metabolite_id"),
                        "analysis_id": row.get("analysis_id"),
                        "voc_id": _map_voc(name),
                    }
                )
        quantified: list[dict[str, Any]] = []
        if isinstance(data, dict):
            for row in data.values():
                if not isinstance(row, dict):
                    continue
                name = row.get("metabolite_name") or row.get("refmet_name") or ""
                raw = row.get("DATA") or {}
                vals: list[float] = []
                if isinstance(raw, dict):
                    for v in raw.values():
                        try:
                            vals.append(float(v))
                        except (TypeError, ValueError):
                            continue
                summary = _summarize_values(vals)
                if not summary:
                    continue
                quantified.append(
                    {
                        "metabolite_name": name,
                        "refmet_name": row.get("refmet_name"),
                        "units": row.get("units"),
                        "analysis_id": row.get("analysis_id"),
                        "voc_id": _map_voc(name),
                        **summary,
                    }
                )
        return {
            "accession": accession,
            "n_metabolites": len(metabolites),
            "n_quantified": len(quantified),
            "metabolites": metabolites,
            "quantified": quantified,
            "mapped_voc_ids": sorted(
                {q["voc_id"] for q in quantified if q.get("voc_id")}
            ),
        }

    def _harvest_metabolights(self, accession: str) -> dict[str, Any]:
        out: dict[str, Any] = {"accession": accession, "files": [], "isa": {}}
        try:
            r = requests.get(
                f"https://www.ebi.ac.uk/metabolights/ws/studies/{accession}",
                timeout=60,
                headers=UA,
            )
            out["study_http_status"] = r.status_code
            if r.ok:
                try:
                    out["study_meta"] = r.json()
                except Exception:  # noqa: BLE001
                    out["study_meta_bytes"] = len(r.content)
        except Exception as e:  # noqa: BLE001
            out["study_error"] = str(e)
        try:
            r = requests.get(
                f"https://www.ebi.ac.uk/metabolights/ws/studies/{accession}/files",
                timeout=60,
                headers=UA,
            )
            if r.ok:
                payload = r.json()
                files = [f.get("file") for f in (payload.get("study") or []) if f.get("file")]
                out["files"] = files
                # Pull investigation + sample ISA headers (public FTP)
                for fname in files:
                    if not (
                        fname.startswith("i_")
                        or fname.startswith("s_")
                        or fname.startswith("a_")
                        or "maf" in fname.lower()
                    ):
                        continue
                    ftp = (
                        "https://ftp.ebi.ac.uk/pub/databases/metabolights/"
                        f"studies/public/{accession}/{fname}"
                    )
                    try:
                        fr = requests.get(ftp, timeout=90, headers=UA)
                        if not fr.ok:
                            continue
                        text = fr.text
                        # Keep a compact preview + metabolite-like tokens
                        lines = text.splitlines()
                        out["isa"][fname] = {
                            "n_lines": len(lines),
                            "preview": lines[:8],
                        }
                        # Extract CHEBI / metabolite names from maf if present
                        if "maf" in fname.lower() or fname.startswith("m_"):
                            names = []
                            for line in lines[1:400]:
                                cols = line.split("\t")
                                if cols:
                                    names.append(cols[0])
                            mapped = []
                            for n in names:
                                vid = _map_voc(n)
                                if vid:
                                    mapped.append({"name": n, "voc_id": vid})
                            out["isa"][fname]["mapped_vocs"] = mapped[:50]
                    except Exception as e:  # noqa: BLE001
                        out.setdefault("ftp_errors", []).append({fname: str(e)})
        except Exception as e:  # noqa: BLE001
            out["files_error"] = str(e)
        return out

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        studies = list(CURATED_STUDIES)
        live_mw: list[dict[str, Any]] = []
        quantified_by_study: dict[str, Any] = {}
        mtbls: list[dict[str, Any]] = []

        if not offline:
            # Discover additional MW breath/exhaled/VOC studies
            seen = {s["accession"] for s in studies}
            for term in ("breath", "exhaled", "VOC"):
                try:
                    r = requests.get(
                        f"https://www.metabolomicsworkbench.org/rest/study/study_title/{term}/summary",
                        timeout=60,
                        headers=UA,
                    )
                    if r.ok:
                        payload = r.json()
                        if isinstance(payload, dict):
                            for row in payload.values():
                                if not isinstance(row, dict):
                                    continue
                                sid = row.get("study_id")
                                if not sid or sid in seen:
                                    continue
                                title = row.get("study_title") or ""
                                # Prefer human breath-ish titles
                                tlow = title.lower()
                                if not any(
                                    k in tlow
                                    for k in (
                                        "breath",
                                        "exhaled",
                                        "volatile",
                                        "ebc",
                                    )
                                ):
                                    continue
                                # Skip obvious plant-only root VOC studies
                                if "arabidopsis" in tlow or "plant root" in tlow:
                                    continue
                                studies.append(
                                    {
                                        "repo": "metabolomics_workbench",
                                        "accession": sid,
                                        "title": title,
                                        "disease_hints": ["pulmonary"],
                                        "url": (
                                            "https://www.metabolomicsworkbench.org/"
                                            f"data/DRCCMetadata.php?StudyID={sid}"
                                        ),
                                        "discovered_via": term,
                                    }
                                )
                                seen.add(sid)
                            live_mw.append({"query": term, "n_hits": len(payload)})
                except Exception as e:  # noqa: BLE001
                    live_mw.append({"query": term, "error": str(e)})

            # Quantified pulls for MW accessions
            for s in studies:
                if s["repo"] != "metabolomics_workbench":
                    continue
                acc = s["accession"]
                try:
                    quantified_by_study[acc] = self._harvest_mw_study(acc)
                except Exception as e:  # noqa: BLE001
                    quantified_by_study[acc] = {"accession": acc, "error": str(e)}

            # MetaboLights ISA / metadata
            for s in studies:
                if s["repo"] != "metabolights":
                    continue
                try:
                    mtbls.append(self._harvest_metabolights(s["accession"]))
                except Exception as e:  # noqa: BLE001
                    mtbls.append({"accession": s["accession"], "error": str(e)})
        else:
            # Reuse prior online harvest if present
            prev_cat = self.out_dir / "breath_study_catalog.json"
            if prev_cat.exists():
                prev_doc = json.loads(prev_cat.read_text())
                if prev_doc.get("studies"):
                    studies = prev_doc["studies"]
                live_mw = prev_doc.get("live_mw_discovery") or []
            prev = self.out_dir / "mw_quantified_studies.json"
            if prev.exists():
                quantified_by_study = json.loads(prev.read_text())
            prev_m = self.out_dir / "metabolights_studies.json"
            if prev_m.exists():
                mtbls = json.loads(prev_m.read_text()).get("studies") or []
            prev_a = self.out_dir / "atlas_mapped_quantified.json"
            # Keep atlas_mapped rows from prior harvest (avoid empty rewrite)
            if prev_a.exists() and not quantified_by_study:
                pass

        catalog = {
            "version": "2.0.0",
            "n_curated_studies": len(studies),
            "studies": studies,
            "live_mw_discovery": live_mw,
            "n_mw_quantified_studies": sum(
                1
                for v in quantified_by_study.values()
                if isinstance(v, dict) and v.get("n_quantified", 0) > 0
            ),
        }
        path = self.write_json("breath_study_catalog.json", catalog)
        q_path = self.write_json("mw_quantified_studies.json", quantified_by_study)
        m_path = self.write_json(
            "metabolights_studies.json", {"version": "1.0.0", "studies": mtbls}
        )

        # Compact atlas-mapped healthy / study summaries
        atlas_rows = []
        for acc, doc in quantified_by_study.items():
            if not isinstance(doc, dict):
                continue
            for q in doc.get("quantified") or []:
                if not q.get("voc_id"):
                    continue
                atlas_rows.append(
                    {
                        "study_id": acc,
                        "voc_id": q["voc_id"],
                        "metabolite_name": q.get("metabolite_name"),
                        "units": q.get("units"),
                        "median": q.get("median"),
                        "mean": q.get("mean"),
                        "q25": q.get("q25"),
                        "q75": q.get("q75"),
                        "n": q.get("n"),
                    }
                )
        a_path = self.write_json(
            "atlas_mapped_quantified.json",
            {"version": "1.0.0", "n_rows": len(atlas_rows), "rows": atlas_rows},
        )

        idx: dict[str, list[str]] = {}
        for s in studies:
            for d in s["disease_hints"]:
                idx.setdefault(d, []).append(s["accession"])
        idx_path = self.write_json("disease_study_index.json", idx)
        man = self.write_manifest(
            n_studies=len(studies),
            n_quantified_rows=len(atlas_rows),
            artifacts=[str(path), str(q_path), str(m_path), str(a_path), str(idx_path)],
        )
        return {
            "catalog": path,
            "quantified": q_path,
            "metabolights": m_path,
            "atlas_mapped": a_path,
            "disease_index": idx_path,
            "manifest": man,
        }

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        cat = json.loads((self.out_dir / "breath_study_catalog.json").read_text())
        atlas = json.loads((self.out_dir / "atlas_mapped_quantified.json").read_text())
        mtbls = json.loads((self.out_dir / "metabolights_studies.json").read_text())
        fused = {
            "version": "2.0.0",
            "n_studies": cat["n_curated_studies"],
            "n_mw_quantified_studies": cat.get("n_mw_quantified_studies"),
            "n_atlas_mapped_rows": atlas.get("n_rows"),
            "studies": cat["studies"],
            "atlas_mapped_quantified": atlas.get("rows") or [],
            "metabolights": mtbls.get("studies") or [],
            "live_mw_discovery": cat.get("live_mw_discovery"),
        }
        out = knowledge_dir / "datasource_metabolomics.json"
        out.write_text(json.dumps(fused, indent=2))

        # Refine healthy baselines from ST003200 (ppbv) when available
        healthy: dict[str, dict[str, float]] = {}
        for row in atlas.get("rows") or []:
            if row.get("study_id") != "ST003200":
                continue
            units = (row.get("units") or "").lower()
            if not (
                "ppb" in units
                or "parts per billion" in units
                or units in {"ppbv", "ppb"}
            ):
                continue
            vid = row["voc_id"]
            healthy[vid] = {
                "median_ppb": float(row["median"]),
                "q25_ppb": float(row["q25"]),
                "q75_ppb": float(row["q75"]),
                "n": float(row["n"]),
                "source_study": "ST003200",
            }
        if healthy:
            hb_path = knowledge_dir / "mw_healthy_breath_baselines.json"
            hb_path.write_text(
                json.dumps(
                    {
                        "version": "1.0.0",
                        "source": "Metabolomics Workbench ST003200",
                        "unit": "ppbv",
                        "vocs": healthy,
                    },
                    indent=2,
                )
            )
            voc_path = knowledge_dir / "voc_catalog.json"
            if voc_path.exists():
                voc_doc = json.loads(voc_path.read_text())
                for v in voc_doc.get("vocs") or []:
                    h = healthy.get(v["voc_id"])
                    if not h:
                        continue
                    v["mw_healthy_ppb_median"] = h["median_ppb"]
                    v["mw_healthy_ppb_q25"] = h["q25_ppb"]
                    v["mw_healthy_ppb_q75"] = h["q75_ppb"]
                    v["mw_healthy_n"] = h["n"]
                voc_path.write_text(json.dumps(voc_doc, indent=2))

        return {
            "path": str(out),
            "n_studies": fused["n_studies"],
            "n_atlas_mapped_rows": fused["n_atlas_mapped_rows"],
            "n_healthy_baselines": len(healthy),
        }
