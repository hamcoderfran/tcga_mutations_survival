"""Priority 1 — Metabolomics Workbench / MetaboLights / MassIVE breath studies."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from .base import DataSource

# Curated public breath / VOC study accessions
CURATED_STUDIES = [
    {
        "repo": "metabolomics_workbench",
        "accession": "ST002449",
        "title": "Valley fever / pneumonia breath VOCs",
        "disease_hints": ["pneumonia_bacterial", "infectious"],
        "url": "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST002449",
    },
    {
        "repo": "metabolights",
        "accession": "MTBLS2400",
        "title": "Epilepsy breath metabolomics (SESI-HRMS)",
        "disease_hints": ["epilepsy"],
        "url": "https://www.ebi.ac.uk/metabolights/MTBLS2400",
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
    {
        "repo": "metabolomics_workbench",
        "accession": "ST001236",
        "title": "Breath VOC related MW study (catalogued)",
        "disease_hints": ["pulmonary"],
        "url": "https://www.metabolomicsworkbench.org/",
    },
]


class MetabolomicsReposSource(DataSource):
    priority = 1
    key = "metabolomics"
    title = "Metabolomics Workbench / MetaboLights / MassIVE"
    description = "Public breath/VOC study catalog for external validation"

    def harvest(self, *, offline: bool = False) -> dict[str, Path]:
        studies = list(CURATED_STUDIES)
        live = []
        if not offline:
            # MW REST list search for breath
            try:
                r = requests.get(
                    "https://www.metabolomicsworkbench.org/rest/study/study_title/breath/fetch",
                    timeout=60,
                    headers={"User-Agent": "ExhalePathAtlas/1.0"},
                )
                if r.ok and r.text.strip():
                    # API may return TSV/JSON depending on endpoint
                    live.append({"query": "breath", "status": r.status_code, "bytes": len(r.content)})
            except Exception as e:  # noqa: BLE001
                live.append({"query": "breath", "error": str(e)})
            try:
                r = requests.get(
                    "https://www.ebi.ac.uk/metabolights/ws/studies",
                    params={"query": "breath VOC", "pageSize": 20},
                    timeout=60,
                    headers={"User-Agent": "ExhalePathAtlas/1.0"},
                )
                live.append({"metabolights": r.status_code, "ok": r.ok})
            except Exception as e:  # noqa: BLE001
                live.append({"metabolights_error": str(e)})

        catalog = {
            "version": "1.0.0",
            "n_curated_studies": len(studies),
            "studies": studies,
            "live_probe": live,
        }
        path = self.write_json("breath_study_catalog.json", catalog)
        # disease → study index
        idx: dict[str, list[str]] = {}
        for s in studies:
            for d in s["disease_hints"]:
                idx.setdefault(d, []).append(s["accession"])
        idx_path = self.write_json("disease_study_index.json", idx)
        man = self.write_manifest(n_studies=len(studies), artifacts=[str(path), str(idx_path)])
        return {"catalog": path, "disease_index": idx_path, "manifest": man}

    def fuse(self, knowledge_dir: Path) -> dict[str, Any]:
        cat = json.loads((self.out_dir / "breath_study_catalog.json").read_text())
        out = knowledge_dir / "datasource_metabolomics.json"
        out.write_text(json.dumps(cat, indent=2))
        return {"path": str(out), "n_studies": cat["n_curated_studies"]}
