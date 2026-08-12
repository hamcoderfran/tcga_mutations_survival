"""Tests for secure fetch + coverage audit."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from exhalepath.ingest.secure_fetch import SecureFetchError, host_allowed, secure_fetch, validate_json_artifact
from exhalepath.ingest.http import request_json
from exhalepath.eval.coverage_audit import (
    build_coverage_universe,
    expand_extended_catalog,
    run_coverage_audit,
)


def test_secure_fetch_blocks_non_allowlisted_host(tmp_path: Path):
    with pytest.raises(SecureFetchError, match="allowlist"):
        secure_fetch("https://evil.example/malware.exe", dest=tmp_path / "x.bin")


def test_secure_fetch_blocks_script_extension(tmp_path: Path):
    with pytest.raises(SecureFetchError):
        secure_fetch(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search.py",
            dest=tmp_path / "x.py",
        )


def test_host_allowed_core_apis():
    assert host_allowed("https://www.metabolomicsworkbench.org/rest/study/ST000001/summary")
    assert host_allowed("https://www.ebi.ac.uk/europepmc/webservices/rest/search")
    assert host_allowed("https://api.gdc.cancer.gov/cases")
    assert host_allowed("https://reactome.org/ContentService/data/query/abc")
    assert not host_allowed("https://evil.example/x")
    assert not host_allowed("file:///etc/passwd")
    # S3 only as redirect target, not as initial fetch origin
    assert not host_allowed("https://s3-eu-west-1.amazonaws.com/bucket/file.csv")
    assert host_allowed(
        "https://s3-eu-west-1.amazonaws.com/bucket/file.csv",
        allow_storage_cdn=True,
    )


def test_request_json_blocks_non_allowlisted():
    with pytest.raises(Exception, match="allowlist"):
        request_json("GET", "https://evil.example/api.json")


def test_validate_json_artifact_ok(tmp_path: Path):
    p = tmp_path / "a.json"
    p.write_text(json.dumps({"compounds": [], "n_compounds": 0}))
    doc = validate_json_artifact(p, required_keys=["compounds"])
    assert "compounds" in doc


def test_validate_json_rejects_executable(tmp_path: Path):
    p = tmp_path / "a.json"
    p.write_bytes(b"\x7fELF" + b"0" * 20)
    with pytest.raises(SecureFetchError):
        validate_json_artifact(p)


def test_validate_json_rejects_code_keys(tmp_path: Path):
    p = tmp_path / "a.json"
    p.write_text(json.dumps({"__import__": "os", "compounds": []}))
    with pytest.raises(SecureFetchError, match="Suspicious"):
        validate_json_artifact(p, required_keys=["compounds"])


def test_coverage_universe_defines_99pct_honestly():
    u = build_coverage_universe()
    assert "epa_volatilome" in {x["id"] for x in u["in_scope"]}
    assert any(x["id"] == "owlstone_voc_atlas" for x in u["out_of_scope_blocked"])
    assert "not 99% of all online" in u["ninety_nine_percent_target"]["note"].lower() or \
        "not" in u["ninety_nine_percent_target"]["note"].lower()


def test_expand_extended_catalog_offline():
    out = expand_extended_catalog(offline=True)
    assert out["n_compounds"] >= 700
    assert out["security"]["downloads_executed"] is False


def test_coverage_audit_offline_meets_open_compound_target():
    report = run_coverage_audit(expand=True, offline=True)
    assert report["metrics"]["open_compound_coverage"] >= 0.99
    assert report["metrics"]["meets_99pct_open_compound_target"] is True
    assert Path("data/knowledge/COVERAGE_AUDIT.md").exists()
    assert Path("data/datasources/INTEGRITY_MANIFEST.json").exists()
    assert report["security"]["no_execution_of_downloads"] is True
    gap_ids = {g["id"] for g in report.get("gaps") or []}
    assert "magdeburg_no_intensity_matrix" in gap_ids
    assert "mh_thin_conditions_unpromoted" in gap_ids
    assert report["metrics"]["literature_panels"] >= 3
