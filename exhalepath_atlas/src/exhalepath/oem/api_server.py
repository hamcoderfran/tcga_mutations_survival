"""Minimal JSON API for procurement / OEM embed (stdlib — no FastAPI required).

  voc serve-api --port 8787

Endpoints:
  GET  /health
  GET  /v1/security
  GET  /v1/kit
  POST /v1/score-sample     {disease_id, vocs, signature_source?}
  POST /v1/verify-split     {manifest_path, study_id?}
  POST /v1/import-matrix    {path, format?, disease_id?}
  POST /v1/demo-close       {disease?, note?}
"""

from __future__ import annotations

import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlparse


SECURITY_STORY = {
    "positioning": "Cut the cost of wrong VOC panels — R&D diligence SaaS/API, not IVD",
    "not_a_medical_device": True,
    "controls": {
        "egress_allowlist": (
            "Open scientific fetches go through exhalepath.ingest.secure_fetch "
            "host allowlist (SSRF guard)."
        ),
        "no_arbitrary_code_in_uploads": (
            "Feature-table import accepts CSV/JSON only; executable extensions rejected "
            "by secure_fetch when downloading."
        ),
        "audit_logging": "Callers should log request_id + content hashes of imported matrices.",
        "data_residency": "Default deploy is customer VPC / air-gapped pip wheel + local API.",
        "secrets": "No patient PHI required; PatientTemplate is research vignette fields only.",
    },
    "soc2_ish_checklist": [
        "Access control: API behind customer IdP / mTLS (deploy-time)",
        "Encryption in transit: TLS terminator in front of voc serve-api",
        "Encryption at rest: customer volume encryption for runs/",
        "Change management: pin voc-breath version; verify wheel hashes",
        "Vendor risk: open-source core + allowlisted egress only",
        "Incident response: customer runbooks; no central PHI store in default mode",
    ],
    "gap_vs_full_soc2": (
        "This document is an architecture story for procurement — not a completed "
        "SOC 2 Type II attestation. Engage your auditor for formal certification."
    ),
}


def _json_response(handler: BaseHTTPRequestHandler, code: int, payload: Any) -> None:
    body = json.dumps(payload, default=str).encode("utf-8")
    handler.send_response(code)
    handler.send_header("Content-Type", "application/json")
    handler.send_header("Content-Length", str(len(body)))
    handler.send_header("X-ExhalePath-Positioning", "cut-cost-of-wrong-voc-panels")
    handler.end_headers()
    handler.wfile.write(body)


def _read_json(handler: BaseHTTPRequestHandler) -> dict[str, Any]:
    length = int(handler.headers.get("Content-Length") or 0)
    raw = handler.rfile.read(length) if length else b"{}"
    if not raw:
        return {}
    return json.loads(raw.decode("utf-8"))


class OemApiHandler(BaseHTTPRequestHandler):
    server_version = "ExhalePathOEM/1.0"

    def log_message(self, fmt: str, *args: Any) -> None:
        # quieter default
        return

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path.rstrip("/") or "/"
        if path == "/health":
            _json_response(self, 200, {"ok": True, "service": "exhalepath-oem"})
            return
        if path == "/v1/security":
            _json_response(self, 200, SECURITY_STORY)
            return
        if path == "/v1/kit":
            from .kit import kit_info

            _json_response(self, 200, kit_info())
            return
        _json_response(self, 404, {"error": "not_found", "path": path})

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path.rstrip("/")
        try:
            data = _read_json(self)
        except Exception as exc:  # noqa: BLE001
            _json_response(self, 400, {"error": "invalid_json", "detail": str(exc)})
            return
        try:
            if path == "/v1/score-sample":
                from .kit import score_sample

                out = score_sample(
                    data["disease_id"],
                    data["vocs"],
                    signature_source=data.get("signature_source") or "hybrid",
                    method=data.get("method") or "cosine",
                )
                _json_response(self, 200, out)
                return
            if path == "/v1/verify-split":
                from .kit import verify_locked_split

                out = verify_locked_split(
                    Path(data["manifest_path"]),
                    study_id=data.get("study_id"),
                )
                _json_response(self, 200, out)
                return
            if path == "/v1/import-matrix":
                from .kit import import_feature_table

                packed = import_feature_table(
                    Path(data["path"]),
                    fmt=data.get("format") or "auto",
                    disease_id=data.get("disease_id") or "malaria",
                    study_id=data.get("study_id") or "PARTNER_IMPORT",
                )
                # do not serialize full matrix in default response
                _json_response(
                    self,
                    200,
                    {k: v for k, v in packed.items() if k != "matrix"},
                )
                return
            if path == "/v1/demo-close":
                from ..eval.closing_demo import run_closing_demo

                out = run_closing_demo(
                    disease=data.get("disease") or "malaria",
                    note=data.get("note"),
                    out_dir=Path(data["out_dir"])
                    if data.get("out_dir")
                    else Path("runs/api_closing_demo"),
                    include_partner_loso=bool(data.get("include_partner_loso", True)),
                )
                _json_response(
                    self,
                    200,
                    {
                        "headline": out.get("headline"),
                        "disease": out.get("disease"),
                        "optimism_gap": out.get("optimism_gap"),
                        "citation_summary": out.get("citation_summary"),
                        "paper_pack": out.get("paper_pack"),
                        "artifacts": out.get("artifacts"),
                        "partner_loso": out.get("partner_loso"),
                    },
                )
                return
        except Exception as exc:  # noqa: BLE001
            _json_response(self, 500, {"error": "server_error", "detail": str(exc)})
            return
        _json_response(self, 404, {"error": "not_found", "path": path})


def serve_api(host: str = "127.0.0.1", port: int = 8787) -> None:
    httpd = ThreadingHTTPServer((host, port), OemApiHandler)
    print(f"ExhalePath OEM API on http://{host}:{port}  (GET /health, /v1/security)")
    httpd.serve_forever()


__all__ = ["SECURITY_STORY", "serve_api"]
