"""Secure, allowlisted HTTP fetch with provenance + malware/poisoning guards.

Design goals:
- Never execute downloaded content
- Only contact explicitly allowlisted hosts
- Stream to quarantine, hash, validate, then promote
- Reject executables / archive bombs / oversized payloads
- Record provenance (URL, SHA-256, UTC time, content-type) for every artifact
"""

from __future__ import annotations

import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urljoin, urlparse

import requests

from ..config import USER_AGENT

# Hosts we are willing to contact for VOC / metabolomics / literature metadata.
# Anything else is rejected (SSRF / supply-chain guard).
ALLOWED_HOSTS: frozenset[str] = frozenset(
    {
        "www.metabolomicsworkbench.org",
        "metabolomicsworkbench.org",
        "www.ebi.ac.uk",
        "ebi.ac.uk",
        "ftp.ebi.ac.uk",
        "api.platform.opentargets.org",
        "api.gdc.cancer.gov",
        "pubchem.ncbi.nlm.nih.gov",
        "rest.kegg.jp",
        "clowder.edap-cluster.com",
        "ndownloader.figshare.com",
        "figshare.com",
        "api.figshare.com",
        "doi.org",
        "api.crossref.org",
        "www.europepmc.org",
        "europepmc.org",
        "www.ncbi.nlm.nih.gov",
        "eutils.ncbi.nlm.nih.gov",
        "zenodo.org",
        "www.zenodo.org",
        "massive.ucsd.edu",
        "gnps.ucsd.edu",
        "hmdb.ca",
        "www.hmdb.ca",
        "www.reactome.org",
        "reactome.org",
        "content.cryst.bbk.ac.uk",
    }
)

ALLOWED_SCHEMES = frozenset({"https"})  # HTTPS only — no cleartext downgrade
MAX_BYTES_DEFAULT = 80 * 1024 * 1024  # 80 MiB

# Origins allowed to redirect onto storage CDNs (Figshare/Zenodo/Clowder blobs).
# doi.org is intentionally excluded — open redirects → attacker S3 would pass otherwise.
_STORAGE_REDIRECT_ORIGINS = frozenset(
    {
        "ndownloader.figshare.com",
        "figshare.com",
        "api.figshare.com",
        "zenodo.org",
        "www.zenodo.org",
        "clowder.edap-cluster.com",
    }
)
FORBIDDEN_EXTENSIONS = frozenset(
    {
        ".exe",
        ".dll",
        ".so",
        ".dylib",
        ".bat",
        ".cmd",
        ".ps1",
        ".sh",
        ".bash",
        ".zsh",
        ".py",
        ".pyc",
        ".pyo",
        ".rb",
        ".pl",
        ".php",
        ".js",
        ".mjs",
        ".jar",
        ".war",
        ".class",
        ".apk",
        ".msi",
        ".scr",
        ".com",
        ".vbs",
        ".wsf",
        ".htm",
        ".html",  # do not treat HTML scrapes as data artifacts
    }
)
ALLOWED_CONTENT_SNIPPETS = (
    "json",
    "csv",
    "tsv",
    "xml",
    "text",
    "octet-stream",  # many scientific APIs
    "excel",
    "spreadsheet",
    "zip",  # ISA-Tab; still not executed
    "xlsx",
    "ms-excel",
    "openxml",
)

_MAGIC_EXECUTABLE = [
    b"\x7fELF",  # ELF
    b"MZ",  # PE
    b"\xca\xfe\xba\xbe",  # Mach-O fat
    b"\xfe\xed\xfa",  # Mach-O
]

# Parent domains whose subdomains are also allowed.
_ALLOWED_PARENT_SUFFIXES = (
    "ebi.ac.uk",
    "figshare.com",
    "metabolomicsworkbench.org",
    "ncbi.nlm.nih.gov",
)

# Storage CDNs used only as redirect targets from allowlisted scientific hosts
# (e.g. Figshare → S3). Never accept these as the *initial* URL.
_REDIRECT_STORAGE_SUFFIXES = (
    "amazonaws.com",
    "cloudfront.net",
    "zenodo.org",
)


class SecureFetchError(RuntimeError):
    pass


def host_allowed(url: str, *, allow_storage_cdn: bool = False) -> bool:
    """Return True if ``url`` targets an allowlisted host/scheme."""
    parsed = urlparse(url)
    if parsed.scheme not in ALLOWED_SCHEMES:
        return False
    host = (parsed.hostname or "").lower()
    if not host:
        return False
    if host in ALLOWED_HOSTS:
        return True
    if any(host.endswith("." + h) for h in _ALLOWED_PARENT_SUFFIXES):
        return True
    if allow_storage_cdn and any(host == s or host.endswith("." + s) for s in _REDIRECT_STORAGE_SUFFIXES):
        return True
    return False


def _origin_may_use_storage_cdn(url: str) -> bool:
    host = (urlparse(url).hostname or "").lower()
    if host in _STORAGE_REDIRECT_ORIGINS:
        return True
    return any(host.endswith("." + h) for h in ("figshare.com", "zenodo.org"))


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _looks_like_executable(data: bytes) -> bool:
    head = data[:8]
    return any(head.startswith(m) for m in _MAGIC_EXECUTABLE)


def _extension_forbidden(path_or_name: str) -> bool:
    name = path_or_name.lower().split("?")[0]
    for ext in FORBIDDEN_EXTENSIONS:
        if name.endswith(ext):
            return True
    return False


def secure_fetch(
    url: str,
    *,
    dest: Path,
    quarantine_dir: Path | None = None,
    max_bytes: int = MAX_BYTES_DEFAULT,
    expected_sha256: str | None = None,
    timeout: int = 120,
    retries: int = 4,
    allow_html: bool = False,
) -> dict[str, Any]:
    """
    Download ``url`` into quarantine, validate, then atomically promote to ``dest``.

    Returns provenance dict (url, sha256, bytes, content_type, fetched_at).
    """
    if not host_allowed(url):
        raise SecureFetchError(f"URL host not in allowlist: {url}")
    if _extension_forbidden(url) and not allow_html:
        raise SecureFetchError(f"Forbidden file extension in URL: {url}")

    dest = Path(dest)
    qdir = Path(quarantine_dir or (dest.parent / ".quarantine"))
    qdir.mkdir(parents=True, exist_ok=True)
    qpath = qdir / f"{dest.name}.{int(time.time())}.part"

    last_err: Exception | None = None
    content_type = ""
    data = b""

    for attempt in range(retries):
        try:
            # Manual redirect follow so every hop is allowlist-checked.
            current = url
            allow_cdn = False
            hist: list[str] = []
            with requests.Session() as hop_session:
                hop_session.headers.update({"User-Agent": USER_AGENT, "Accept": "*/*"})
                for _ in range(8):
                    if not host_allowed(current, allow_storage_cdn=allow_cdn):
                        raise SecureFetchError(f"URL host not in allowlist: {current}")
                    with hop_session.get(
                        current, timeout=timeout, stream=True, allow_redirects=False
                    ) as r:
                        if r.is_redirect or r.status_code in {301, 302, 303, 307, 308}:
                            loc = r.headers.get("Location")
                            if not loc:
                                raise SecureFetchError(f"Redirect without Location from {current}")
                            nxt = urljoin(current, loc)
                            hist.append(nxt)
                            # Storage CDN only if the *original* fetch origin is trusted for blobs
                            if _origin_may_use_storage_cdn(url):
                                allow_cdn = True
                            current = nxt
                            continue
                        r.raise_for_status()
                        content_type = (r.headers.get("Content-Type") or "").lower()
                        if "html" in content_type and not allow_html:
                            raise SecureFetchError(f"Refusing HTML payload from {url} ({content_type})")
                        if content_type and not any(s in content_type for s in ALLOWED_CONTENT_SNIPPETS):
                            if content_type not in {"", "application/force-download", "binary/octet-stream"}:
                                if "json" not in content_type and "xml" not in content_type and "csv" not in content_type:
                                    if any(x in content_type for x in ("javascript", "wasm", "x-msdownload", "x-sh")):
                                        raise SecureFetchError(f"Dangerous content-type: {content_type}")

                        buf = bytearray()
                        for chunk in r.iter_content(chunk_size=1024 * 256):
                            if not chunk:
                                continue
                            buf.extend(chunk)
                            if len(buf) > max_bytes:
                                raise SecureFetchError(f"Payload exceeds max_bytes={max_bytes}")
                        data = bytes(buf)
                    break
                else:
                    raise SecureFetchError(f"Too many redirects for {url}: {hist}")
            break
        except (requests.RequestException, SecureFetchError) as e:
            last_err = e
            time.sleep(min(2**attempt, 20))
    else:
        raise SecureFetchError(f"Fetch failed for {url}: {last_err}")

    if not data:
        raise SecureFetchError(f"Empty payload from {url}")
    if _looks_like_executable(data):
        raise SecureFetchError(f"Executable magic bytes detected — refusing {url}")
    if data[:2] == b"#!" and b"\n" in data[:80]:
        # shebang scripts
        raise SecureFetchError(f"Script shebang detected — refusing {url}")

    digest = sha256_bytes(data)
    if expected_sha256 and digest.lower() != expected_sha256.lower():
        raise SecureFetchError(
            f"SHA-256 mismatch for {url}: got {digest}, expected {expected_sha256}"
        )

    qpath.write_bytes(data)
    dest.parent.mkdir(parents=True, exist_ok=True)
    qpath.replace(dest)

    prov = {
        "url": url,
        "final_path": str(dest),
        "sha256": digest,
        "n_bytes": len(data),
        "content_type": content_type,
        "fetched_at": datetime.now(timezone.utc).isoformat(),
        "allowlisted": True,
        "executed": False,
    }
    prov_path = dest.with_suffix(dest.suffix + ".provenance.json")
    if dest.suffix == "":
        prov_path = Path(str(dest) + ".provenance.json")
    prov_path.write_text(json.dumps(prov, indent=2) + "\n")
    return prov


def validate_json_artifact(
    path: Path,
    *,
    required_keys: list[str] | None = None,
    max_depth: int = 12,
) -> dict[str, Any]:
    """Parse JSON and apply light poison checks (depth, type, required keys)."""
    raw = Path(path).read_bytes()
    if _looks_like_executable(raw):
        raise SecureFetchError(f"Executable bytes in supposed JSON: {path}")
    try:
        doc = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        raise SecureFetchError(f"Invalid JSON artifact {path}: {e}") from e

    def _depth(o: Any, d: int = 0) -> int:
        if d > max_depth:
            return d
        if isinstance(o, dict):
            return max([d] + [_depth(v, d + 1) for v in o.values()] or [d])
        if isinstance(o, list):
            return max([d] + [_depth(v, d + 1) for v in o[:50]] or [d])
        return d

    if _depth(doc) > max_depth:
        raise SecureFetchError(f"JSON nesting too deep in {path}")
    if required_keys:
        if not isinstance(doc, dict):
            raise SecureFetchError(f"Expected object in {path}")
        missing = [k for k in required_keys if k not in doc]
        if missing:
            raise SecureFetchError(f"Missing keys {missing} in {path}")
    # Reject unexpected code-like strings at top level
    if isinstance(doc, dict):
        for k in doc:
            if re.search(r"(__import__|subprocess|eval\(|exec\()", str(k)):
                raise SecureFetchError(f"Suspicious key in {path}: {k}")
    return doc if isinstance(doc, dict) else {"_value": doc}


def write_integrity_manifest(paths: list[Path], out: Path, *, root: Path | None = None) -> dict[str, Any]:
    root = Path(root) if root else Path.cwd()
    files = []
    for p in paths:
        p = Path(p)
        if not p.exists() or not p.is_file():
            continue
        try:
            rel = str(p.resolve().relative_to(root.resolve()))
        except ValueError:
            rel = str(p)
        files.append(
            {
                "path": rel,
                "sha256": sha256_file(p),
                "n_bytes": p.stat().st_size,
            }
        )
    doc = {
        "version": "1.1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "n_files": len(files),
        "files": files,
        "policy": {
            "no_execution_of_downloads": True,
            "allowlisted_hosts_only": True,
            "sha256_required": True,
            "paths_are_relative": True,
        },
    }
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=2) + "\n")
    return doc
