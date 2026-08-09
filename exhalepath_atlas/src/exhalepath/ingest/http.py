from __future__ import annotations

import json
import time
from typing import Any

import requests

from ..config import USER_AGENT
from .secure_fetch import SecureFetchError, host_allowed

_session = requests.Session()
_session.headers.update({"User-Agent": USER_AGENT})


def request_json(
    method: str,
    url: str,
    *,
    params: dict[str, Any] | None = None,
    json_body: dict[str, Any] | None = None,
    retries: int = 5,
    timeout: int = 90,
) -> dict[str, Any]:
    """JSON API helper with host allowlist (SSRF / supply-chain guard)."""
    if not host_allowed(url):
        raise SecureFetchError(f"URL host not in allowlist: {url}")
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            r = _session.request(
                method,
                url,
                params=params,
                json=json_body,
                timeout=timeout,
            )
            # Reject redirects off allowlist (storage CDNs OK as redirect targets)
            if not host_allowed(r.url, allow_storage_cdn=True):
                raise SecureFetchError(f"Redirected to non-allowlisted host: {r.url}")
            r.raise_for_status()
            ct = (r.headers.get("Content-Type") or "").lower()
            if any(x in ct for x in ("javascript", "wasm", "x-msdownload", "x-sh")):
                raise SecureFetchError(f"Dangerous content-type from API: {ct}")
            raw = r.content
            if raw[:4] == b"\x7fELF" or raw[:2] == b"MZ":
                raise SecureFetchError("Executable payload disguised as JSON API response")
            return r.json()
        except (requests.RequestException, json.JSONDecodeError, SecureFetchError) as e:
            last_err = e
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"HTTP {method} {url} failed after {retries} retries: {last_err}")
