from __future__ import annotations

import json
import time
from typing import Any

import requests

from ..config import USER_AGENT

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
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, json.JSONDecodeError) as e:
            last_err = e
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"HTTP {method} {url} failed after {retries} retries: {last_err}")
