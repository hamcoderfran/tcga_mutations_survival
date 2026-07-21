"""Optional very-light LLM backends for QuerySlots extraction (Ollama / OpenAI-compatible)."""

from __future__ import annotations

import json
import os
import re
from typing import Any

import requests

from .rules import parse_rules
from .slots import QuerySlots

_SYSTEM = """You extract structured fields for an exhaled VOC biomarker tool.
Return ONLY compact JSON with keys:
disease (string), comorbidities (string array), location (string|null),
age_years (number|null), sex ("male"|"female"|"other"|null),
stage (string|null), genes (string array), smoking ("never"|"former"|"current"|null),
mode ("hybrid"|"physiology"|"legacy"), top (int).
Use atlas disease names when possible (e.g. depression, schizophrenia, lung adenocarcinoma, obesity, heart_disease).
No markdown, no commentary."""


def _coerce_slots(data: dict[str, Any], *, text: str, method: str) -> QuerySlots:
    sex = data.get("sex")
    if sex not in {"male", "female", "other", None}:
        sex = None
    smoking = data.get("smoking")
    if smoking not in {"never", "former", "current", None}:
        smoking = None
    mode = data.get("mode") or "hybrid"
    if mode not in {"physiology", "hybrid", "legacy"}:
        mode = "hybrid"
    genes = data.get("genes") or []
    if isinstance(genes, str):
        genes = [g.strip() for g in genes.split(",") if g.strip()]
    comorb = data.get("comorbidities") or []
    if isinstance(comorb, str):
        comorb = [c.strip() for c in comorb.split(",") if c.strip()]
    age = data.get("age_years")
    try:
        age_f = float(age) if age is not None else None
    except (TypeError, ValueError):
        age_f = None
    top = data.get("top") or 20
    try:
        top_i = int(top)
    except (TypeError, ValueError):
        top_i = 20
    return QuerySlots(
        disease=data.get("disease"),
        comorbidities=[str(c) for c in comorb],
        location=data.get("location"),
        age_years=age_f,
        sex=sex,
        stage=data.get("stage"),
        genes=[str(g).upper() for g in genes],
        smoking=smoking,
        mode=mode,  # type: ignore[arg-type]
        top=top_i,
        source_text=text,
        parse_method=method,
    )


def _extract_json(text: str) -> dict[str, Any]:
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    # Find first {...}
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        raise ValueError("no JSON object in LLM response")
    return json.loads(m.group(0))


def parse_ollama(
    text: str,
    *,
    model: str | None = None,
    host: str | None = None,
    timeout: float = 60.0,
) -> QuerySlots:
    host = (host or os.environ.get("VOC_OLLAMA_HOST") or "http://127.0.0.1:11434").rstrip(
        "/"
    )
    model = model or os.environ.get("VOC_OLLAMA_MODEL") or "qwen2.5:0.5b"
    payload = {
        "model": model,
        "stream": False,
        "format": "json",
        "options": {"temperature": 0, "num_predict": 256},
        "messages": [
            {"role": "system", "content": _SYSTEM},
            {"role": "user", "content": text},
        ],
    }
    r = requests.post(f"{host}/api/chat", json=payload, timeout=timeout)
    r.raise_for_status()
    content = (r.json().get("message") or {}).get("content") or ""
    data = _extract_json(content)
    return _coerce_slots(data, text=text, method=f"llm:ollama:{model}")


def parse_openai_compatible(
    text: str,
    *,
    model: str | None = None,
    base_url: str | None = None,
    api_key: str | None = None,
    timeout: float = 60.0,
) -> QuerySlots:
    base_url = (
        base_url
        or os.environ.get("VOC_OPENAI_BASE_URL")
        or os.environ.get("OPENAI_BASE_URL")
        or "https://api.openai.com/v1"
    ).rstrip("/")
    api_key = api_key or os.environ.get("VOC_OPENAI_API_KEY") or os.environ.get(
        "OPENAI_API_KEY"
    )
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY / VOC_OPENAI_API_KEY not set")
    model = model or os.environ.get("VOC_OPENAI_MODEL") or "gpt-4o-mini"
    headers = {"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"}
    payload = {
        "model": model,
        "temperature": 0,
        "response_format": {"type": "json_object"},
        "messages": [
            {"role": "system", "content": _SYSTEM},
            {"role": "user", "content": text},
        ],
    }
    r = requests.post(
        f"{base_url}/chat/completions", headers=headers, json=payload, timeout=timeout
    )
    r.raise_for_status()
    content = r.json()["choices"][0]["message"]["content"]
    data = _extract_json(content)
    return _coerce_slots(data, text=text, method=f"llm:openai:{model}")


def parse_with_optional_llm(
    text: str,
    *,
    llm: str = "auto",
    disease_catalog: list[dict[str, Any]] | None = None,
) -> QuerySlots:
    """
    llm:
      - off|rules|none → rules only
      - ollama → local Ollama tiny model
      - openai → OpenAI-compatible API
      - auto → try ollama if VOC_LLM=ollama or daemon up; else openai if key; else rules
    """
    llm = (llm or "auto").lower().strip()
    env = (os.environ.get("VOC_LLM") or "").lower().strip()
    if llm == "auto" and env:
        llm = env

    if llm in {"off", "rules", "none", "false", "0"}:
        return parse_rules(text, disease_catalog=disease_catalog)

    if llm == "ollama":
        try:
            return parse_ollama(text)
        except Exception:
            return parse_rules(text, disease_catalog=disease_catalog)

    if llm in {"openai", "openai-compatible"}:
        try:
            return parse_openai_compatible(text)
        except Exception:
            return parse_rules(text, disease_catalog=disease_catalog)

    # auto
    # 1) Prefer explicit Ollama if reachable
    host = (os.environ.get("VOC_OLLAMA_HOST") or "http://127.0.0.1:11434").rstrip("/")
    try:
        ping = requests.get(f"{host}/api/tags", timeout=1.5)
        if ping.ok:
            try:
                return parse_ollama(text)
            except Exception:
                pass
    except Exception:
        pass
    # 2) OpenAI key
    if os.environ.get("OPENAI_API_KEY") or os.environ.get("VOC_OPENAI_API_KEY"):
        try:
            return parse_openai_compatible(text)
        except Exception:
            pass
    # 3) Rules fallback (always works, zero deps beyond atlas)
    return parse_rules(text, disease_catalog=disease_catalog)
