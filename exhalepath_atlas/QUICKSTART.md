# Quick start — install & run in 30 seconds

```bash
pip install "voc-breath @ git+https://github.com/hamcoderfran/tcga_mutations_survival.git#subdirectory=exhalepath_atlas"
```

Or from a downloaded wheel:

```bash
pip install voc_breath-*.whl
```

Then:

```bash
# Structured flags (disease-first)
voc "depression" -l brain -c obesity --age 24 --sex male
voc "schizophrenia" -l brain -c heart_disease --age 18 --sex male
voc "lung adenocarcinoma" -l "left lower lobe" --stage II --genes KRAS,TP53

# Interactive question fields (disease, comorbidities, location, age, …)
voc ask

# Natural language → slots (optional very-light LLM; else rules; else ask)
voc nl "24yo obese male with depression"
voc nl "stage II LUAD left lower lobe KRAS TP53" --llm rules --yes
voc ask --nl "depression with obesity" --llm auto

voc list-diseases
voc --help
```

**NL backends** (`--llm`):
- `auto` (default) — try local Ollama (`qwen2.5:0.5b`) → OpenAI-compatible API → rules
- `rules` — zero-weight pattern parser (always works)
- `ollama` — `VOC_OLLAMA_HOST` / `VOC_OLLAMA_MODEL` (default `qwen2.5:0.5b`)
- `openai` — `OPENAI_API_KEY` or `VOC_OPENAI_API_KEY`

No subcommand needed for structured predictions — `voc "<disease>" …` runs the biomarker engine.

Research / hypothesis tool only — not a medical device.
