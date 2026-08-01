# Architecture

```
Query → 16 StackModel.predict() → ModelOutput[]
      → fuse_vocs (weighted log2fc + RRF + sign agreement)
      → fuse_aspects (genes/pathways/microbes/cells/…)
      → STACK_REPORT.html + dashboard + JSON
```

ExhalePath hybrid runs first and caches `biomarker_report` for calibrator, Census, OPERA, and PBPK heads.
