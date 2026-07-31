"""Visual dashboard pack for one-line voc usage."""

from __future__ import annotations

from pathlib import Path

from exhalepath.biomarker import ExhaleBiomarkerEngine
from exhalepath.viz.dashboard import (
    default_report_dir,
    print_comprehensive_console,
    save_visual_dashboard,
)


def test_default_report_dir_slug():
    p = default_report_dir("Maple syrup urine disease", "systemic")
    assert p.parts[0] == "runs"
    assert "maple_syrup" in p.name
    assert "systemic" in p.name


def test_save_visual_dashboard_writes_pack(tmp_path: Path):
    engine = ExhaleBiomarkerEngine(use_opentargets=False, reload_knowledge=True)
    report = engine.predict("depression", location="brain", top_n=10, explain=True)
    out = tmp_path / "pack"
    paths = save_visual_dashboard(report, out)
    assert paths["report_html"].is_file()
    assert paths["dashboard"].is_file()
    assert paths["report_md"].is_file()
    assert paths["manifest"].is_file()
    html = paths["report_html"].read_text()
    assert report.disease_name.lower() in html.lower()
    assert "data:image/png" in html
    assert paths["dashboard"].stat().st_size > 1000
    print_comprehensive_console(report, top_display=5)
