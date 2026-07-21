"""Packaging / easy-install smoke tests."""

from pathlib import Path

from exhalepath.config import DATA_DIR, KNOWLEDGE_DIR
from exhalepath import __version__


def test_packaged_or_repo_knowledge_resolves():
    assert (KNOWLEDGE_DIR / "voc_catalog.json").exists()
    assert (KNOWLEDGE_DIR / "disease_voc_priors.json").exists()
    assert DATA_DIR.exists()


def test_version_is_set():
    assert __version__
    parts = __version__.split(".")
    assert len(parts) >= 2


def test_disease_first_cli_dispatch(monkeypatch):
    import sys
    from exhalepath import cli

    called = {}

    def fake_app():
        called["argv"] = list(sys.argv)

    monkeypatch.setattr(cli, "app", fake_app)
    monkeypatch.setattr(sys, "argv", ["voc", "depression", "-l", "brain", "-c", "obesity"])
    cli.main(None)
    # main() rewrites sys.argv when disease-first
    assert called["argv"][1] == "biomarker"
    assert called["argv"][2] == "depression"
