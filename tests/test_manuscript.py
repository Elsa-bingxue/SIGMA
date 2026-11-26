import json

from sigma_spatial.manuscript import (
    FIGURE_BUILDERS, MANUSCRIPT_PROGRAMS, manuscript_manifest,
    write_manuscript_manifest,
)


def test_manuscript_manifest_freezes_current_programs():
    assert MANUSCRIPT_PROGRAMS["HPD_A1"] == 3
    assert MANUSCRIPT_PROGRAMS["HPD_B1"] == 2
    assert MANUSCRIPT_PROGRAMS["HPD_C1"] == 1
    assert MANUSCRIPT_PROGRAMS["HBC515"] == 3
    assert MANUSCRIPT_PROGRAMS["GBM"] == 1
    assert set(FIGURE_BUILDERS) == {f"Fig{x}" for x in range(2, 9)}


def test_write_manuscript_manifest(tmp_path):
    path = write_manuscript_manifest(tmp_path / "manifest.json")
    payload = json.loads(path.read_text())
    assert payload == manuscript_manifest()
    assert payload["package_version"] == "0.3.2"
