import json

import pytest

from sigma_spatial.config import load_analysis_config


def test_load_analysis_config_resolves_relative_paths(tmp_path):
    path = tmp_path / "run.json"
    path.write_text(json.dumps({
        "schema_version": 1,
        "input_h5ad": "data/example.h5ad",
        "output_dir": "results/example",
        "workflow": "direct_pathology",
        "random_state": 0,
    }))
    config = load_analysis_config(path)
    assert config.input_h5ad == tmp_path / "data/example.h5ad"
    assert config.output_dir == tmp_path / "results/example"
    assert config.options["random_state"] == 0
    assert config.options["report_level"] == "complete"


def test_load_analysis_config_rejects_unknown_keys(tmp_path):
    path = tmp_path / "run.json"
    path.write_text(json.dumps({
        "input_h5ad": "example.h5ad", "output_dir": "results", "guess": True,
    }))
    with pytest.raises(ValueError, match="unknown analysis configuration keys"):
        load_analysis_config(path)
