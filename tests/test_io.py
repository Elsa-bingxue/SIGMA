import pandas as pd
import pytest

from sigma_spatial.io import load_spatial_metabolomics


def test_csv_loader_aligns_rows_by_spot_id(tmp_path):
    matrix = pd.DataFrame({"mz1": [1, 2], "mz2": [3, 4]}, index=["s1", "s2"])
    coords = pd.DataFrame({"x": [20, 10], "y": [2, 1]}, index=["s2", "s1"])
    anchors = pd.DataFrame({"anchor": [0, 1]}, index=["s2", "s1"])
    mp, cp, ap = tmp_path/"m.csv", tmp_path/"c.csv", tmp_path/"a.csv"
    matrix.to_csv(mp); coords.to_csv(cp); anchors.to_csv(ap)
    a = load_spatial_metabolomics(mp, cp, ap)
    assert a.obs_names.tolist() == ["s1", "s2"]
    assert a.obsm["spatial"].tolist() == [[10.0, 1.0], [20.0, 2.0]]
    assert a.obs["sigma_anchor"].tolist() == [1, 0]


def test_csv_loader_rejects_mismatched_spots(tmp_path):
    matrix = pd.DataFrame({"mz1": [1]}, index=["s1"])
    coords = pd.DataFrame({"x": [1], "y": [2]}, index=["s2"])
    mp, cp = tmp_path/"m.csv", tmp_path/"c.csv"
    matrix.to_csv(mp); coords.to_csv(cp)
    with pytest.raises(ValueError, match="do not match"):
        load_spatial_metabolomics(mp, cp)
