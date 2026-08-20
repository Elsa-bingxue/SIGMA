import anndata as ad
import numpy as np

from sigma_spatial import get_lambda_profile_preset, rank_lambda_profiles


def test_lambda_profile_ranks_boundary_decay_without_frozen_results():
    distance = np.linspace(-10, 10, 400)
    decay = 5 * np.exp(-np.abs(distance) / 2)
    region = 2 * (distance < 0)
    obj = ad.AnnData(np.column_stack([decay, region]))
    obj.var_names = ["mz_100.1", "mz_200.2"]
    obj.var["feature_type"] = "SM"
    obj.obs["sigma_d_signed"] = distance
    table = rank_lambda_profiles(obj, var_top=2, min_valid=100)
    assert int(table.iloc[0].j) == 0
    assert table.iloc[0].mz == 100.1
    assert table.iloc[0]["lambda"] > 0
    assert table.iloc[0].r2_logI > .05


def test_workflow_presets_are_shared_not_slice_specific():
    hcc = get_lambda_profile_preset("transferred_pathology")
    disease = get_lambda_profile_preset("region_defined_disease")
    assert hcc.already_log is True
    assert hcc.min_r2 == .005
    assert disease.ranking_col == "interface_score"
    assert disease.min_detect_rate == .01
