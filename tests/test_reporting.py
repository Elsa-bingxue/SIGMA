import numpy as np
import pandas as pd
import anndata as ad

from sigma_spatial.reporting import (
    anisotropy_report, disease_profile_representatives,
    program_report_data, representative_metabolites,
    stroma_near_far_statistics,
)
from sigma_spatial.preprocessing import resolve_sm_matrix


def _example():
    rng = np.random.default_rng(4)
    xy = rng.normal(size=(240, 2))
    x = rng.gamma(2, 1, size=(240, 8))
    obj = ad.AnnData(x)
    obj.obsm["spatial"] = xy
    obj.obs["sigma_d_signed"] = xy[:, 0]
    assignments = pd.DataFrame({
        "j": np.arange(8), "mz": 100 + np.arange(8),
        "cluster": np.repeat(np.arange(4), 2), "r2_logI": np.linspace(.1, .8, 8),
    })
    return obj, assignments


def test_frozen_program_report_does_not_recluster():
    obj, assignments = _example()
    report = program_report_data(obj, assignments, n_bins=10, min_bin_points=1)
    assert sorted(report["scores"]) == [0, 1, 2, 3]
    assert all(score.shape == (obj.n_obs,) for score in report["scores"].values())
    assert report["assignments"]["cluster"].tolist() == assignments["cluster"].tolist()


def test_representatives_respect_requested_cluster():
    _, assignments = _example()
    selected = representative_metabolites(assignments, cluster=3, n=2)
    assert selected["cluster"].eq(3).all()
    assert selected["r2_logI"].is_monotonic_decreasing


def test_disease_representatives_preserve_original_profile_ranking():
    assignments = pd.DataFrame({
        "j": np.arange(7), "mz": 100 + np.arange(7),
        "cluster": [0, 0, 0, 0, 0, 1, 1],
        "r2_logI": [.31, .28, .24, .20, .19, .90, .80],
        "slope": [-.1, -.2, -.3, -.4, -.5, -.1, -.1],
        "lambda": [10, 20, 30, 40, 10000, 10, 20],
        "boundary_enrichment_ratio": [1.1, 1.2, 1.3, 1.4, 2.0, 2.0, 2.0],
    })
    selected = disease_profile_representatives(assignments, cluster=0, n=4)
    assert selected["cluster"].eq(0).all()
    assert selected["mz"].tolist() == [100, 101, 102, 103]
    assert selected["representative_selection"].eq(
        "original_hpd_distance_profile"
    ).all()


def test_region_defined_disease_uses_original_r2_ranking():
    from sigma_spatial.lambda_profile import get_lambda_profile_preset

    assert get_lambda_profile_preset("region_defined_disease").ranking_col == "r2_logI"


def test_anisotropy_report_has_all_sectors():
    obj, _ = _example()
    report = anisotropy_report(obj, 0, min_points=1)
    assert report["table"]["sector"].tolist() == list(range(8))
    assert "anisotropy_cv" in report


def test_near_far_selects_program_from_metabolic_scores_only():
    obj, _ = _example()
    obj.obs["annotation"] = "Stroma"
    # Ensure both preserved distance windows contain observations.
    obj.obs["sigma_d_signed"] = np.linspace(1, 220, obj.n_obs)
    near_signal = -obj.obs["sigma_d_signed"].to_numpy()
    programs = {0: np.zeros(obj.n_obs), 1: near_signal, 2: near_signal * .1, 3: near_signal * 2}
    signatures = {"plasma_cell_score": near_signal * 100}
    table, leading = stroma_near_far_statistics(
        obj, programs, signature_scores=signatures, spot_spacing=10, near_width=5, far_start=8
    )
    assert leading == 3
    assert "plasma_cell_score" in table["score_col"].tolist()


def test_joint_st_sm_reporting_uses_sm_local_indices():
    n = 40
    # First two columns are ST and deliberately have unrelated large values.
    x = np.column_stack([
        np.full(n, 1000.0), np.full(n, 2000.0),
        np.linspace(0, 1, n), np.linspace(1, 0, n),
    ])
    obj = ad.AnnData(x)
    obj.var["feature_type"] = ["ST", "ST", "SM", "SM"]
    obj.obsm["spatial"] = np.column_stack([np.arange(n), np.zeros(n)])
    obj.obs["sigma_d_signed"] = np.linspace(-2, 2, n)
    assignments = pd.DataFrame({
        "sigma_sm_index": [0, 1], "mz": [100.1, 200.2], "cluster": [0, 1],
    })
    report = program_report_data(obj, assignments, n_bins=8, min_bin_points=1)
    assert report["scores"][0][0] < report["scores"][0][-1]
    assert report["scores"][1][0] > report["scores"][1][-1]


def test_explicit_matrix_source_does_not_silently_prefer_raw():
    obj = ad.AnnData(np.array([[1., 2.], [3., 4.]]))
    obj.layers["raw"] = np.array([[10., 20.], [30., 40.]])
    x_matrix, _ = resolve_sm_matrix(obj, "X")
    raw_matrix, _ = resolve_sm_matrix(obj, "raw")
    np.testing.assert_array_equal(x_matrix, obj.X)
    np.testing.assert_array_equal(raw_matrix, obj.layers["raw"])


def test_sm_var_source_maps_joint_object_by_feature_name():
    x = np.array([[10., 20., 1., 2.], [30., 40., 3., 4.]])
    obj = ad.AnnData(x)
    obj.var_names = ["ST_A", "ST_B", "SM_100", "SM_200"]
    obj.uns["SM_features"] = np.array(["SM_100", "SM_200"])
    obj.uns["SM_mz"] = np.array([100.1, 200.2])
    matrix, names = resolve_sm_matrix(obj, "sm_var")
    np.testing.assert_array_equal(matrix, x[:, 2:])
    assert names.tolist() == ["100.1", "200.2"]
