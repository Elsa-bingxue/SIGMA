import anndata as ad
import numpy as np
import pandas as pd
import pytest

from sigma_spatial.programs import discover_programs


def example():
    rng = np.random.default_rng(12)
    distance = np.linspace(-10, 10, 240)
    shapes = [
        np.exp(-(distance / 2) ** 2),
        1 / (1 + np.exp(-distance)),
        1 / (1 + np.exp(distance)),
        distance / 10,
    ]
    columns = []
    for shape in shapes:
        for _ in range(5):
            columns.append(shape + rng.normal(0, .04, len(distance)))
    obj = ad.AnnData(np.column_stack(columns))
    obj.obs["sigma_d_signed"] = distance
    ranking = pd.DataFrame({
        "j": np.arange(20), "mz": 100 + np.arange(20),
        "candidate_score": np.linspace(1, 0, 20),
    })
    return obj, ranking


def test_program_discovery_is_deterministic_and_complete():
    obj, ranking = example()
    first = discover_programs(obj, ranking, n_bins=20, min_bin_points=3)
    second = discover_programs(obj, ranking, n_bins=20, min_bin_points=3)
    assert first.assignments.cluster.tolist() == second.assignments.cluster.tolist()
    assert len(first.assignments) == 20
    assert set(first.scores) == {0, 1, 2, 3}
    assert all(len(value) == obj.n_obs for value in first.scores.values())
    assert first.centroid_profiles.groupby("program").size().eq(20).all()


def test_reference_mode_preserves_final_table_order():
    obj, ranking = example()
    shuffled = ranking.sample(frac=1, random_state=3).reset_index(drop=True)
    result = discover_programs(
        obj, shuffled, selection_mode="reference", top_k=8,
        n_programs=2, n_bins=20, min_bin_points=3,
    )
    assert result.assignments.j.tolist() == shuffled.head(8).j.tolist()
    assert result.config["ranking_col"] is None


def test_standardized_mode_requires_an_explicit_or_supported_ranking_score():
    obj, ranking = example()
    with pytest.raises(KeyError, match="ranking_col"):
        discover_programs(
            obj, ranking.drop(columns="candidate_score"), n_programs=2,
            n_bins=20, min_bin_points=3,
        )


def test_joint_matrix_uses_sm_local_feature_indices():
    obj, ranking = example()
    joint = ad.AnnData(np.c_[np.full((obj.n_obs, 2), 99.0), obj.X])
    joint.obs["sigma_d_signed"] = obj.obs["sigma_d_signed"].to_numpy()
    joint.var["feature_type"] = ["ST", "ST"] + ["SM"] * obj.n_vars
    result = discover_programs(joint, ranking, n_bins=20, min_bin_points=3)
    assert len(result.assignments) == 20
