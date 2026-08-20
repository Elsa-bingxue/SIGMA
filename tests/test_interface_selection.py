import anndata as ad
import numpy as np
import pandas as pd

from sigma_spatial.interface_selection import (
    _mz_value, rank_interface_metabolites, select_representative_metabolites,
    summarize_interface_programs,
)
from sigma_spatial.selection import side_near_far_statistics


def _transition_example():
    rng = np.random.default_rng(12)
    n = 400
    d = np.linspace(-10, 10, n)
    xy = np.column_stack([np.arange(n), np.zeros(n)])
    interface = 3 * np.exp(-np.abs(d)) + rng.normal(0, .08, n)
    region_only = 2 * (d < 0) + rng.normal(0, .08, n)
    sparse = np.zeros(n); sparse[rng.choice(n, 8, replace=False)] = 10
    noise = rng.normal(1, .5, n)
    x = np.column_stack([interface, region_only, sparse, noise])
    obj = ad.AnnData(np.maximum(x, 0))
    obj.var_names = ["100.1", "200.2", "300.3", "400.4"]
    obj.var["feature_type"] = "SM"
    obj.obsm["spatial"] = xy
    obj.obs["sigma_d_signed"] = d
    obj.obs["sigma_inside"] = d < 0
    return obj


def test_rank_interface_metabolites_separates_interface_from_region_and_sparse():
    table = rank_interface_metabolites(_transition_example())
    interface = table[table.j.eq(0)].iloc[0]
    region = table[table.j.eq(1)].iloc[0]
    sparse = table[table.j.eq(2)].iloc[0]
    assert interface.interface_effect > region.interface_effect
    assert interface.region_adjusted_delta_r2 > region.region_adjusted_delta_r2
    assert bool(interface.gradient_near_boundary)
    assert sparse.detect_rate < .05
    assert not bool(sparse.eligible_representative)


def test_program_summary_and_representatives_use_feature_qc():
    obj = _transition_example()
    ranking = rank_interface_metabolites(obj)
    assignments = pd.DataFrame({
        "sigma_sm_index": [0, 1, 2, 3], "mz": [100.1, 200.2, 300.3, 400.4],
        "cluster": [0, 1, 0, 1],
    })
    scores = {0: obj.X[:, 0], 1: obj.X[:, 1]}
    stats = side_near_far_statistics(scores, obj.obs.sigma_d_signed.to_numpy())
    summary, merged = summarize_interface_programs(assignments, ranking, stats)
    assert summary.iloc[0].program == 0
    reps = select_representative_metabolites(assignments, ranking, program=0, n=1)
    assert int(reps.iloc[0].sigma_sm_index) == 0
    assert reps.iloc[0].representative_qc == "pass"


def test_common_prefixed_mz_names_are_parsed():
    assert _mz_value("mz_585.49561") == 585.49561
    assert _mz_value("SM_585.49561") == 585.49561
    assert _mz_value("m/z 585.49561") == 585.49561
