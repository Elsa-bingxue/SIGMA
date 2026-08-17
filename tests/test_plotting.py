import anndata as ad
import matplotlib
import numpy as np

matplotlib.use("Agg")

from sigma_spatial.plotting import (  # noqa: E402
    plot_boundary, plot_interface_feature, plot_region_probability, plot_signed_distance,
)


def plotted_adata():
    a = ad.AnnData(np.arange(12, dtype=float).reshape(4, 3))
    a.var_names = ["mz1", "mz2", "mz3"]
    a.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    a.obs["sigma_region_probability"] = [0, .25, .75, 1]
    a.obs["sigma_boundary"] = [False, True, True, False]
    a.obs["sigma_d_signed"] = [-1, 0, 0, 1]
    return a


def test_public_plots_do_not_modify_anndata():
    a = plotted_adata(); before = a.copy()
    axes = [plot_region_probability(a), plot_boundary(a), plot_signed_distance(a),
            plot_interface_feature(a, "mz2")]
    assert all(ax is not None for ax in axes)
    np.testing.assert_array_equal(a.X, before.X)
    assert a.obs.equals(before.obs)
