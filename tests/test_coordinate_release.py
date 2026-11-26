"""Release regressions for coordinate integration and near/far validity."""
import anndata as ad
import numpy as np
import pytest

from sigma_spatial import run_sigma, run_sigma_weak_anchor
from sigma_spatial.interface_selection import rank_interface_metabolites
from sigma_spatial.reporting import stroma_near_far_statistics


def sample():
    rng = np.random.default_rng(31)
    x, y = np.meshgrid(np.arange(6.), np.arange(6.))
    a = ad.AnnData(rng.gamma(2., 1., size=(36, 8)))
    a.obsm['spatial'] = np.c_[x.ravel(), y.ravel()]
    a.obsm['X_harmony'] = rng.normal(size=(36, 3))
    a.obs['sigma_anchor'] = np.where(x.ravel() <= 1, 0., np.where(x.ravel() >= 4, 1., np.nan))
    return a


def test_continuous_readout_does_not_change_fitted_model():
    a = sample()
    options = dict(n_components=4, n_neighbors=4, epochs=2, random_state=7, device='cpu')
    legacy = run_sigma(a, **options)
    continuous = run_sigma(a, distance_mode='continuous', **options)
    for name in ['sigma_region_probability', 'sigma_inside', 'sigma_boundary']:
        np.testing.assert_array_equal(legacy.obs[name], continuous.obs[name])
    np.testing.assert_array_equal(legacy.obsm['X_sigma_corrected'], continuous.obsm['X_sigma_corrected'])
    assert legacy.uns['sigma']['distance_mode'] == 'legacy'
    assert continuous.uns['sigma_coordinate']['mode'] == 'continuous'
    assert np.sum(legacy.obs.sigma_d_signed == 0) > np.sum(continuous.obs.sigma_d_signed == 0)


def test_weak_anchor_continuous_aliases_and_h5ad_roundtrip(tmp_path):
    a = run_sigma_weak_anchor(sample(), n_components=4, n_neighbors=4,
                              epochs=2, random_state=7, distance_mode='continuous')
    np.testing.assert_array_equal(a.obs.d_signed, a.obs.sigma_d_signed)
    np.testing.assert_allclose(a.obs.d_abs, np.abs(a.obs.sigma_d_signed))
    a.write_h5ad(tmp_path/'continuous.h5ad')
    read = ad.read_h5ad(tmp_path/'continuous.h5ad')
    np.testing.assert_array_equal(read.obs.sigma_inside, a.obs.sigma_inside)
    np.testing.assert_allclose(read.uns['sigma_coordinate']['segments'], a.uns['sigma_coordinate']['segments'])


def test_ion_ranking_does_not_test_overlapping_quantile_groups():
    a = sample()
    a.obs['sigma_d_signed'] = np.r_[-np.ones(18), np.ones(18)]
    a.obs['sigma_inside'] = np.r_[np.ones(18, dtype=bool), np.zeros(18, dtype=bool)]
    result = rank_interface_metabolites(a)
    assert result.negative_near_minus_far.isna().all()
    assert result.positive_near_minus_far.isna().all()
    assert result.p_value.isna().all()
    assert not result.eligible_representative.any()


@pytest.mark.parametrize('near,far,spacing', [(5, 5, 1), (8, 5, 1), (5, 8, 0)])
def test_fixed_window_comparison_rejects_overlap(near, far, spacing):
    with pytest.raises(ValueError, match='near_width < far_start'):
        stroma_near_far_statistics(sample(), {}, near_width=near, far_start=far, spot_spacing=spacing)


def test_stromal_fixed_window_keeps_explicit_zero_side():
    a = ad.AnnData(np.ones((20, 2)))
    a.obs['sigma_d_signed'] = np.r_[np.zeros(10), np.ones(10)*10]
    a.obs['sigma_inside'] = np.r_[np.ones(5, dtype=bool), np.zeros(15, dtype=bool)]
    scores = {0: np.r_[np.ones(5)*100, np.ones(5)*2, np.zeros(10)]}
    table, leading = stroma_near_far_statistics(a, scores, spot_spacing=1)
    assert table.n_near.item() == 5
    assert table.near_minus_far.item() == 2
    assert leading == 0
