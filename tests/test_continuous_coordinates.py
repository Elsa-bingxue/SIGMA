import numpy as np
import pytest
import anndata as ad

from sigma_spatial.coordinates import (
    continuous_signed_distance, interface_side_masks, update_interface_coordinate,
    distances_to_segments,
)
from sigma_spatial.boundary import build_boundary_from_field, signed_distance_from_boundary_points
from sigma_spatial.selection import side_near_far_statistics


def grid():
    x, y = np.meshgrid(np.arange(-3., 4), np.arange(-3., 4))
    xy = np.c_[x.ravel(), y.ravel()]
    p = (xy[:, 0] + 4) / 8
    return xy, p


def test_continuous_distance_matches_analytic_line_and_rescales():
    xy, p = grid()
    d, inside, segments = continuous_signed_distance(xy, p, level=.55)
    np.testing.assert_allclose(d, .4 - xy[:, 0], atol=1e-12)
    np.testing.assert_array_equal(inside, p >= .55)
    d2, _, _ = continuous_signed_distance(xy * 21, p, level=.55)
    np.testing.assert_allclose(d2, d * 21, atol=1e-12)
    b, side = build_boundary_from_field(xy, p, level=.55)
    legacy = signed_distance_from_boundary_points(xy, b, side)
    assert (legacy == 0).sum() > (d == 0).sum()


def test_side_labels_include_both_signed_zeros_without_epsilon_hack():
    d = np.array([-2., -0., 0., 2.])
    masks = interface_side_masks(d, [1, 1, 0, 0])
    np.testing.assert_array_equal(masks['negative'], [True, True, False, False])
    np.testing.assert_array_equal(masks['positive'], [False, False, True, True])
    with pytest.warns(UserWarning, match='no side label'):
        masks = interface_side_masks(d)
    assert sum(m.sum() for m in masks.values()) == 2
    with pytest.raises(ValueError, match='disagree'):
        interface_side_masks(d, [0, 1, 0, 0])


def test_exact_contour_zero_has_explicit_side():
    xy, p = grid()
    d, inside, _ = continuous_signed_distance(xy, p)
    masks = interface_side_masks(d, inside)
    assert (d == 0).sum() == 7
    assert np.all(masks['negative'][d == 0])


def test_no_contour_or_gap_bridge_is_rejected():
    xy, p = grid()
    with pytest.raises(ValueError, match='No separating contour'):
        continuous_signed_distance(xy, np.ones(len(xy)))
    keep = abs(xy[:, 0]) >= 2
    with pytest.raises(ValueError, match='gap masking'):
        continuous_signed_distance(xy[keep], p[keep])


def test_segment_distance_projects_to_interior_and_clamps_endpoints():
    d = distances_to_segments([[1, 3], [-1, 0], [3, 1]], [[[0, 0], [2, 0]]])
    np.testing.assert_allclose(d, [3, 1, np.sqrt(2)])


def test_frozen_upgrade_preserves_input_and_legacy_values():
    xy, p = grid()
    obj = ad.AnnData(np.ones((len(xy), 1)))
    obj.obsm['spatial'] = xy
    obj.obs['sigma_region_probability'] = p
    b, inside = build_boundary_from_field(xy, p)
    old = signed_distance_from_boundary_points(xy, b, inside)
    obj.obs['sigma_d_signed'] = old
    obj.obs['sigma_inside'] = inside
    upgraded = update_interface_coordinate(obj)
    np.testing.assert_array_equal(obj.obs.sigma_d_signed, old)
    np.testing.assert_array_equal(upgraded.obs.sigma_d_signed_before_continuous, old)
    np.testing.assert_allclose(upgraded.obs.sigma_d_signed, -xy[:, 0])
    assert upgraded.uns['sigma_coordinate']['downstream_requires_recomputation']


def test_near_far_retains_zero_spots_on_each_correct_side():
    d = np.array([-4., -3., -2., -1., -0., 0., 1., 2., 3., 4.])
    y = np.array([0., 0., 0., 1., 10., 20., 1., 0., 0., 0.])
    table = side_near_far_statistics({'p': y}, d, inside_mask=np.arange(10) < 5,
        near_quantile=.25, far_quantile=.75, min_group_size=2).set_index('orientation')
    assert table.loc['negative', 'near_minus_far'] == 5.5
    assert table.loc['positive', 'near_minus_far'] == 10.5


def test_frozen_upgrade_synchronizes_existing_weak_anchor_aliases():
    xy, p = grid()
    obj = ad.AnnData(np.ones((len(xy), 1)))
    obj.obsm['spatial'] = xy
    obj.obs['sigma_region_probability'] = p
    b, inside = build_boundary_from_field(xy, p)
    old = signed_distance_from_boundary_points(xy, b, inside)
    obj.obs['sigma_d_signed'] = old
    obj.obs['sigma_inside'] = inside
    obj.obs['d_signed'] = old
    obj.obs['d_abs'] = np.abs(old)
    upgraded = update_interface_coordinate(obj)
    np.testing.assert_allclose(upgraded.obs.d_signed, -xy[:, 0])
    np.testing.assert_allclose(upgraded.obs.d_abs, np.abs(xy[:, 0]))
    np.testing.assert_array_equal(upgraded.obs.d_signed_before_continuous, old)
    np.testing.assert_array_equal(upgraded.obs.d_abs_before_continuous, np.abs(old))
    np.testing.assert_array_equal(obj.obs.d_signed, old)


def test_tied_near_far_groups_do_not_test_the_same_observations_twice():
    d = np.r_[-np.ones(20), np.ones(20)]
    y = np.arange(40, dtype=float)
    table = side_near_far_statistics({'p': y}, d, min_group_size=2)
    assert table.n_overlap.eq(20).all()
    assert table.near_minus_far.isna().all()
    assert table.p_value.isna().all()
    assert table.fdr.isna().all()
