import hashlib
from pathlib import Path

import numpy as np

from sigma_spatial.boundary import build_boundary_from_binary_mask, build_boundary_from_field
from sigma_spatial.graph import gaussian_knn_graph, gaussian_label_smoothing, smooth_embedding


ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "HBC_515" / "HBC_515_SIGMA.ipynb"
ORIGINAL_SHA256 = "58cdc1c0d18eddcbeab4718fe9fc5748ae852b3f5ff7dc485d49c2400c4f45d4"
ORIGINAL_HBC525 = ROOT / "HBC_525" / "HBC-525_SIGMA.ipynb"
ORIGINAL_HBC525_SHA256 = "d085a2d3a236b99095a3b18c4b5763adeb5a4cab127a08e22533c31cd1a3c52e"
ORIGINAL_HLC091 = ROOT / "HLC_091" / "HLC_091_SIGMA.ipynb"
ORIGINAL_HLC091_SHA256 = "f7336555559e3e58195440a4993ffe3c12d961ce5723953ed64171a108a85eec"
ORIGINAL_HLC276 = ROOT / "HLC_276" / "HLC_276_SIGMA.ipynb"
ORIGINAL_HLC276_SHA256 = "36186865a1d8a2f7e0226c31c8764fb944d7c19e3e32ea59e62fcf2097b129ea"
ORIGINAL_CCRCC = ROOT / "ccRCC_Y27T" / "ccRCC_SIGMA.ipynb"
ORIGINAL_CCRCC_SHA256 = "b2df2184f159b5e8e2d183bce995f9c62992c1f0fbfdd365bf7838382885abe0"
ORIGINAL_GBM = ROOT / "GBM" / "GBM_SIGMA.ipynb"
ORIGINAL_GBM_SHA256 = "33a5ec01c4c028730ca3c025d115ccd65313d696966764e363e62268a507787a"
ADDITIONAL_NOTEBOOKS = {
    ROOT / "HCC_P1/HCC_spectral GCN202606012.ipynb": "ea728cdb9b19924178c57b11a0eaa10666ebe0479875b8850086b6f142147db9",
    ROOT / "HCC_P4/HCC_spectral GCN202606012.ipynb": "b51ed71e63fc246402ab6e66c3a715cf413cabba111fc96f22dbfd1ef23a8c1c",
    ROOT / "HPD_A1/HPDA1_spectral GCN20260626.ipynb": "b4309b936e552d289af881d8a659afe120f5747222635111878148c0c3b1e177",
    ROOT / "HPD_B1/HPDB1_spectral GCN20260626.ipynb": "6a50a434eb00c7222e1d316bf33b2dc82b13a68631d39027fdb107d823e502a0",
    ROOT / "HPD_C1/HPDC1_spectral GCN20260626.ipynb": "43f1ab0f3347f26211e921bfd367f505a2bbaf93d44d57c1beaf413e0df6d1ad",
}


def test_original_hbc515_notebook_is_unchanged():
    assert hashlib.sha256(ORIGINAL.read_bytes()).hexdigest() == ORIGINAL_SHA256


def test_original_hbc525_notebook_is_unchanged():
    assert hashlib.sha256(ORIGINAL_HBC525.read_bytes()).hexdigest() == ORIGINAL_HBC525_SHA256


def test_original_hlc_notebooks_are_unchanged():
    assert hashlib.sha256(ORIGINAL_HLC091.read_bytes()).hexdigest() == ORIGINAL_HLC091_SHA256
    assert hashlib.sha256(ORIGINAL_HLC276.read_bytes()).hexdigest() == ORIGINAL_HLC276_SHA256


def test_original_joint_omics_notebooks_are_unchanged():
    assert hashlib.sha256(ORIGINAL_CCRCC.read_bytes()).hexdigest() == ORIGINAL_CCRCC_SHA256
    assert hashlib.sha256(ORIGINAL_GBM.read_bytes()).hexdigest() == ORIGINAL_GBM_SHA256


def test_additional_dataset_notebooks_are_unchanged():
    for path, expected in ADDITIONAL_NOTEBOOKS.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected


def test_reference_directed_gaussian_graph_values():
    xy = np.array([[0., 0.], [1., 0.], [0., 2.], [3., 0.]], dtype=np.float32)
    edge_index, weights, sigma = gaussian_knn_graph(xy, k=2, symmetrize=False)
    expected_edges = np.array([[0, 0, 1, 1, 2, 2, 3, 3],
                               [1, 2, 0, 3, 0, 1, 1, 0]])
    expected_dist = np.array([1., 2., 1., 2., 2., np.sqrt(5), 2., 3.])
    expected_sigma = np.median(expected_dist.reshape(4, 2))
    np.testing.assert_array_equal(edge_index, expected_edges)
    np.testing.assert_allclose(sigma, expected_sigma)
    np.testing.assert_allclose(weights, np.exp(-(expected_dist ** 2) /
                                                (2 * expected_sigma ** 2)), rtol=1e-6)


def test_label_and_embedding_diffusion_match_frozen_values():
    edges = np.array([[0, 0, 1, 1, 2, 2], [1, 2, 0, 2, 0, 1]])
    weights = np.ones(6, dtype=np.float32)
    y = np.array([1, -1, 0])
    mask = y >= 0
    field = gaussian_label_smoothing(edges, weights, y, mask, n_iter=2, alpha=0.5)
    np.testing.assert_allclose(field, [1.0, 0.25, 0.0], rtol=0, atol=1e-7)
    embedding = np.array([[1., 0.], [0., 1.], [0., 0.]], dtype=np.float32)
    smoothed = smooth_embedding(edges, weights, embedding, n_iter=1, alpha=0.5)
    np.testing.assert_allclose(smoothed, [[0.5, 0.25], [0.25, 0.5], [0.25, 0.25]])


def test_boundary_requires_same_and_opposite_neighbours():
    xy = np.array([[0., 0.], [1., 0.], [2., 0.], [3., 0.]])
    inside = np.array([True, True, False, False])
    np.testing.assert_array_equal(
        build_boundary_from_binary_mask(xy, inside, k_nn=2),
        np.array([True, True, True, True]),
    )


def test_robust_boundary_removes_tiny_island_without_changing_default():
    xx, yy = np.meshgrid(np.arange(8), np.arange(8))
    xy = np.column_stack([xx.ravel(), yy.ravel()])
    field = np.zeros(64)
    main = (xy[:, 0] >= 2) & (xy[:, 0] <= 5) & (xy[:, 1] >= 2) & (xy[:, 1] <= 5)
    field[main] = 1
    field[(xy[:, 0] == 7) & (xy[:, 1] == 7)] = 1
    default_boundary, default_inside = build_boundary_from_field(xy, field, k_nn=4)
    legacy_boundary, legacy_inside = build_boundary_from_field(
        xy, field, k_nn=4, boundary_mode="legacy"
    )
    robust_boundary, robust_inside = build_boundary_from_field(
        xy, field, k_nn=4, boundary_mode="robust",
        min_component_size=2, max_hole_size=0,
    )
    np.testing.assert_array_equal(default_inside, legacy_inside)
    np.testing.assert_array_equal(default_boundary, legacy_boundary)
    assert default_inside.sum() == robust_inside.sum() + 1
    assert robust_boundary.sum() <= default_boundary.sum()
