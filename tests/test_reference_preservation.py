import hashlib
from pathlib import Path

import numpy as np

from sigma_spatial.boundary import build_boundary_from_binary_mask
from sigma_spatial.graph import gaussian_knn_graph, gaussian_label_smoothing, smooth_embedding


ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "HBC_515" / "HBC_515_SIGMA.ipynb"
ORIGINAL_SHA256 = "58cdc1c0d18eddcbeab4718fe9fc5748ae852b3f5ff7dc485d49c2400c4f45d4"


def test_original_hbc515_notebook_is_unchanged():
    assert hashlib.sha256(ORIGINAL.read_bytes()).hexdigest() == ORIGINAL_SHA256


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
