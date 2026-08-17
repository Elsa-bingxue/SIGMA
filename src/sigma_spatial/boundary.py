import numpy as np
from sklearn.neighbors import NearestNeighbors


def build_boundary_from_binary_mask(xy, inside_mask, k_nn=10):
    """Return the reference notebook's same-and-opposite-neighbour boundary."""
    xy = np.asarray(xy, dtype=float)
    inside = np.asarray(inside_mask, dtype=bool)
    k_nn = min(int(k_nn), xy.shape[0] - 1)
    nn = NearestNeighbors(n_neighbors=k_nn + 1).fit(xy)
    idx = nn.kneighbors(xy, return_distance=False)[:, 1:]
    same = np.any(inside[idx] == inside[:, None], axis=1)
    opposite = np.any(inside[idx] != inside[:, None], axis=1)
    return same & opposite


def build_boundary_from_field(xy, field, level=0.5, k_nn=10):
    """Threshold a continuous field and construct the reference boundary."""
    field = np.asarray(field, dtype=float)
    if not np.all(np.isfinite(field)):
        raise ValueError("The boundary-defining field contains non-finite values.")
    inside = field >= float(level)
    boundary = build_boundary_from_binary_mask(xy, inside, k_nn=k_nn)
    return boundary, inside


def signed_distance_from_boundary_points(xy, boundary_mask, inside_mask):
    """Distance to nearest boundary point; tumor/inside side negative, outside positive."""
    xy = np.asarray(xy, dtype=float)
    b = np.asarray(boundary_mask, dtype=bool)
    inside = np.asarray(inside_mask, dtype=bool)
    if not np.any(b):
        raise ValueError("No boundary points were detected.")
    nn = NearestNeighbors(n_neighbors=1).fit(xy[b])
    dist = nn.kneighbors(xy, return_distance=True)[0][:, 0]
    return np.where(inside, -dist, dist)
