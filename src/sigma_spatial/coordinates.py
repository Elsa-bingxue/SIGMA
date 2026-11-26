"""Continuous interface coordinates, separate from discrete side membership.

Distances are in the units of ``xy``. Zero denotes the interpolated contour,
not the band of spots with neighbours on both sides. Legacy point-band
distances remain available for reproducing the original notebooks.
"""
from __future__ import annotations

import warnings
import numpy as np
from sklearn.neighbors import NearestNeighbors


def interface_side_masks(signed_distance, inside_mask=None):
    """Return negative/positive side masks without assigning zero to either side.

    An explicit 0/1 inside mask includes zero-distance observations correctly.
    Without it, zero is ambiguous and is excluded with a warning. Signed zero
    is not used as metadata because file formats need not preserve its sign.
    """
    d = np.asarray(signed_distance, dtype=float)
    if d.ndim != 1:
        raise ValueError("signed_distance must be one-dimensional")
    finite = np.isfinite(d)
    if inside_mask is None:
        if np.any(finite & (d == 0)):
            warnings.warn("Zero-distance observations have no side label; supply inside_mask. "
                          "They are excluded from side-specific comparisons.", UserWarning,
                          stacklevel=2)
        return {"negative": finite & (d < 0), "positive": finite & (d > 0)}
    raw = np.asarray(inside_mask)
    if raw.shape != d.shape or not np.isin(raw, [0, 1]).all():
        raise ValueError("inside_mask must be a complete boolean/0-1 vector matching distance")
    inside = raw.astype(bool)
    if np.any(finite & (((d < 0) & ~inside) | ((d > 0) & inside))):
        raise ValueError("Nonzero signed distances disagree with inside_mask")
    return {"negative": finite & inside, "positive": finite & ~inside}


def contour_segments(xy, field, *, level=0.5, max_edge_factor=1.8):
    """Linearly interpolate a contour on a gap-masked spatial triangulation.

    Triangles with an edge longer than max_edge_factor times the median nearest
    neighbour spacing are excluded. No segments are drawn across those gaps.
    A missing contour raises an error rather than producing a constant distance.
    """
    from matplotlib.tri import Triangulation

    xy = np.asarray(xy, dtype=float)
    value = np.asarray(field, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or value.shape != (len(xy),):
        raise ValueError("Require xy shape (n, 2) and field shape (n,)")
    if len(xy) < 3 or not np.isfinite(xy).all() or not np.isfinite(value).all():
        raise ValueError("At least three finite spatial points and field values are required")
    if len(np.unique(xy, axis=0)) != len(xy):
        raise ValueError("Duplicate spatial positions must be resolved before contouring")
    if not np.isfinite(level) or not np.isfinite(max_edge_factor) or max_edge_factor <= 0:
        raise ValueError("level must be finite and max_edge_factor must be positive")
    if not (np.min(value) < level < np.max(value)):
        raise ValueError("No separating contour: level must lie strictly within the field range")
    spacing = np.median(NearestNeighbors(n_neighbors=2).fit(xy).kneighbors(xy)[0][:, 1])
    try:
        triangles = Triangulation(xy[:, 0], xy[:, 1]).triangles
    except (ValueError, RuntimeError) as exc:
        raise ValueError("Contour construction requires non-collinear spatial positions") from exc
    points = xy[triangles]
    keep = np.linalg.norm(points - np.roll(points, 1, axis=1), axis=2).max(1) <= max_edge_factor * spacing
    points, v = points[keep], value[triangles[keep]] - level
    segments = []
    for p, z in zip(points, v):
        hits = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            if (z[a] >= 0) != (z[b] >= 0):
                hits.append(p[a] + (-z[a] / (z[b] - z[a])) * (p[b] - p[a]))
        if len(hits) == 2 and np.linalg.norm(hits[1] - hits[0]) > spacing * 1e-12:
            segments.append(hits)
    if not segments:
        raise ValueError("No contour remains within sampled tissue after gap masking")
    return np.asarray(segments, dtype=float)


def distances_to_segments(xy, segments, *, chunk_size=128):
    """Exact nearest point-to-segment distances; no vertex sampling approximation."""
    xy = np.asarray(xy, dtype=float)
    seg = np.asarray(segments, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or seg.ndim != 3 or seg.shape[1:] != (2, 2):
        raise ValueError("Require xy (n, 2) and segments (m, 2, 2)")
    if not len(seg) or not np.isfinite(xy).all() or not np.isfinite(seg).all():
        raise ValueError("Finite points and at least one finite segment are required")
    if int(chunk_size) < 1:
        raise ValueError("chunk_size must be positive")
    start, delta = seg[:, 0], seg[:, 1] - seg[:, 0]
    length2 = np.sum(delta ** 2, axis=1)
    if np.any(length2 <= 0):
        raise ValueError("Zero-length contour segments are invalid")
    out = np.empty(len(xy))
    for lo in range(0, len(xy), int(chunk_size)):
        rel = xy[lo:lo + int(chunk_size), None] - start
        t = np.clip(np.sum(rel * delta, axis=2) / length2, 0, 1)
        out[lo:lo + int(chunk_size)] = np.sqrt(np.min(np.sum((rel - t[..., None] * delta) ** 2, axis=2), axis=1))
    return out


def continuous_signed_distance(xy, field, *, level=0.5, max_edge_factor=1.8):
    """Return (distance, inside, segments), with field >= level on the negative side."""
    seg = contour_segments(xy, field, level=level, max_edge_factor=max_edge_factor)
    inside = np.asarray(field, dtype=float) >= level
    distance = distances_to_segments(xy, seg) * np.where(inside, -1.0, 1.0)
    return distance, inside, seg


def distance_readout(xy, field, boundary, inside, *, mode="continuous",
                     level=0.5, max_edge_factor=1.8):
    """Shared distance adapter for fitted models; morphology must match the field."""
    from .boundary import signed_distance_from_boundary_points
    if mode == "legacy":
        return signed_distance_from_boundary_points(xy, boundary, inside), None
    if mode != "continuous":
        raise ValueError("distance_mode must be 'continuous' or 'legacy'")
    if not np.array_equal(np.asarray(inside, bool), np.asarray(field) >= level):
        raise ValueError("Continuous distance requires the unmodified threshold mask; "
                         "morphological boundary cleanup defines a different interface")
    d, _, seg = continuous_signed_distance(xy, field, level=level, max_edge_factor=max_edge_factor)
    return d, seg


def update_interface_coordinate(adata, *, copy=True, field_key="sigma_region_probability",
                                spatial_key="spatial", level=0.5, max_edge_factor=1.8):
    """Upgrade a frozen probability field without training or changing its boundary.

    Existing distance/side columns are retained with a ``_before_continuous``
    suffix. Existing weak-anchor aliases ``d_signed`` and ``d_abs`` are kept
    consistent with the upgraded coordinate. Downstream candidates and patterns must be recomputed explicitly.
    The neighbour-based ``sigma_boundary`` remains a legacy display mask; use
    ``uns['sigma_coordinate']['segments']`` for the actual continuous interface.
    """
    d, inside, seg = continuous_signed_distance(adata.obsm[spatial_key],
        adata.obs[field_key].to_numpy(float), level=level, max_edge_factor=max_edge_factor)
    if "sigma_inside" in adata.obs and not np.array_equal(adata.obs.sigma_inside.to_numpy(bool), inside):
        raise ValueError("Stored side mask differs from the probability threshold; resolve provenance first")
    out = adata.copy() if copy else adata
    for key in ("sigma_d_signed", "sigma_inside", "d_signed", "d_abs"):
        backup = key + "_before_continuous"
        if key in out.obs and backup not in out.obs:
            out.obs[backup] = out.obs[key].copy()
    out.obs["sigma_d_signed"] = d
    out.obs["sigma_inside"] = inside
    if "d_signed" in out.obs:
        out.obs["d_signed"] = d
    if "d_abs" in out.obs:
        out.obs["d_abs"] = np.abs(d)
    out.uns["sigma_coordinate"] = dict(mode="continuous", level=float(level),
        max_edge_factor=float(max_edge_factor), field_key=field_key, spatial_key=spatial_key,
        units="spatial coordinate units", segments=seg, model_refitted=False,
        downstream_requires_recomputation=True)
    out.uns.setdefault("sigma", {})["distance_mode"] = "continuous"
    return out
