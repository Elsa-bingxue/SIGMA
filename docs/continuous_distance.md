# Continuous interface coordinates

Version 0.3.2 adds a coordinate readout to original SIGMA; the trainable model and
loss are unchanged. The default remains `distance_mode="legacy"`.

## Geometry and units

`continuous_signed_distance(xy, field)` returns `(distance, inside, segments)`.
The interface is the 0.5 contour of the field, linearly interpolated over a spatial
triangulation. Triangles with an edge longer than `max_edge_factor` (default 1.8)
times the median nearest-neighbour spacing are excluded. Distances are exact
nearest distances to the remaining segments, in the units of `xy`.

The negative side has `field >= 0.5`. The positive side has `field < 0.5`.
Use the separate boolean side mask for exact zeros. Duplicate or collinear spatial
points, absent contours and contours removed entirely by gap masking raise errors.
Coordinate meaning depends on the supplied field; the function does not validate
whether that field is an independently established tumour boundary.

## Fitting

Both `run_sigma(..., distance_mode="continuous")` and
`run_sigma_weak_anchor(..., distance_mode="continuous")` support the new readout.
Continuous distance requires `boundary_mode="legacy"`, which preserves the
unmodified threshold mask. Morphological mask cleanup would define a different
interface and is rejected with this readout.

## Saved fields

```python
from sigma_spatial import update_interface_coordinate
updated = update_interface_coordinate(adata, copy=True)
```

The function requires `obsm['spatial']` and `obs['sigma_region_probability']` by
default. If a stored `sigma_inside` disagrees with the field threshold, it raises
an error rather than silently changing side provenance.

- Existing `sigma_d_signed`, `sigma_inside`, `d_signed` and `d_abs` columns receive
  a `_before_continuous` backup if no backup exists already.
- The canonical distance and side columns are updated; existing aliases are synchronized.
- `uns['sigma_coordinate']` records segments, threshold, gap setting and units.
- `model_refitted=False` and `downstream_requires_recomputation=True` describe the update.
- `sigma_boundary` remains a neighbour-based display mask. Use the stored segments
  to draw the continuous contour.

The first backups are preserved on repeated calls. Downstream results are not
silently modified: recompute selection, distance profiles and clustering explicitly.

## Near/far statistics

Supply `inside_mask` to `side_near_far_statistics` and
`side_near_far_signature_statistics`. Nonzero distance signs must agree with this
mask; exact zeros are allocated using the mask. Without a mask, zeros are excluded
with a warning. High-level downstream analysis forwards `obs['sigma_inside']` when
available.

Programme comparisons with tied quantiles that cause overlap return unavailable
(`NaN`) effects, p-values and FDR, and expose `n_overlap`. Ion-ranking comparisons
also return unavailable effects and p-values for overlap. Fixed-width stromal
comparisons require `0 <= near_width < far_start` and positive spot spacing.

These tests describe within-side contrasts; spatial dependence and independent
biological replication still require appropriate validation designs.
