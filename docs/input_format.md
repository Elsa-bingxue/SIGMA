# Input format

SIGMA accepts an `AnnData` object with:

- a non-negative spot-by-feature matrix in `adata.X`;
- two-dimensional coordinates in `adata.obsm["spatial"]`;
- weak anchors in `adata.obs["sigma_anchor"]`, encoded as `1` for tumor,
  `0` for non-tumor, and missing for unknown spots;
- a validated auxiliary representation such as matched RNA PCs in
  `adata.obsm["X_harmony"]`.

CSV inputs can be loaded with `load_spatial_metabolomics`. Rows are matched by
spot identifier, not by their existing order.

```python
from sigma_spatial import load_spatial_metabolomics

adata = load_spatial_metabolomics(
    "metabolites.csv",
    "coordinates.csv",
    "anchors.csv",
)
```

## SM-only weak-anchor workflow

For SM-only data with weak tumor/non-tumor anchors, use the separate HCC-derived
workflow. This path does not add or substitute an RNA loss:

```python
from sigma_spatial import run_sigma_weak_anchor

result = run_sigma_weak_anchor(
    adata,
    anchor_key="sigma_anchor",
    already_log=True,
    random_state=0,
)
```

The defaults preserve the executed HCC P1/P4 notebook branch. Set
`already_log=False` only when the supplied intensities have not already been
log-transformed.

## Frozen matrix source for downstream analysis

Downstream feature ranking never guesses whether a matrix is raw or processed.
Pass the source used by the dataset's reference analysis:

```python
run_downstream_analysis(..., matrix_source="X")
run_downstream_analysis(..., matrix_source="raw")
run_downstream_analysis(..., matrix_source="layer:normalized")
run_downstream_analysis(..., matrix_source="msi_uns")
```

For joint ST+SM storage, `feature_type == "SM"` (or `type == "SM"`) is applied
before selecting `X` or a layer. This only selects the metabolomics columns; it
does not imply that ST was used to define the interface. The selected source is
recorded in the downstream summary JSON.
