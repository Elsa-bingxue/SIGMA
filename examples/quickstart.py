"""Minimal SIGMA workflow for a user-supplied AnnData file."""

import scanpy as sc

from sigma_spatial import run_sigma
from sigma_spatial.plotting import plot_boundary, plot_region_probability


adata = sc.read_h5ad("sample.h5ad")
result = run_sigma(
    adata,
    anchor_key="sigma_anchor",
    representation_key="X_harmony",
    random_state=0,
)
result.write_h5ad("sample_sigma.h5ad")

plot_region_probability(result)
plot_boundary(result)
