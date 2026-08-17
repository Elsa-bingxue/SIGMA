"""SM-only weak-anchor SIGMA, preserving the HCC P1/P4 workflow."""

import scanpy as sc

from sigma_spatial import run_sigma_weak_anchor
from sigma_spatial.plotting import plot_boundary, plot_signed_distance


adata = sc.read_h5ad("sm_only_sample.h5ad")
result = run_sigma_weak_anchor(
    adata,
    anchor_key="sigma_anchor",  # 1=tumor, 0=non-tumor, NaN=unknown
    already_log=True,
    random_state=0,
)
result.write_h5ad("sm_only_sample_sigma.h5ad")

plot_boundary(result)
plot_signed_distance(result)
