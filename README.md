# SIGMA (`sigma-omics`)

SIGMA is a spectral interface graph model for quantitative metabolic field
analysis of tumor interfaces. It constructs a spatial Gaussian graph from
tissue coordinates, transforms tumor/non-tumor annotations into an
interface-aware prior field, and decomposes spatial metabolomics
representations into low-frequency tissue background and high-frequency
boundary-associated residual components. A spectral residual graph
convolutional network then learns boundary-enriched metabolic fields under weak
pathological supervision. Matched spatial transcriptomic profiles can be used
for auxiliary alignment and biological interpretation but are not required for
the separate SM-only weak-anchor workflow.

SIGMA returns a continuous tumor-associated field, an inferred boundary, and a
signed-distance coordinate for ranking and quantifying interface-associated
metabolic features. The scientific definitions and default parameters are
preserved from the HBC515 and HCC reference implementations.

## Install

After the first PyPI release:

```bash
pip install sigma-omics
```

The distribution name is `sigma-omics`; the Python import remains
`sigma_spatial`.

## Install from source

```bash
pip install -e .
```

## Quick start

```python
import scanpy as sc
from sigma_spatial import run_sigma

adata = sc.read_h5ad("sample.h5ad")
result = run_sigma(
    adata,
    anchor_key="sigma_anchor",        # 1=tumor, 0=non-tumor, NaN=unknown
    representation_key="X_harmony",  # validated auxiliary representation
    random_state=0,
)

# Main outputs
result.obs[["sigma_region_probability", "sigma_boundary", "sigma_d_signed"]]
```

For SM-only data, use the separately validated HCC-derived entry point:

```python
from sigma_spatial import run_sigma_weak_anchor

result = run_sigma_weak_anchor(
    adata,
    anchor_key="sigma_anchor",
    already_log=True,
    random_state=0,
)
```

This path preserves the executed HCC P1/P4 loss without introducing an RNA
target. See `docs/input_format.md` and `docs/output_schema.md`.

## Scope

The package should contain reusable SIGMA computation only. Manuscript-specific simulation, benchmarking, GO enrichment, plotting, and sample-specific analyses should live under `examples/` or a separate analysis repository.
