# SIGMA (`sigma-omics`)

SIGMA is a Python toolkit for spatial interface analysis integrating mass-spectrometry imaging (MSI) and spatial transcriptomic representations.

The scientific definitions and default parameters are preserved from the
HBC515 reference implementation.

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
