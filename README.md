# SIGMA (`sigma-omics`)

SIGMA is a Python toolkit for spatial interface analysis integrating mass-spectrometry imaging (MSI) and spatial transcriptomic representations.

> This package is currently an alpha release. Its scientific definitions are
> preserved from the HBC515 reference implementation while the API and
> multi-dataset validation are being completed.

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

## Minimal usage

```python
import scanpy as sc
from sigma_spatial import SIGMA

adata = sc.read_h5ad("sample.h5ad")
SIGMA(seed=0).fit(
    adata,
    annotation_key="annotation",
    tumor_label="Tumor",
    stroma_label="Stroma",
    rna_key="X_harmony",
)

# Main outputs
adata.obs[["sigma_region_probability", "sigma_boundary", "sigma_d_signed"]]
```

## Scope

The package should contain reusable SIGMA computation only. Manuscript-specific simulation, benchmarking, GO enrichment, plotting, and sample-specific analyses should live under `examples/` or a separate analysis repository.
