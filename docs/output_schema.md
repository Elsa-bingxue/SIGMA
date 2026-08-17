# Output schema

`run_sigma` returns a new `AnnData` object by default and writes:

| Location | Key | Meaning |
|---|---|---|
| `obs` | `sigma_gaussian_anchor` | spatially diffused weak-anchor field |
| `obs` | `sigma_region_probability_raw` | unnormalized learned probability |
| `obs` | `sigma_region_probability` | reference min-max normalized probability |
| `obs` | `sigma_inside` | probability threshold at the preserved level 0.5 |
| `obs` | `sigma_boundary` | kNN boundary mask |
| `obs` | `sigma_d_signed` | negative inside/tumor, positive outside |
| `obsm` | `X_sigma_msi` | MSI embedding |
| `obsm` | `X_sigma_gauss` | graph-smoothed embedding |
| `obsm` | `X_sigma_residual` | learned spectral residual |
| `obsm` | `X_sigma_corrected` | smoothed embedding plus residual |

The estimator configuration and public input keys are recorded in
`adata.uns["sigma"]`.

The SM-only `run_sigma_weak_anchor` path also writes the legacy HCC keys
`X_msi_embed`, `E_gauss`, `R_spectral`, `E_hat_gauss_spectral`, `f_gauss`,
`p_tumor_gauss_spectral`, `is_boundary`, and `d_signed`. These aliases permit
direct comparison with the frozen P1/P4 notebooks; the `sigma_*` keys remain
the public output schema.
