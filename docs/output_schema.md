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

## Downstream report files

`run_downstream_analysis` writes a provenance-tracked result directory. With
prefix `<sample>`, its core outputs are:

| File | Meaning |
|---|---|
| `<sample>_program_assignments.csv` | metabolite-to-program assignments and SM-local indices |
| `<sample>_program_centroid_profiles.csv` | program profiles along signed distance |
| `<sample>_program_scores.csv` | spot-level program scores |
| `<sample>_interface_program_statistics.csv` | near-versus-far tests on both interface sides |
| `<sample>_selected_interface_programs.csv` | program/side pairs passing the requested criteria |
| `<sample>_representative_metabolites.csv` | representatives from the leading interface program |
| `<sample>_metabolic_patterns_and_distance_profiles.*` | program maps and profiles |
| `<sample>_interface_program_enrichment.*` | two-sided interface-program summary |
| `<sample>_boundary_program_<k>_metabolites.*` | representative metabolite maps |
| `<sample>_anisotropy.*` | sector-wise influence-range result |
| `<sample>_downstream_summary.json` | parameters and output status |

If independent matched-ST scores are supplied, the directory also contains
`<sample>_st_program_associations.csv`, `<sample>_st_program_heatmap.*`,
`<sample>_st_interface_statistics.csv`, and
`<sample>_st_interface_lollipop.*`.
# Unified run provenance

`run_analysis()` writes `<prefix>_run_config.json` and records the same payload
in `adata.uns['sigma_analysis_run']`. The payload includes the resolved evidence
workflow, core route, selection mode, report level, matrix source, package
version, and random seed.

Plot generation is controlled by `report_level` and does not alter tables,
program assignments, representative ranking, or fitted SIGMA fields.
