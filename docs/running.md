# Running SIGMA 0.3.2

## 1. Check the input

The paired path requires non-negative MSI intensities, 2D coordinates, both tumour
and outside anchors, and matched RNA features in the same spot order. Registration
and RNA preprocessing are performed upstream. See [input format](input_format.md).
`run_sigma` does not perform registration and does not infer tissue anchors from RNA.

```python
import scanpy as sc
from sigma_spatial import run_sigma, run_downstream_analysis
adata = sc.read_h5ad("sample.h5ad")
```

## 2. Fit the interface

```python
result = run_sigma(
    adata,
    anchor_key="sigma_anchor",
    representation_key="X_harmony",
    distance_mode="continuous",
    random_state=0,
)
result.write_h5ad("sample_sigma.h5ad")
```

The paired default is 1,000 epochs. The generated demo reduces training for a quick
installation check. Continuous distance uses the 0.5 field contour and requires the
unmodified threshold mask (`boundary_mode="legacy"`). If the field has no contour,
the continuous-distance routine raises an error rather than inventing an interface.

Without matched RNA, call `run_sigma_weak_anchor` with explicit anchors instead.
Set `already_log=False` for non-negative unlogged MSI, or `True` only when the input
was already log-transformed. This path omits the RNA loss.

## 3. Discover programmes and compare near/far positions

```python
report = run_downstream_analysis(
    result,
    output_dir="results/sample",
    prefix="sample",
    matrix_source="X",
    selection_mode="standardized",
    report_level="standard",
    random_state=0,
)
```

Choose `matrix_source` explicitly for the intended MSI matrix (`X`, `raw`,
`layer:<name>` or `msi_uns`). In this API, `raw` means `layers["raw"]`, not
`AnnData.raw`. The example uses the MSI matrix in `X` throughout. This call keeps the
fitted interface fixed, ranks ions, groups distance profiles and writes tables and
figures. It forwards `sigma_inside` for same-side near/far comparisons.

The default exploratory ranking (`standardized`) and archived lambda-profile
rankings are different analyses. To reproduce a manuscript P0/P3 identity, supply
the correct frozen ranking and assignment table with `run_reference_analysis`.
Do not label new clusters with historical identities merely because their numbers
match. See [workflow and reproduction notes](workflows_and_reproducibility.md).

## 4. Update a saved field

```python
from sigma_spatial import update_interface_coordinate
saved = sc.read_h5ad("sample_sigma_legacy.h5ad")
updated = update_interface_coordinate(saved, copy=True)
updated.write_h5ad("sample_sigma_continuous.h5ad")
```

Run step 3 again on `updated`; old profiles and memberships are not automatically
recalculated. Exact-zero side labels come from `sigma_inside`. Overlapping groups
produce unavailable statistics rather than apparently significant contrasts.

## Command-line and configuration runs

`sigma-omics report` is a legacy-default one-command report. It does not expose a
continuous-distance command-line option in version 0.3.2. Use the Python calls above
when choosing distance mode explicitly.

The GitHub checkout also provides `python scripts/run_dataset.py config.json`.
This is a repository helper using the public `run_analysis` API; `sigma-omics run-config`
is not a command in the published 0.3.2 wheel. Set
`"core_kwargs": {"distance_mode": "continuous"}` in a configuration for a new fit.
Already fitted objects require the explicit migration in step 4.

## Manuscript reproduction

The installed package retains a historical figure manifest for compatibility.
Its numbering is not the current six-figure layout. The manuscript archive's root
README, `Figure_Legends` and `Source_Data/INDEX.md` identify the current figures.
The generated demo neither reproduces those figures nor downloads their raw data.
