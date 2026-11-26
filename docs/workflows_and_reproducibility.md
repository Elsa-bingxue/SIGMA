# Workflows and reproducibility

This document contains the workflow and manuscript-reproduction details kept
out of the package landing page.

## Evidence workflows

`run_analysis()` supports four evidence workflows:

- `direct_pathology`: same-section pathological annotation;
- `multiomics_inferred`: annotation supported by matched spatial data;
- `transferred_pathology`: pathology transferred from an adjacent section;
- `region_defined_disease`: named disease regions or compartments.

With `workflow="auto"`, SIGMA reads recorded provenance or the supplied
evidence flags: `same_section_pathology`, `matched_st`,
`adjacent_section_pathology` and `region_defined`.

## Report levels

`report_level` controls output generation without changing the fitted model:

- `none`: tables and provenance;
- `standard`: program maps, distance profiles, enrichment and representative
  metabolites;
- `complete`: standard outputs plus influence-range and anisotropy panels;
- `manuscript`: publication-layout report with complete diagnostics.

## Program selection

The default selection is workflow-aware. Direct, multi-omics-inferred and
transferred-pathology analyses use same-side near-versus-far comparisons.
Region-defined disease analyses use a two-stage boundary-localized rule: an
initial effect-size and FDR gate is followed by ranking on median feature-level
interface support. Matched-ST signatures are supplied as contextual readouts rather than
programme-selection criteria. They are not independent validation if RNA was used
to fit the field or construct anchors.

Program-selection metrics and the resolved strategy are written beside each
report. Sparse signed-distance grids may use fewer bins than requested; both
values are recorded in the summary.

## Reference analyses

Reference releases may include the core `AnnData`, a frozen feature ranking and
a frozen program-assignment table. Use all three when manuscript cluster
identities must remain fixed:

```python
import scanpy as sc
from sigma_spatial import run_reference_analysis

adata = sc.read_h5ad("HBC515_SIGMA_core.h5ad")
report = run_reference_analysis(
    adata,
    dataset="HBC515",
    ranking_path="BC515_lambda_ranking.csv",
    assignments_path="BC515_metabolic_program_assignments.csv",
    output_dir="results/HBC515_reference",
)
```

A frozen representative table can be supplied with `representatives_path=`.
For other frozen analyses, use `selection_mode="reference"` and provide the
final ranking table. An explicit `leading_program=` may be used when the
manuscript program identity is already fixed.

## Current and historical figures

The package's built-in manuscript manifest is retained for historical compatibility.
It does not identify the current six-figure manuscript layout. Use the manuscript
archive README, complete figure legends and source-data index for current panels.
Raw data and frozen fits are not included in the wheel. A new fit is not expected
to preserve historical programme numbers without the frozen assignment table.

## Provenance

`run_analysis()` writes `<prefix>_run_config.json` and stores the same
information in `adata.uns["sigma_analysis_run"]`. The record includes the
resolved workflow, core route, selection mode, report level, matrix source,
package version and random seed.
