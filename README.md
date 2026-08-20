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

For an unfamiliar input, inspect the declared modality before choosing a
workflow:

```python
from sigma_spatial import inspect_input

assessment = inspect_input(adata)  # uses var['feature_type']: SM, ST, or both
assessment.print_report()
```

SM-only data are routed to the explicitly anchored weak-anchor workflow;
joint SM+ST data use SM as the target and an ST representation for auxiliary
alignment. ST-only data may be used for interface/signature validation, but
SIGMA does not label ST-only outputs as metabolic programs.

## Reproduce a reference dataset

Reference data releases contain the core AnnData, frozen feature ranking, and
frozen program assignment table. These are all required because re-clustering
a ranking is not guaranteed to preserve manuscript cluster identities.

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

The preset fixes the manuscript matrix source, program, distance windows, and
random seed. A frozen representative table can additionally be supplied with
`representatives_path=` when a manuscript panel used a curated final display
set.

## Generate the downstream result set

For a new user, the recommended entry point is the combined workflow:

```python
import sigma_spatial as sigma

sigma.assess_input(adata).print_report()
result = sigma.run_analysis(
    adata,
    output_dir="results/sample_01",
    prefix="sample_01",
    workflow="auto",
    evidence={"same_section_pathology": True},
    report_level="standard",
    random_state=0,
)
print(result.workflow)
print(result.downstream.leading_program)
```

For scientific safety, `workflow="auto"` reads existing provenance or requires
explicit evidence. It does not infer the annotation source from whichever
workflow gives the strongest result. Supported evidence flags are
`same_section_pathology`, `matched_st`, `adjacent_section_pathology`, and
`region_defined`.

`report_level` controls plotting without changing numerical analysis:

- `none`: tables and provenance only;
- `standard`: program patterns, profiles, lollipop, and representative metabolites;
- `complete`: standard report plus influence-range and anisotropy panels;
- `manuscript`: publication-layout report with the complete diagnostics.

After SIGMA has produced a signed-distance field, one call discovers the
metabolic programs and writes the standard figures and tables:

```python
from sigma_spatial import run_downstream_analysis

report = run_downstream_analysis(
    result,
    "results/sample_01",
    prefix="sample_01",
    selection_mode="lambda_profile",
    workflow="direct_pathology",   # one of the four evidence workflows
    signature_scores=None,           # optional matched-ST score arrays
)
print(report.leading_program)
print(report.representatives[["mz", "cluster"]])
```

The report includes spatial program maps, signed-distance profiles,
two-sided near-versus-far interface statistics, five representative metabolite
maps, and directional anisotropy. When independent ST signature scores are
provided, it additionally creates the program-signature heatmap and ST
interface-enrichment plot. The downstream call does not refit the SIGMA model.

`workflow` is selected once from `direct_pathology`, `multiomics_inferred`,
`transferred_pathology`, or `region_defined_disease`. Each workflow has one
fixed lambda-profile preset shared by its samples; users do not tune thresholds
separately for every section. For sparse signed-distance grids, SIGMA can
deterministically reduce the number of bins and records both the requested and
effective bin counts. The exact preset is saved beside every result.

For exact reproduction of a frozen manuscript analysis, set
`selection_mode="reference"` and pass its frozen final ranking table. If the
manuscript program was selected using dataset-specific biological criteria,
also pass `leading_program=<frozen cluster>`; otherwise the leading program is
selected automatically from the two-sided near-versus-far analysis. The JSON
summary records both the automatic result and any explicit reference override.

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

The package contains reusable SIGMA computation and standardized downstream
reporting. Manuscript-specific simulation, benchmarking, GO enrichment, and
sample-specific interpretation remain under `examples/` or the analysis
repository.
