# SIGMA (`sigma-omics`)

**SIGMA resolves pathological metabolic transitions with interface-aware
spectral graph learning.**

SIGMA is a Python framework for identifying pathological interfaces from
spatial metabolomics. It combines tissue coordinates, metabolite intensities
and pathology-derived or weak spatial anchors to recover continuous interface
fields and boundary-associated metabolic programs.

Spatial transcriptomics is optional. When matched data are available, it can
support auxiliary alignment and independent biological validation.

## Install

```bash
pip install sigma-omics
```

The distribution is named `sigma-omics`; the Python package is
`sigma_spatial`.

## Inputs

SIGMA uses an `AnnData` object containing:

- a spot-by-metabolite intensity matrix;
- two-dimensional tissue coordinates;
- tumor/non-tumor annotations, transferred pathology labels or weak spatial
  anchors;
- an optional matched transcriptomic representation.

See the [input format guide](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/input_format.md)
for supported keys and loading examples.

## Quick start

```python
import scanpy as sc
from sigma_spatial import run_analysis

adata = sc.read_h5ad("sample.h5ad")

result = run_analysis(
    adata,
    output_dir="results/sample",
    prefix="sample",
    workflow="auto",
    evidence={"same_section_pathology": True},
    report_level="standard",
    random_state=0,
)

print(result.workflow)
print(result.downstream.leading_program)
```

For an unfamiliar dataset, inspect its declared modalities first:

```python
from sigma_spatial import inspect_input

inspect_input(adata).print_report()
```

The same standard report is available from the command line:

```bash
sigma-omics report sample.h5ad results/sample \
  --workflow auto \
  --same-section-pathology \
  --report-level standard \
  --seed 0
```

For repeated or multi-dataset runs, use a portable JSON configuration:

```bash
sigma-omics run-config examples/configs/direct_pathology.json
```

The GitHub repository includes configuration templates and thin notebooks
under `examples/`.

## Outputs

SIGMA returns:

- a continuous pathological-region probability field;
- an inferred interface and signed-distance coordinate;
- spatial metabolic programs and their distance profiles;
- ranked interface-associated metabolites;
- interface influence ranges and directional anisotropy;
- optional matched-ST validation summaries.

Core fields are stored in `AnnData`; downstream tables and publication-quality
figures are written to the selected output directory. See
[output schema](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/output_schema.md)
for the complete list.

## Supported workflows

SIGMA supports same-section pathology, joint SM+ST inference, transferred
pathology and region-defined disease analyses. Spatial metabolomics remains the
target modality in every workflow. An explicitly anchored SM-only workflow is
also available for datasets without matched transcriptomics.

```python
from sigma_spatial import run_sigma_weak_anchor

result = run_sigma_weak_anchor(
    adata,
    anchor_key="sigma_anchor",
    already_log=True,
    random_state=0,
)
```

## Reproducibility

Every run records its workflow, parameters, matrix source, package version and
random seed. Reference analyses can additionally use frozen rankings and
program assignments to reproduce manuscript panels.

Detailed workflow rules, program selection and manuscript reproduction are
documented in the
[workflow and reproducibility guide](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/workflows_and_reproducibility.md).

## Data availability

Reference data are not bundled with the package. The study uses published
breast, lung, ccRCC, GBM, HCC and Parkinson's disease spatial-omics datasets.
Download links, source papers, SIGMA sample names and preprocessing notes are
listed in the
[reference-data guide](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/reference_data.md).

## Links

- [GitHub repository](https://github.com/Elsa-bingxue/SIGMA)
- [Input format](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/input_format.md)
- [Output schema](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/output_schema.md)
- [Weak-anchor validation](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/weak_anchor_validation.md)
- [Result generation](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/result_generation.md)
- [Reference data](https://github.com/Elsa-bingxue/SIGMA/blob/main/docs/reference_data.md)
- [Issues](https://github.com/Elsa-bingxue/SIGMA/issues)

## Scope

The package contains reusable SIGMA computation and standardized downstream
reporting. Large reference datasets, manuscript simulations, benchmarks, GO
enrichment and sample-specific interpretation remain in the analysis
repository.
