# SIGMA

SIGMA identifies metabolic programmes around tumour–stroma interfaces by integrating
spatial metabolomics, tissue coordinates and pathological or weak anchors. Matched
spatial RNA provides an auxiliary representation when available.

The analysis asks where a programme changes, how far its signal extends, and which
parts of an interface differ. The main outputs are a learned tissue field, a signed
distance coordinate, ranked ions and spatial programmes with distance profiles.

[![PyPI](https://img.shields.io/pypi/v/sigma-omics)](https://pypi.org/project/sigma-omics/)

## Install

Python **3.10 or later**; Python **3.11** was used for release testing.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --prefer-binary "sigma-omics==0.3.2"
python -c "import sigma_spatial; print(sigma_spatial.__version__)"
```

Use `.venv\Scripts\activate` on Windows. The distribution is `sigma-omics`; import
it as `sigma_spatial`. CPU execution is supported. PyTorch and PyTorch Geometric are
installed as dependencies. See [installation and upgrades](docs/installation.md).

## Quick start

Provide an MSI `AnnData` with non-negative intensities in `X`, coordinates in
`obsm["spatial"]`, and anchors in `obs["sigma_anchor"]` (1: tumour, 0: outside,
NaN: unknown). The paired workflow also requires an aligned RNA representation in
`obsm["X_harmony"]`. The same spot IDs must refer to the same positions in both modalities.

```python
import scanpy as sc
from sigma_spatial import run_sigma, run_downstream_analysis

adata = sc.read_h5ad("sample.h5ad")
result = run_sigma(
    adata,
    anchor_key="sigma_anchor",
    representation_key="X_harmony",
    distance_mode="continuous",
    random_state=0,
)
result.write_h5ad("sample_sigma.h5ad")

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

This is a new exploratory fit. Reproducing a historical programme such as P3 requires
its saved ranking and membership table; a newly fitted cluster number is not a fixed
biological label. The paired model applies its own `log1p` transformation: do not
silently substitute scaled or already logged data for the documented input.

## Workflows and examples

| Input | Entry point | Example |
|---|---|---|
| Paired MSI and RNA with anchors | `run_sigma` | [Paired workflow](examples/quickstart.py) |
| MSI with anchors, without matched RNA | `run_sigma_weak_anchor` | [MSI-only workflow](examples/weak_anchor_usage.py) |
| Previously fitted probability field | `update_interface_coordinate` | [Coordinate migration](docs/continuous_distance.md) |
| Fitted coordinate and MSI | `run_downstream_analysis` | [Programme analysis](docs/running.md) |
| No downloaded dataset | Small generated example | [Synthetic demo](examples/synthetic_demo.py) |

Run the generated example from this directory after installing the package:

```bash
python examples/synthetic_demo.py --output results/synthetic_demo
```

It exercises fitting, continuous distance, programme discovery and report generation.
Its small generated dataset is an installation check, not the manuscript benchmark.
For real data, see [input requirements](docs/input_format.md) and
[running an analysis](docs/running.md).

## Outputs

| Output | Location |
|---|---|
| Learned tissue field | `obs["sigma_region_probability"]` |
| Signed distance and explicit side | `obs["sigma_d_signed"]`, `obs["sigma_inside"]` |
| Continuous contour segments | `uns["sigma_coordinate"]["segments"]` |
| Programme members, scores and profiles | `<prefix>_program_*.csv` |
| Same-side near/far comparisons | `<prefix>_interface_program_statistics.csv` |
| Maps, profiles and settings | Figure files and `<prefix>_downstream_summary.json` |

See the [complete output schema](docs/output_schema.md). Negative distance means
inside; positive means outside. Distance units follow the supplied coordinates.
Use `sigma_inside` to assign exact-zero positions to their tissue side.

## Changes in 0.3.2

- Adds continuous distance to the interpolated 0.5 contour.
- Preserves explicit tissue sides at zero distance.
- Updates coordinate aliases and keeps old columns when migrating a saved field.
- Marks overlapping near/far groups as unavailable rather than testing the same
  positions against themselves.

The original trainable model is retained. `distance_mode="legacy"` remains the
default for historical reproduction; the example above explicitly selects
`"continuous"`. Recompute rankings, profiles and programmes after changing a saved
coordinate. See [migration details](docs/continuous_distance.md).

## Data and reproducibility

Raw datasets and trained models are not distributed in the Python wheel. Data sources
and preparation notes are in the [reference-data guide](docs/reference_data.md).
Use explicit matrix sources and preserve the spot order, anchors, seeds and output
settings. Complete manuscript reproduction also requires the corresponding saved
rankings, memberships and reference-data files.

The public 0.3.2 wheel passed 93 tests; six tests requiring external historical data
were skipped. The final package test used existing scientific dependencies; a full
clean dependency installation was not completed on the release machine.

RNA used for fitting or anchor construction supplies cellular context, not an
independent validation of that fit. Spatial signal extent alone does not establish
diffusion, metabolic flux or tumour invasion.

## Source and support

- [GitHub repository](https://github.com/Elsa-bingxue/SIGMA)
- [PyPI 0.3.2](https://pypi.org/project/sigma-omics/0.3.2/)
- [Issue tracker](https://github.com/Elsa-bingxue/SIGMA/issues)
