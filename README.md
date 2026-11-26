# SIGMA

SIGMA identifies metabolic programmes around tumour–stroma interfaces. It uses
spatial MSI, tissue coordinates and pathological or weak anchors; matched RNA is
included in the paired workflow. Ions are grouped by their distance profiles to
show where metabolic signals change and how far they extend within the tissue.

[![PyPI](https://img.shields.io/pypi/v/sigma-omics)](https://pypi.org/project/sigma-omics/)

## Install

Use Python 3.11 in a separate environment:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --prefer-binary "sigma-omics==0.3.2"
python -m pip check
python -c "import sigma_spatial; print(sigma_spatial.__version__)"
```

On Windows, activate with `.venv\Scripts\activate`. The distribution is
`sigma-omics`; the import name is `sigma_spatial`. CPU execution is supported.
See [installation notes](docs/installation.md) for upgrades and source installation.

## Inputs

Prepare an `AnnData` object with the following fields:

| Field | Contents |
|---|---|
| `X` | Non-negative MSI intensities, before log transformation |
| `obsm["spatial"]` | Two-dimensional tissue coordinates |
| `obs["sigma_anchor"]` | 1: tumour; 0: outside; NaN: unlabelled |
| `obsm["X_harmony"]` | Aligned RNA representation for the paired workflow |

Registration and RNA preprocessing are performed before running SIGMA. Rows must
refer to the same tissue positions across modalities. The paired model applies
`log1p` to the MSI matrix. See [input format](docs/input_format.md).

## Run

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

The first call fits the interface. The second ranks ions, groups their distance
profiles and writes programme maps, profiles, membership tables and same-side
near/far comparisons. Parameters and output status are saved with the results.
See [running an analysis](docs/running.md) and [output fields](docs/output_schema.md).

## Examples

| Data | Example |
|---|---|
| Generated data, no download required | [Synthetic example](examples/synthetic_demo.py) |
| Paired MSI and RNA | [Paired analysis](examples/quickstart.py) |
| MSI with anchors, without matched RNA | [MSI-only analysis](examples/weak_anchor_usage.py) |
| A saved SIGMA probability field | [Update the distance coordinate](docs/continuous_distance.md) |

After downloading or cloning this repository, run the small example:

```bash
python examples/synthetic_demo.py --output results/synthetic_demo
```

It uses 400 generated spots and 80 ions to check fitting and report generation.
The full example and the default 1,000-epoch workflow are tested separately.
See the [0.3.2 installation results](docs/installation_test_0.3.2.md).

## Continuous distance and near/far comparisons

Version 0.3.2 adds distance to the interpolated 0.5 contour. Select it with
`distance_mode="continuous"`; the default remains `"legacy"` for earlier analyses.
Negative distance denotes the inside and positive distance the outside. Units
follow the supplied coordinates. Exact-zero positions retain their side through
`obs["sigma_inside"]`.

Near/far comparisons are made within each side. If tied thresholds produce
overlapping groups, the effect, P value and FDR are recorded as unavailable.
Updating a saved coordinate preserves the earlier columns; rankings, profiles
and programme assignments must then be recalculated. See the
[migration guide](docs/continuous_distance.md) and [change log](CHANGELOG.md).

## Data and reproduction

[Data sources](docs/reference_data.md) lists the cancer datasets used in the study.
Raw data and trained models are obtained separately. To reproduce a historical
programme such as P3, use its saved ranking and membership table. New fits can
produce different cluster labels.

RNA used for fitting provides cellular context. Independent tissue measurements
are needed to validate biological interpretations. Spatial signal extent is not
a measurement of diffusion or metabolic flux.

## Links

- [PyPI release](https://pypi.org/project/sigma-omics/0.3.2/)
- [Workflow and reproduction notes](docs/workflows_and_reproducibility.md)
- [Issues](https://github.com/Elsa-bingxue/SIGMA/issues)
