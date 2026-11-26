# Installation: sigma-omics 0.3.2

## New environment

Use Python 3.10 or later; Python 3.11 was used for release testing.

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --prefer-binary "sigma-omics==0.3.2"
python -m pip check
python -c "import sigma_spatial; print(sigma_spatial.__version__, sigma_spatial.__file__)"
```

On Windows use `.venv\Scripts\activate`. No GPU is required. PyTorch and PyTorch
Geometric are dependencies. For CUDA, install the appropriate PyTorch build for the
machine before installing SIGMA. Do not mix packages from unrelated environments.

## Upgrade

```bash
python -m pip install --upgrade --prefer-binary "sigma-omics==0.3.2"
```

Restart active notebooks and check the version and import path. If an old editable
checkout is loaded, use a new environment or uninstall that checkout first.

## Run

Follow the [quick analysis guide](running.md). From the repository root:

```bash
python examples/synthetic_demo.py --output results/synthetic_demo
```

In the manuscript archive, use `Code_Release/examples/synthetic_demo.py` instead.
The expected files include `demo_sigma.h5ad`, programme tables, maps and profiles.

## Source install

From a checked-out source directory or the extracted 0.3.2 source archive:

```bash
python -m pip install --prefer-binary -e ".[dev]"
MPLBACKEND=Agg python -m pytest
```

External-data preservation tests skip when their historical datasets are absent.
Building a local wheel is optional: `python -m build`, then
`python -m twine check dist/*`. Local builds are not the verified PyPI artifacts.

## Bundled installation files

From the manuscript archive root:

```bash
python -m pip install --prefer-binary Code_Release/public_0.3.2/dist/sigma_omics-0.3.2-py3-none-any.whl
```

Dependencies may still require internet access. The supplied wheel and source
archive are the original public files; their SHA-256 hashes are recorded alongside
them. Their historical compatibility code is preserved unchanged.

## Distance migration

The default remains `legacy`. Select `distance_mode="continuous"` for a new fit,
or use `update_interface_coordinate` for a saved probability field. See
[continuous distance](continuous_distance.md). Changing the coordinate does not
recompute old rankings or programme membership.

`--prefer-binary` selects compatible prebuilt dependencies when available. An earlier
Intel macOS attempt compiled Numba/LLVM and failed; this option avoids that compilation
when compatible wheels exist. It is not a guarantee for every operating system.
