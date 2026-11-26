# Installation

Use Python 3.11 in a new environment. The package name on PyPI is `sigma-omics`;
the Python import is `sigma_spatial`.

## PyPI

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --prefer-binary "sigma-omics==0.3.2"
python -m pip check
python -c "import sigma_spatial; print(sigma_spatial.__version__, sigma_spatial.__file__)"
```

On Windows, activate with `.venv\Scripts\activate`. PyTorch and PyTorch Geometric
are installed as dependencies. CPU execution is supported.
Version 0.3.2 passed clean-install checks on Intel macOS 15 and Ubuntu 24.04;
see the [test results](installation_test_0.3.2.md).

## GitHub

Create and activate a separate environment, then run:

```bash
git clone https://github.com/Elsa-bingxue/SIGMA.git
cd SIGMA
git checkout v0.3.2
python -m pip install --upgrade pip
python -m pip install --prefer-binary -e ".[dev]"
python -m pip check
```

The tag fixes the source version. To use the latest documentation and installation
tests, use the main branch. Model code is unchanged between these versions.

## Run the example

From the repository directory:

```bash
python examples/synthetic_demo.py --output results/synthetic_demo
```

The example writes `demo_sigma.h5ad`, programme membership and score tables,
signed-distance profiles, near/far statistics and figures. See
[running an analysis](running.md) for your own data and
[continuous distance](continuous_distance.md) for saved-field migration.

For source-code tests:

```bash
MPLBACKEND=Agg python -m pytest
```

Tests requiring unbundled historical datasets are skipped when those files are absent.

## Upgrade an existing installation

```bash
python -m pip install --upgrade --prefer-binary "sigma-omics==0.3.2"
```

Restart notebooks after upgrading. Check the version and import path to make sure
Python is loading the intended environment.

## Manuscript archive

From the extracted archive directory:

```bash
python -m pip install --prefer-binary Code_Release/public_0.3.2/dist/sigma_omics-0.3.2-py3-none-any.whl
python Code_Release/examples/synthetic_demo.py --output results/synthetic_demo
```

Dependencies may require internet access. The bundled installation files match the
public PyPI files and retain their recorded SHA-256 hashes.
