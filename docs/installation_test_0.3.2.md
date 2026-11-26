# Installation test: 0.3.2

Tested on 26 November 2026. [Run logs and outputs](https://github.com/Elsa-bingxue/SIGMA/actions/runs/45824514430) are available on GitHub Actions.

## Environments

Each job started on a new hosted virtual machine and created a Python 3.11
environment without system packages or a dependency cache. Before installation,
the environments contained only pip and setuptools.

| Platform | Installation | Dependency check | Workflow | Targeted tests |
|---|---|---|---|---|
| macOS 15, Intel | PyPI 0.3.2 | Passed | Passed | 20 passed |
| macOS 15, Intel | GitHub source | Passed | Passed | 20 passed |
| Ubuntu 24.04, x86-64 | PyPI 0.3.2 | Passed | Passed | 20 passed |
| Ubuntu 24.04, x86-64 | GitHub source | Passed | Passed | 20 passed |

The GitHub installations used commit
`782407119f0a725d1bfdcc0ae60d22bb005d5f94`, the source tagged `v0.3.2`.
PyPI installations downloaded the published package and dependencies directly
from the package index. All analysis ran on CPU.

## Continuous distance

For a probability field with a known straight 0.5 contour, the largest absolute
distance error was `4.44e-16` coordinate units in every environment. Exact-zero
positions retained their inside/outside assignment. Coordinate export and H5AD
reload also passed.

## Near/far comparisons

Same-side comparisons recovered the two prescribed effects, 5.5 and 10.5.
When tied thresholds caused near and far groups to overlap, the overlap count
was reported and the effect, P value and FDR were unavailable.

## Pattern analysis

Each environment ran the README example with 400 generated spots, 80 ions and
50 training epochs. A separate run used the paired model's default 1,000 epochs
with continuous distance and the default downstream analysis. All runs produced
the expected programme maps, profiles, membership tables, near/far statistics
and H5AD files. The default runs assigned all 80 ions to four programmes and
reported eight same-side comparisons.

Within each platform, PyPI and GitHub installations produced identical programme
membership after matching cluster labels (adjusted Rand index = 1). Their
near/far group counts and effects also agreed, within floating-point precision.
The two macOS distance fields were identical. The independent Ubuntu fits
differed by at most 0.0383 generated-coordinate units, with correlation 0.999994.
Thus the learned fields were not identical in every run, despite matching
programme membership in this example.

## Scope

These checks establish installation and execution on the listed CPU environments.
They do not test Windows, Apple Silicon, GPU execution or performance on patient
data. Dependency versions, input packages, import paths and output checks are
recorded in the run artifacts. The 20 targeted tests cover continuous coordinates,
coordinate migration and programme analysis; they are not the full test suite.
