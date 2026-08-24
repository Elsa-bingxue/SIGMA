# Result generation

SIGMA separates computation from presentation. Python entry points produce the
versioned result tables and figures; notebooks remain short, inspectable views
of the same calls.

## One dataset

Copy one of the JSON templates in `examples/configs`, update the input and
output paths, and run:

```bash
sigma-omics run-config my_dataset.json
```

or:

```bash
python scripts/run_dataset.py my_dataset.json
```

Relative paths are resolved against the JSON file. The resolved paths, seed,
workflow and selection settings are written to
`analysis_config.resolved.json` in the result directory.

## Reference analyses

Manuscript cluster identities require the versioned ranking and assignment
tables distributed with the reference data. They are not silently inferred
from a fresh exploratory run. Use `run_reference_analysis` as documented in
`workflows_and_reproducibility.md`.

## Manuscript figures

The package stores the figure manifest and orchestration, while restricted or
large data remain in the reference project:

```bash
sigma-omics reproduce-manuscript /path/to/SIGMA-Code-Clean
sigma-omics reproduce-manuscript /path/to/SIGMA-Code-Clean --figures Fig2 Fig7
```

Fig. 2 contains the simulation, benchmark and core ablation presentation;
Fig. 7 contains the cross-cancer analysis. The underlying source tables should
be regenerated before plotting whenever the model or scientific defaults
change. A formatting-only release can reuse source tables if its provenance
matches the recorded package version and seed.

## Expected outputs

A complete dataset report contains the fitted core fields, program assignments,
signed-distance profiles, interface statistics, representative metabolite
table, influence-range summaries, anisotropy quality control, figures and run
provenance. Matched-ST validation outputs are written only when signature
scores are supplied.
