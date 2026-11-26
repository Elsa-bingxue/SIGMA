# Result generation

Follow [Running SIGMA 0.3.2](running.md) for fitting and programme reports.
Use [workflow and reproduction notes](workflows_and_reproducibility.md) for saved
rankings and memberships.

JSON configuration runs use the repository helper:

```bash
python scripts/run_dataset.py examples/configs/direct_pathology.json
```

Edit the file paths before running. They resolve relative to the configuration
file. This helper is not part of the wheel command-line interface.
