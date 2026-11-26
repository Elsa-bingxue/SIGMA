# Data

`section_results/` contains compact saved results for the section notebooks:

- `display.npz`: coordinates, the labelled distance coordinate and programme scores.
- `members.csv`: molecular membership for the displayed programmes.
- `profiles.csv`: mean score and spot count in 32 equal-width distance bins.

The files do not contain the full MSI/RNA matrices or trained models. Sources and
raw-data requirements are listed in [reference data](../docs/reference_data.md).
SHA-256 hashes and section-specific notes are in `notebooks/sections.json`.
