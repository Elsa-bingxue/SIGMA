# Section notebooks

Start Jupyter from the repository directory with `jupyter lab notebooks`.
Each notebook shows saved maps, a member table and distance profiles. The small
files in `data/section_results/` are included, so this viewing mode needs no raw-data
download and does not train a model. File checksums are verified on loading.

## Recalculate

The final cell is optional. Set `INPUT_H5AD` to a full fitted AnnData object and
`RUN_ANALYSIS=True`. It uses sigma-omics 0.3.2 to update continuous distance and
rerun programme analysis. Results go to `results/<section>/continuous/`; the saved
preview files are untouched. To fit from prepared input data, follow the
[analysis guide](../docs/running.md).

The preview coordinate is labelled in each notebook. Breast previews use the
common pathology-reference coordinate from Figure 3. Lung previews use saved
continuous SIGMA distances. Other archived previews retain their original
distance values. Viewing an archived result does not apply the coordinate fix.

## Interpretation

Programme IDs are specific to a section and fit. BC515 S1/S2 are two sections from
the same case. CRC B_2/C_2 are tissue regions. Glioma574 and SMINT P1 are kept
separate from GBM, and their independence as patients is unverified.

These notebooks show section-level programmes. Current ROI contrasts, held-out
validation and composition-adjusted analyses remain in the manuscript's source
data and supplementary material. No PD examples are included in this index.
