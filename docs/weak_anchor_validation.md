# SM-only weak-anchor validation

The public weak-anchor module preserves the executed HCC P1/P4 notebook path:
64-component MSI SVD, directed Gaussian kNN graph (`k=15`), 50-step anchor
diffusion (`alpha=0.85`), 10-step embedding diffusion (`alpha=0.9`), and the
1000-epoch loss `loss_rec + 0.05 * loss_sup`. No RNA loss is used.

Full CPU reruns with seed 0 produced:

| Sample | Probability correlation | Probability MAE | Inside agreement | Boundary Dice | Signed-distance correlation |
|---|---:|---:|---:|---:|---:|
| HCC P1 | 0.999407 | 0.006984 | 0.999093 | 0.978224 | 0.966162 |
| HCC P4 | 0.999772 | 0.004785 | 1.000000 | 1.000000 | 1.000000 |

The MSI embeddings were identical to the frozen outputs. Small learned-field
differences are attributable to the current PyTorch/PyG numerical environment;
P4 boundary and distance outputs were identical, while P1 changed side for four
of 4,409 spots near the 0.5 decision level.

Detailed values are stored in
`cross_cancer/results/nature_methods/tables/HCC_weak_anchor_rerun_validation.csv`.
