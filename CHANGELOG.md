# Changelog

## 0.3.2 — 26 November 2026

This patch adds continuous interface coordinates and corrects side-specific near/far comparisons. The original SIGMA training architecture and loss remain unchanged.

- Add `continuous_signed_distance` and `update_interface_coordinate`. Distances are measured to an interpolated contour, with long triangulation edges masked to avoid bridging unsampled gaps.
- Add `distance_mode="continuous"` to paired and weak-anchor workflows. The default remains `"legacy"` for reproducibility.
- Store tissue side separately from distance. Zero-distance observations use an explicit inside mask; without one they are excluded from side-specific comparisons with a warning.
- Preserve original coordinate columns when upgrading a saved field and synchronize existing `d_signed` and `d_abs` aliases. Downstream ranking and programmes must be recomputed explicitly.
- Mark near/far statistics as unavailable when tied distance thresholds put observations in both groups. Programme statistics include `n_overlap`.
- Apply the overlap guard to ion ranking and reject overlapping fixed-width near/far windows.
- Update installation, coordinate migration and output documentation; add regression and installation tests.

This patch does not introduce experimental raw-MSI/logistic model branches or change reference-dataset quantile presets.
