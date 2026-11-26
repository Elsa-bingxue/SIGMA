"""Small generated installation check for sigma-omics 0.3.2, not a benchmark."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import torch

from sigma_spatial import __version__, run_downstream_analysis, run_sigma


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('results/synthetic_demo'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)
    rng = np.random.default_rng(0)
    gx, gy = np.meshgrid(np.arange(20), np.arange(20))
    xy = np.column_stack([gx.ravel(), gy.ravel()]).astype(float)
    d = np.linalg.norm(xy - 9.5, axis=1) - 5.2
    templates = np.column_stack([
        np.exp(-np.abs(d) / 1.5),
        1 - np.exp(-np.abs(d) / 2),
        1 / (1 + np.exp(d)),
        1 / (1 + np.exp(-d)),
    ])
    x = np.column_stack([
        0.1 + (2 + rng.random()) * templates[:, j % 4]
        + rng.gamma(2, 0.08, len(d)) for j in range(80)
    ]).astype(np.float32)
    a = ad.AnnData(x)
    a.var_names = [f'{100 + j * 1.013:.4f}' for j in range(x.shape[1])]
    a.var['feature_type'] = 'SM'
    a.obsm['spatial'] = xy
    a.obsm['X_harmony'] = np.column_stack([
        templates[:, j % 4] + rng.normal(0, 0.1, len(d)) for j in range(8)
    ]).astype(np.float32)
    anchors = np.full(len(d), np.nan)
    anchors[d < -2] = 1
    anchors[d > 2] = 0
    a.obs['sigma_anchor'] = anchors

    fitted = run_sigma(
        a, distance_mode='continuous', random_state=0,
        n_components=8, n_neighbors=6, z_dim=8, epochs=50, device='cpu',
    )
    fitted.write_h5ad(args.output / 'demo_sigma.h5ad')
    report = run_downstream_analysis(
        fitted, output_dir=args.output, prefix='demo', matrix_source='X',
        selection_mode='standardized', report_level='standard',
        top_k=40, n_programs=4, n_bins=16, min_bin_points=3,
        min_group_size=8, random_state=0,
    )
    summary = {
        'package_version': __version__, 'example': 'generated installation check',
        'n_spots': a.n_obs, 'n_ions': a.n_vars, 'epochs': 50,
        'distance_mode': fitted.uns['sigma']['distance_mode'],
        'distance_units': 'synthetic grid units',
        'programmes': len(report.programs.scores),
    }
    (args.output / 'demo_run.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    print(f'Outputs: {args.output.resolve()}')


if __name__ == '__main__':
    main()
