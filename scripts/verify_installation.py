"""Check an installed release independently of its source-tree tests."""
from __future__ import annotations

import argparse
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import site
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch
import sigma_spatial as sigma
from sigma_spatial.coordinates import interface_side_masks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    torch.set_num_threads(1)
    assert sys.prefix != sys.base_prefix, 'Not a virtual environment'
    cfg = (Path(sys.prefix) / 'pyvenv.cfg').read_text().lower()
    assert 'include-system-site-packages = false' in cfg
    assert not site.ENABLE_USER_SITE
    before = json.loads((root / 'packages_before.json').read_text())
    assert {p['name'].lower() for p in before} <= {'pip', 'setuptools', 'wheel'}
    assert sigma.__version__ == '0.3.2'
    assert Path(metadata.distribution('sigma-omics')._path).is_relative_to(Path(sys.prefix))
    source = os.environ['INSTALL_SOURCE']
    imported = Path(sigma.__file__).resolve()
    if source == 'pypi':
        assert imported.is_relative_to(Path(sys.prefix)), imported
    else:
        expected = Path(os.environ['RUNNER_TEMP']) / 'sigma-source/src'
        assert imported.is_relative_to(expected), imported
    evidence = {
        'install_source': source, 'python': sys.version, 'platform': platform.platform(),
        'architecture': platform.machine(), 'prefix': sys.prefix,
        'base_prefix': sys.base_prefix, 'import_path': str(imported),
        'version': sigma.__version__, 'clean_environment': True,
    }

    # A line between sampled columns has known distances in coordinate units.
    x, y = np.meshgrid(np.arange(-3., 4), np.arange(-3., 4))
    xy = np.column_stack([x.ravel(), y.ravel()])
    field = (xy[:, 0] + 4) / 8
    d, inside, segments = sigma.continuous_signed_distance(xy, field, level=.55)
    np.testing.assert_allclose(d, .4 - xy[:, 0], atol=1e-12)
    np.testing.assert_array_equal(inside, field >= .55)
    scaled, _, _ = sigma.continuous_signed_distance(xy * 21, field, level=.55)
    np.testing.assert_allclose(scaled, d * 21, atol=1e-12)
    evidence['analytic_distance_max_error'] = float(np.max(np.abs(d - (.4 - xy[:, 0]))))

    zeros = np.array([-2., -0., 0., 2.])
    masks = interface_side_masks(zeros, [True, True, False, False])
    np.testing.assert_array_equal(masks['negative'], [True, True, False, False])
    np.testing.assert_array_equal(masks['positive'], [False, False, True, True])
    distances = np.array([-4., -3., -2., -1., -0., 0., 1., 2., 3., 4.])
    values = np.array([0., 0., 0., 1., 10., 20., 1., 0., 0., 0.])
    nearfar = sigma.side_near_far_statistics(
        {'known': values}, distances, inside_mask=np.arange(10) < 5,
        near_quantile=.25, far_quantile=.75, min_group_size=2,
    ).set_index('orientation')
    assert nearfar.loc['negative', 'near_minus_far'] == 5.5
    assert nearfar.loc['positive', 'near_minus_far'] == 10.5
    nearfar.to_csv(root / 'known_near_far.csv')
    tied = sigma.side_near_far_statistics(
        {'known': np.arange(40.)}, np.r_[-np.ones(20), np.ones(20)], min_group_size=2,
    )
    assert tied.n_overlap.eq(20).all()
    assert tied[['near_minus_far', 'p_value', 'fdr']].isna().all().all()
    tied.to_csv(root / 'overlap_guard.csv', index=False)
    evidence.update(zero_side_preserved=True, near_far_known_effects=[5.5, 10.5],
                    overlapping_groups_unavailable=True)

    # Run the README calls with their default 1,000-epoch fit, using the generated
    # input rather than a private dataset. Clear earlier outputs before refitting.
    demo = ad.read_h5ad(root / 'demo/demo_sigma.h5ad')
    a = ad.AnnData(demo.X.copy(), obs=demo.obs[['sigma_anchor']].copy(), var=demo.var.copy())
    a.obsm['spatial'] = demo.obsm['spatial'].copy()
    a.obsm['X_harmony'] = demo.obsm['X_harmony'].copy()
    fitted = sigma.run_sigma(a, anchor_key='sigma_anchor', representation_key='X_harmony',
                             distance_mode='continuous', random_state=0)
    assert np.isfinite(fitted.obs.sigma_d_signed).all()
    assert fitted.obs.sigma_inside.nunique() == 2
    assert fitted.uns['sigma_coordinate']['segments'].size > 0
    side = fitted.obs.sigma_inside.to_numpy(bool)
    dist = fitted.obs.sigma_d_signed.to_numpy()
    assert (dist[side] <= 0).all() and (dist[~side] >= 0).all()
    assert np.count_nonzero(dist) > .95 * len(dist)
    result_dir = root / 'readme_default'
    result_dir.mkdir(exist_ok=True)
    fitted.write_h5ad(result_dir / 'sample_sigma.h5ad')
    reread = ad.read_h5ad(result_dir / 'sample_sigma.h5ad')
    np.testing.assert_array_equal(reread.obs.sigma_inside, fitted.obs.sigma_inside)
    np.testing.assert_allclose(reread.obs.sigma_d_signed, fitted.obs.sigma_d_signed)
    report = sigma.run_downstream_analysis(
        reread, output_dir=result_dir, prefix='sample', matrix_source='X',
        selection_mode='standardized', report_level='standard', random_state=0,
    )
    assert len(report.programs.scores) == 4
    for score in report.programs.scores.values():
        assert len(score) == a.n_obs and np.isfinite(score).all()
    assignments = pd.read_csv(result_dir / 'sample_program_assignments.csv')
    assert assignments['cluster'].nunique() == 4
    assert len(assignments) == a.n_vars
    stats = pd.read_csv(result_dir / 'sample_interface_program_statistics.csv')
    assert stats['n_overlap'].eq(0).all()
    assert np.isfinite(stats['near_minus_far']).all()
    for name in ['sample_program_scores.csv', 'sample_program_centroid_profiles.csv',
                 'sample_metabolic_patterns_and_distance_profiles.pdf',
                 'sample_interface_program_enrichment.pdf', 'sample_downstream_summary.json']:
        assert (result_dir / name).stat().st_size > 0, name
    evidence.update(readme_default_fit_passed=True, h5ad_roundtrip_passed=True,
                    n_spots=a.n_obs, n_ions=a.n_vars, n_programmes=4,
                    assigned_ions=len(assignments), same_side_comparisons=len(stats),
                    pattern_outputs_complete=True)
    (root / 'verification.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(evidence, indent=2))


if __name__ == '__main__':
    main()
