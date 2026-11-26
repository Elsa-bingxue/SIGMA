"""Read the small, saved section results supplied with these notebooks."""
from pathlib import Path
import hashlib
import json

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


def repository_root():
    for path in [Path.cwd(), *Path.cwd().parents]:
        if (path / 'notebooks/sections.json').exists():
            return path
    raise FileNotFoundError('Start Jupyter from the SIGMA repository or notebooks directory.')


def load_section(section):
    root = repository_root()
    config = json.loads((root / 'notebooks/sections.json').read_text())[section]
    folder = root / 'data/section_results' / section
    for name, expected in config['sha256'].items():
        actual = hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Saved file changed: ' + str(folder / name))
    with np.load(folder / 'display.npz', allow_pickle=False) as source:
        data = {key: source[key] for key in source.files}
    members = pd.read_csv(folder / 'members.csv')
    profiles = pd.read_csv(folder / 'profiles.csv')
    n = len(data['xy'])
    if data['scores'].shape != (n, len(config['programmes'])):
        raise ValueError('Score columns do not match programme labels.')
    if len(data['distance']) != n or not np.isfinite(data['xy']).all():
        raise ValueError('Invalid spatial coordinates.')
    return config, data, members, profiles


def plot_section(config, data, members, profiles):
    """Show every saved position; do not smooth spatial observations."""
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    cmap = LinearSegmentedColormap.from_list('programme', ['#91B98C', '#FBFAF7', '#86428E'])
    k = len(config['programmes'])
    fig, axes = plt.subplots(2, k + 1, figsize=(3.0 * (k + 1), 5.8),
                             gridspec_kw={'height_ratios': [2.0, 1.0]}, constrained_layout=True)
    xy = data['xy']
    distance = data['distance']
    limit = max(float(np.nanmax(abs(distance))), 1e-6)
    im = axes[0, 0].scatter(*xy.T, c=distance, s=8, linewidths=0,
                             cmap='RdBu_r', vmin=-limit, vmax=limit, rasterized=True)
    axes[0, 0].set_title(config['distance_label'])
    fig.colorbar(im, ax=axes[0, 0], shrink=.7, label=config['distance_units'])
    axes[1, 0].axis('off')
    axes[1, 0].text(0, .98, f"{len(xy):,} positions\n{len(members):,} listed ions\nSaved results; no refitting",
                    transform=axes[1, 0].transAxes, va='top', fontsize=10)
    for j, programme in enumerate(config['programmes']):
        count = int((members['programme'] == programme).sum())
        image = axes[0, j + 1].scatter(*xy.T, c=data['scores'][:, j], s=8,
                    linewidths=0, cmap=cmap, vmin=-2.5, vmax=2.5, rasterized=True)
        axes[0, j + 1].set_title(f'{programme} · {count} ions')
        subset = profiles[profiles.programme == programme].sort_values('distance')
        axes[1, j + 1].plot(subset.distance, subset['mean'], color='#207F73', lw=1.6)
        axes[1, j + 1].axvline(0, color='#7C8589', lw=.7, ls='--')
        axes[1, j + 1].axhline(0, color='#D8DDDE', lw=.6)
        axes[1, j + 1].set_xlabel('Signed distance')
        axes[1, j + 1].set_ylabel('Mean score')
    for axis in axes[0]:
        axis.set_aspect('equal')
        axis.axis('off')
    fig.colorbar(image, ax=list(axes[0, 1:]), shrink=.7, label='Programme score (z)', extend='both')
    fig.suptitle(config['title'], fontsize=14, fontweight='bold')
    return fig
