"""Shared, layout-safe plotting helpers for SIGMA analyses."""

from pathlib import Path
import string

import matplotlib as mpl
import matplotlib.pyplot as plt


FIGURE_SIZES = {
    "single": (3.35, 2.7),
    "single_square": (3.35, 3.35),
    "double": (7.0, 4.2),
    "double_square": (7.0, 7.0),
    "full_page": (7.0, 9.0),
    "slides_wide": (10.0, 4.8),
}


def set_publication_style(font_family="Arial", base_font_size=9):
    """Apply a consistent vector-friendly manuscript style."""
    mpl.rcParams.update({
        "font.family": font_family,
        "font.size": base_font_size,
        "axes.titlesize": base_font_size + 1,
        "axes.labelsize": base_font_size,
        "xtick.labelsize": base_font_size - 1,
        "ytick.labelsize": base_font_size - 1,
        "legend.fontsize": base_font_size - 1,
        "figure.titlesize": base_font_size + 2,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 300,
        "savefig.bbox": None,
    })


def figure_size(name="double"):
    """Return a standard manuscript figure size."""
    try:
        return FIGURE_SIZES[name]
    except KeyError as exc:
        raise ValueError(f"Unknown figure size {name!r}") from exc


def add_panel_labels(axes, labels=None, x=-0.12, y=1.06, fontsize=11):
    """Add consistent panel labels in axes coordinates."""
    axes = list(axes)
    labels = list(labels or string.ascii_uppercase[:len(axes)])
    for ax, label in zip(axes, labels):
        ax.text(x, y, label, transform=ax.transAxes, ha="left", va="top",
                fontsize=fontsize, fontweight="bold", clip_on=False)


def legend_outside(ax, *, loc="upper left", anchor=(1.02, 1.0), **kwargs):
    """Place a legend outside the data axes."""
    return ax.legend(loc=loc, bbox_to_anchor=anchor, borderaxespad=0,
                     frameon=False, **kwargs)


def shared_colorbar(fig, mappable, axes, *, label=None, location="right", **kwargs):
    """Create one colorbar for a collection of axes."""
    colorbar = fig.colorbar(mappable, ax=list(axes), location=location, **kwargs)
    if label:
        colorbar.set_label(label)
    return colorbar


def style_spatial_axis(ax):
    """Remove decorations that are not meaningful for tissue coordinates."""
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


def save_figure(fig, path, *, formats=("pdf", "svg"), dpi=300,
                transparent=False, close=False):
    """Save a figure in consistent vector and optional raster formats.

    ``path`` is a stem or a filename.  Layout is expected to be solved by
    constrained layout/GridSpec rather than repaired with ``bbox_inches``.
    """
    path = Path(path)
    stem = path.with_suffix("") if path.suffix else path
    stem.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for fmt in formats:
        output = stem.with_suffix(f".{fmt.lstrip('.')}")
        fig.savefig(output, dpi=dpi, transparent=transparent, facecolor=fig.get_facecolor())
        written.append(output)
    if close:
        plt.close(fig)
    return written
