"""Side-aware selection of spatial interface programs."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini--Hochberg correction without an extra dependency."""
    p = np.asarray(p_values, float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    values = p[valid]
    if not len(values):
        return out
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted); restored[order] = np.clip(adjusted, 0, 1)
    out[valid] = restored
    return out


def side_near_far_statistics(
    program_scores: Mapping[object, np.ndarray],
    signed_distance: np.ndarray,
    *,
    side_names: Mapping[str, str] | None = None,
    near_quantile: float = 0.20,
    far_quantile: float = 0.80,
    min_group_size: int = 10,
) -> pd.DataFrame:
    """Compare near and far values independently on both interface sides.

    Negative and positive refer only to signed-distance orientation. Biological
    names such as ``CI / dopamine-low`` and ``Cd / dopamine-high`` are supplied
    explicitly through ``side_names``.
    """
    if not 0 <= near_quantile < far_quantile <= 1:
        raise ValueError("Require 0 <= near_quantile < far_quantile <= 1")
    d = np.asarray(signed_distance, float)
    names = {"negative": "negative", "positive": "positive"}
    if side_names is not None:
        names.update(side_names)
    rows = []
    for orientation, side_mask in (
        ("negative", d < 0), ("positive", d > 0),
    ):
        finite_side = side_mask & np.isfinite(d)
        distances = np.abs(d[finite_side])
        if len(distances) < 2 * min_group_size:
            raise ValueError(f"Too few finite spots on {orientation} side")
        near_threshold = float(np.quantile(distances, near_quantile))
        far_threshold = float(np.quantile(distances, far_quantile))
        for program, values in program_scores.items():
            y = np.asarray(values, float)
            if y.shape != d.shape:
                raise ValueError(f"Program {program!r} score length does not match distance")
            valid = finite_side & np.isfinite(y)
            near = valid & (np.abs(d) <= near_threshold)
            far = valid & (np.abs(d) >= far_threshold)
            if near.sum() < min_group_size or far.sum() < min_group_size:
                effect = p_value = np.nan
            else:
                effect = float(y[near].mean() - y[far].mean())
                p_value = float(mannwhitneyu(y[near], y[far], alternative="two-sided").pvalue)
            rows.append({
                "program": program, "orientation": orientation,
                "side_name": names[orientation], "near_minus_far": effect,
                "p_value": p_value, "n_near": int(near.sum()), "n_far": int(far.sum()),
                "near_threshold": near_threshold, "far_threshold": far_threshold,
            })
    result = pd.DataFrame(rows)
    # The legacy HCC/PD notebooks tested the two biological sides separately.
    # Preserve that family of hypotheses instead of allowing one side to alter
    # the multiple-testing correction on the other side.
    result["fdr"] = np.nan
    for orientation in result["orientation"].unique():
        mask = result["orientation"].eq(orientation)
        result.loc[mask, "fdr"] = _bh_fdr(result.loc[mask, "p_value"].to_numpy())
    return result


def select_bidirectional_programs(
    statistics: pd.DataFrame,
    *,
    fdr_max: float = 0.05,
    effect_min: float = 0.0,
) -> pd.DataFrame:
    """Return significant positive near-versus-far programs on either side.

    The function intentionally reports candidates on both sides and does not
    collapse them into a single universal 'positive-side' winner.
    """
    required = {"program", "orientation", "side_name", "near_minus_far", "fdr"}
    missing = required.difference(statistics.columns)
    if missing:
        raise KeyError(f"Missing columns: {sorted(missing)}")
    selected = statistics[
        (statistics.fdr <= fdr_max) & (statistics.near_minus_far >= effect_min)
    ].copy()
    return selected.sort_values(
        ["orientation", "near_minus_far"], ascending=[True, False]
    ).reset_index(drop=True)
