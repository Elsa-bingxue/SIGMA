import numpy as np

import pandas as pd

from sigma_spatial.selection import (
    side_near_far_statistics, select_bidirectional_programs,
    select_leading_program,
)


def test_bidirectional_statistics_preserve_named_side_semantics():
    d = np.r_[-10:-0:1, 1:11]
    # Program A rises toward the boundary on CI/negative; B does so on Cd/positive.
    scores = {
        "A": np.r_[-np.abs(d[:10]), np.zeros(10)],
        "B": np.r_[np.zeros(10), -np.abs(d[10:])],
    }
    table = side_near_far_statistics(
        scores, d,
        side_names={"negative": "CI / dopamine-low", "positive": "Cd / dopamine-high"},
        near_quantile=.2, far_quantile=.8, min_group_size=2,
    )
    assert set(table.side_name) == {"CI / dopamine-low", "Cd / dopamine-high"}
    ci = table[(table.program == "A") & (table.orientation == "negative")].iloc[0]
    cd = table[(table.program == "B") & (table.orientation == "positive")].iloc[0]
    assert ci.near_minus_far > 0
    assert cd.near_minus_far > 0


def test_selection_reports_both_sides_instead_of_positive_only():
    d = np.r_[-100:-0:1, 1:101]
    scores = {
        "CI_program": np.r_[-np.abs(d[:100]), np.zeros(100)],
        "Cd_program": np.r_[np.zeros(100), -np.abs(d[100:])],
    }
    table = side_near_far_statistics(scores, d, near_quantile=.2, far_quantile=.8,
                                     min_group_size=20)
    selected = select_bidirectional_programs(table, fdr_max=.05, effect_min=.1)
    assert set(selected.program) == {"CI_program", "Cd_program"}
    assert set(selected.orientation) == {"negative", "positive"}


def test_boundary_localized_selection_uses_two_stage_rule():
    statistics = pd.DataFrame({
        "program": [0, 1, 2, 3],
        "orientation": ["positive"] * 4,
        "side_name": ["disease"] * 4,
        "near_minus_far": [.61, .64, .83, .64],
        "fdr": [1e-20, 1e-15, 1e-30, 1e-16],
    })
    assignments = pd.DataFrame({
        "cluster": np.repeat([0, 1, 2, 3], 3),
        "interface_score": np.repeat([.55, .65, .33, .74], 3),
    })
    leading, audit = select_leading_program(
        statistics, assignments, orientation="positive",
        strategy="boundary_localized",
    )
    assert leading == 3
    assert audit.loc[audit.selected, "program"].item() == 3
    assert audit.passes_stage1.all()


def test_near_far_selection_preserves_original_rule():
    statistics = pd.DataFrame({
        "program": [0, 1, 2, 3],
        "orientation": ["positive"] * 4,
        "side_name": ["disease"] * 4,
        "near_minus_far": [.61, .64, .83, .64],
        "fdr": [1e-20, 1e-15, 1e-30, 1e-16],
    })
    assignments = pd.DataFrame({
        "cluster": np.repeat([0, 1, 2, 3], 2),
        "interface_score": np.repeat([.55, .65, .33, .74], 2),
    })
    leading, _ = select_leading_program(
        statistics, assignments, orientation="positive", strategy="near_far",
    )
    assert leading == 2
