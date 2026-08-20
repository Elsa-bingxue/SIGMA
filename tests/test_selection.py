import numpy as np

from sigma_spatial.selection import side_near_far_statistics, select_bidirectional_programs


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
