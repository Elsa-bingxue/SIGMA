import numpy as np

from sigma_spatial.simulation import simulate_interface_features


def test_simulation_is_seeded_and_retains_reference_return_type():
    distance = np.linspace(0, 10, 12)
    coords = np.column_stack([distance, np.zeros_like(distance)])
    kwargs = dict(d_abs=distance, coords=coords, n_positive=2, n_negative=1,
                  n_region_only=1, n_spatial_background=1, n_noise=1,
                  effect_size=0.5, noise_sd=0.2, lam=100, random_state=7)
    x1, meta1 = simulate_interface_features(**kwargs)
    x2, meta2 = simulate_interface_features(**kwargs)
    np.testing.assert_array_equal(x1, x2)
    assert meta1 == meta2
    assert isinstance(meta1, list)  # Preserve HBC_515 cell 46 behaviour.
    assert x1.shape == (12, 6)
    assert [row["sim_type"] for row in meta1] == [
        "positive_boundary", "positive_boundary", "negative_boundary",
        "region_only", "spatial_background", "noise",
    ]
