import anndata as ad
import numpy as np
import pytest

from sigma_spatial.validation import validate_sigma_input, validate_sigma_output


def valid_adata():
    a = ad.AnnData(np.ones((4, 3)))
    a.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    a.obsm["X_harmony"] = np.arange(8, dtype=float).reshape(4, 2)
    a.obs["sigma_anchor"] = [1, np.nan, 0, np.nan]
    return a


def test_validate_sigma_input_accepts_binary_weak_anchors():
    assert validate_sigma_input(
        valid_adata(), anchor_key="sigma_anchor", representation_key="X_harmony"
    ).n_obs == 4


def test_validate_sigma_input_reports_missing_representation():
    a = valid_adata(); del a.obsm["X_harmony"]
    with pytest.raises(KeyError, match="SM-only training is not yet exposed"):
        validate_sigma_input(a, anchor_key="sigma_anchor", representation_key="X_harmony")


def test_validate_sigma_input_rejects_invalid_coordinates_and_anchors():
    a = valid_adata(); a.obsm["spatial"][0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        validate_sigma_input(a, anchor_key="sigma_anchor", representation_key="X_harmony")
    a = valid_adata(); a.obs["sigma_anchor"] = [1, np.nan, 1, np.nan]
    with pytest.raises(ValueError, match="both 0"):
        validate_sigma_input(a, anchor_key="sigma_anchor", representation_key="X_harmony")


def test_validate_sigma_output_checks_schema_and_probability_range():
    a = valid_adata()
    for key in ("sigma_gaussian_anchor", "sigma_region_probability_raw", "sigma_region_probability",
                "sigma_inside", "sigma_boundary", "sigma_d_signed"):
        a.obs[key] = np.zeros(a.n_obs)
    for key in ("X_sigma_msi", "X_sigma_gauss", "X_sigma_residual", "X_sigma_corrected"):
        a.obsm[key] = np.zeros((a.n_obs, 2))
    assert validate_sigma_output(a) is a
    a.obs["sigma_region_probability"] = 2
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        validate_sigma_output(a)
