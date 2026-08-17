import hashlib
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import torch

from sigma_spatial.weak_anchor import run_sigma_weak_anchor, weak_anchor_msi_embedding


ROOT = Path(__file__).resolve().parents[2]
HCC_REFERENCES = {
    ROOT / "HCC_P1/clean/results/HCC_P1_SIGMA_core.h5ad": {
        "sigma_region_probability": "1e19fda570e530e0f85a026b8fe58bf32b1a36c0f383a1cf5076eddb7c73f87f",
        "sigma_boundary": "18caebe6fe908fa8f9ebb2d67091e12a5185e681959a42524d067d6ae92f58ad",
        "sigma_d_signed": "1c23a19a9358cd2ded4597dd98c36c30ae587a17f304628d2e1ee9e60d7d4d5b",
    },
    ROOT / "HCC_P4/clean/results/HCC_P4_SIGMA_core.h5ad": {
        "sigma_region_probability": "2385f93f9a9500501943dc26ddc316da2a40c858408ae0b49ab30fd83bb89089",
        "sigma_boundary": "242c27d84996a389eb7f7a99d309a8d411fef5702155a061865f6ea034c0ebbd",
        "sigma_d_signed": "8964e56c7366e842ba838ddc05e88548b4a0591c5a8ed4dbaa2169444eaa4f3d",
    },
}


def digest(values, key):
    dtype = "u1" if key == "sigma_boundary" else "<f8"
    return hashlib.sha256(np.asarray(values).astype(dtype).tobytes()).hexdigest()


def test_frozen_hcc_p1_p4_outputs_are_unchanged():
    for path, expected in HCC_REFERENCES.items():
        a = ad.read_h5ad(path, backed="r")
        for key, checksum in expected.items():
            assert digest(a.obs[key], key) == checksum


def test_weak_anchor_embedding_preserves_already_log_switch():
    matrix = np.arange(30, dtype=float).reshape(6, 5)
    direct = weak_anchor_msi_embedding(matrix, n_components=3, random_state=0, already_log=True)
    logged = weak_anchor_msi_embedding(matrix, n_components=3, random_state=0, already_log=False)
    assert direct.shape == (6, 3)
    assert not np.allclose(direct, logged)


def test_weak_anchor_smoke_and_legacy_aliases():
    torch.set_num_threads(1)
    rng = np.random.default_rng(3)
    side = 5
    xx, yy = np.meshgrid(np.arange(side), np.arange(side))
    a = ad.AnnData(rng.gamma(2, 1, size=(side * side, 7)))
    a.obsm["spatial"] = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    anchor = np.full(a.n_obs, np.nan)
    anchor[xx.ravel() <= 0] = 0
    anchor[xx.ravel() >= 4] = 1
    a.obs["sigma_anchor"] = anchor
    result = run_sigma_weak_anchor(a, n_components=4, n_neighbors=4, epochs=1, random_state=0)
    assert result.uns["sigma"]["rna_loss"] is False
    np.testing.assert_allclose(result.obs["sigma_region_probability"], result.obs["p_tumor_gauss_spectral"])
    np.testing.assert_allclose(result.obs["sigma_d_signed"], result.obs["d_signed"])
    np.testing.assert_allclose(result.obsm["X_sigma_residual"], result.obsm["R_spectral"])


def test_weak_anchor_requires_both_anchor_classes():
    a = ad.AnnData(np.ones((4, 3)))
    a.obsm["spatial"] = np.arange(8).reshape(4, 2)
    a.obs["sigma_anchor"] = [1, 1, np.nan, np.nan]
    with pytest.raises(ValueError, match="must contain"):
        run_sigma_weak_anchor(a, epochs=1)
