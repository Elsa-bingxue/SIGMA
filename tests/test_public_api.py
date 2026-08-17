import anndata as ad
import numpy as np

from sigma_spatial.api import run_sigma


def test_run_sigma_preserves_input_and_returns_standard_schema(monkeypatch):
    a = ad.AnnData(np.ones((4, 3)))
    a.obsm["spatial"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    a.obsm["X_harmony"] = np.ones((4, 2))
    a.obs["sigma_anchor"] = [1, np.nan, 0, np.nan]

    def fake_fit(self, result, **kwargs):
        for key in ("sigma_gaussian_anchor", "sigma_region_probability_raw", "sigma_region_probability",
                    "sigma_inside", "sigma_boundary", "sigma_d_signed"):
            result.obs[key] = np.zeros(result.n_obs)
        for key in ("X_sigma_msi", "X_sigma_gauss", "X_sigma_residual", "X_sigma_corrected"):
            result.obsm[key] = np.zeros((result.n_obs, 2))
        result.uns["sigma"] = {"seed": self.seed}
        return self

    monkeypatch.setattr("sigma_spatial.api.SIGMA.fit", fake_fit)
    result = run_sigma(a, random_state=13)
    assert result is not a
    assert "sigma_boundary" not in a.obs
    assert "_sigma_public_anchor" not in result.obs
    assert result.uns["sigma"]["seed"] == 13
    assert result.uns["sigma"]["public_api"]["anchor_key"] == "sigma_anchor"
