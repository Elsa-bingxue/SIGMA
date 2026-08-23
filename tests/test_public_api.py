import anndata as ad
import numpy as np
import pandas as pd

from sigma_spatial.api import (
    DownstreamAnalysis, run_analysis, run_downstream_analysis, run_sigma,
)


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


def test_downstream_analysis_writes_complete_report(tmp_path):
    rng = np.random.default_rng(4)
    n = 320
    distance = np.linspace(-12, 12, n)
    patterns = np.column_stack([
        np.exp(-np.abs(distance) / scale) + rng.normal(0, .03, n)
        for scale in (1.2, 2.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0)
    ])
    obj = ad.AnnData(patterns)
    obj.obs["sigma_d_signed"] = distance
    obj.obsm["spatial"] = np.column_stack([np.arange(n), np.zeros(n)])
    obj.var["feature_type"] = "SM"
    ranking = pd.DataFrame({
        "j": np.arange(patterns.shape[1]),
        "mz": np.arange(patterns.shape[1]) + 100.0,
        "candidate_score": np.arange(patterns.shape[1], 0, -1),
    })
    signatures = {"hypoxia_score": np.exp(-np.abs(distance) / 2)}
    result = run_downstream_analysis(
        obj, ranking, tmp_path, prefix="demo", signature_scores=signatures,
        top_k=8, n_programs=4, n_bins=20, min_bin_points=3,
        min_group_size=10, anisotropy_min_points=5,
    )
    assert result.leading_program in result.programs.scores
    assert len(result.representatives) <= 5
    for name in (
        "demo_program_assignments.csv",
        "demo_interface_program_statistics.csv",
        "demo_representative_metabolites.csv",
        "demo_metabolic_patterns_and_distance_profiles.pdf",
        "demo_interface_program_enrichment.pdf",
        "demo_st_program_heatmap.pdf",
        "demo_st_interface_lollipop.pdf",
        "demo_downstream_summary.json",
    ):
        assert (tmp_path / name).exists(), name


def test_run_analysis_routes_precomputed_sigma_from_explicit_evidence(tmp_path, monkeypatch):
    obj = ad.AnnData(np.ones((8, 3)))
    obj.obsm["spatial"] = np.column_stack([np.arange(8), np.zeros(8)])
    obj.var["feature_type"] = "SM"
    obj.obs["sigma_region_probability"] = np.linspace(0, 1, 8)
    obj.obs["sigma_boundary"] = False
    obj.obs["sigma_d_signed"] = np.linspace(-2, 2, 8)

    captured = {}
    expected = DownstreamAnalysis(None, None, None, 2, pd.DataFrame())

    def fake_downstream(result, **kwargs):
        captured.update(kwargs)
        return expected

    monkeypatch.setattr("sigma_spatial.api.run_downstream_analysis", fake_downstream)
    result = run_analysis(
        obj, tmp_path, workflow="auto",
        evidence={"same_section_pathology": True}, report_level="none",
    )
    assert result.workflow == "direct_pathology"
    assert result.downstream is expected
    assert captured["selection_mode"] == "lambda_profile"
    assert captured["program_selection"] == "auto"
    assert captured["report_level"] == "none"
    assert (tmp_path / "sigma_run_config.json").exists()


def test_run_analysis_auto_does_not_guess_evidence(tmp_path):
    obj = ad.AnnData(np.ones((4, 2)))
    obj.obsm["spatial"] = np.ones((4, 2))
    with np.testing.assert_raises_regex(ValueError, "explicit evidence"):
        run_analysis(obj, tmp_path, workflow="auto")
