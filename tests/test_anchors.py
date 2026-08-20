import anndata as ad
import numpy as np
import pytest

from sigma_spatial.anchors import (
    AnchorSpec, anchors_from_preset, build_anchors, summarize_anchors,
)


def test_listed_policy_preserves_unlisted_as_unknown():
    spec = AnchorSpec("test", "annotation", "direct_annotation", ("T",), ("N",), ("M",))
    anchors = build_anchors(["T", "N", "M", "other"], spec)
    np.testing.assert_equal(anchors, [1.0, 0.0, np.nan, np.nan])


def test_complement_policy_reproduces_legacy_gbm_semantics():
    spec = AnchorSpec("GBM", "region", "marker_informed", ("core",), negative_policy="complement")
    anchors = build_anchors(["core", "margin", "mixed", "neural"], spec)
    np.testing.assert_equal(anchors, [1.0, 0.0, 0.0, 0.0])


def test_preset_records_auditable_provenance():
    a = ad.AnnData(np.ones((5, 2)))
    a.obs["annotation"] = ["Tumor", "Normal cells", "Airway", "other", "Tumor"]
    anchors_from_preset(a, "HLC091", scheme="legacy")
    assert summarize_anchors(a.obs["sigma_anchor"]) == {
        "positive": 2, "negative": 2, "unknown": 1, "total": 5,
    }
    assert a.uns["sigma_anchor_provenance"]["scheme"] == "legacy"
    assert a.uns["sigma_anchor_provenance"]["source_column"] == "annotation"


@pytest.mark.parametrize(
    "dataset,column,labels,expected",
    [
        (
            "ccRCC_Y27T", "tumor_region_by_marker_in_GT",
            ["tumor_region", "TME_or_mixed_region", "tumor_transition", "mixed_region"],
            [1.0, 0.0, np.nan, np.nan],
        ),
        (
            "GBM", "gbm_region_by_marker_in_GT",
            ["tumor_core_like", "microenvironment_like", "neural_like",
             "invasive_margin_like", "mixed_transition"],
            [1.0, 0.0, 0.0, np.nan, np.nan],
        ),
    ],
)
def test_core_anchored_presets_leave_transition_unknown(dataset, column, labels, expected):
    a = ad.AnnData(np.ones((len(labels), 2)))
    a.obs[column] = labels
    anchors = anchors_from_preset(a, dataset, scheme="core_anchored")
    np.testing.assert_equal(anchors, expected)
    assert a.uns["sigma_anchor_provenance"]["scheme"] == "core_anchored"


@pytest.mark.parametrize(
    "dataset,column,labels",
    [
        ("ccRCC_Y27T", "tumor_region_by_marker_in_GT",
         ["tumor_region", "tumor_transition", "TME_or_mixed_region", "mixed_region"]),
        ("GBM", "gbm_region_by_marker_in_GT",
         ["tumor_core_like", "invasive_margin_like", "microenvironment_like", "mixed_transition"]),
    ],
)
def test_core_vs_rest_makes_every_noncore_label_negative(dataset, column, labels):
    a = ad.AnnData(np.ones((len(labels), 2)))
    a.obs[column] = labels
    anchors = anchors_from_preset(a, dataset, scheme="core_vs_rest")
    np.testing.assert_equal(anchors, [1.0, 0.0, 0.0, 0.0])


def test_invalid_label_overlap_is_rejected():
    spec = AnchorSpec("bad", "x", "direct_annotation", ("A",), ("A",))
    with pytest.raises(ValueError, match="disjoint"):
        build_anchors(["A", "B"], spec)


@pytest.mark.parametrize("dataset", ["HPD_A1", "HPD_B1", "HPD_C1"])
def test_pd_named_region_presets_keep_non_ci_cd_regions_unknown(dataset):
    a = ad.AnnData(np.ones((6, 2)))
    a.obs["RegionLoupe_str"] = ["Cd", "CI", "ACB", "NA", "unk", "nan"]
    anchors = anchors_from_preset(a, dataset, scheme="named_regions")
    np.testing.assert_equal(anchors, [0.0, 1.0, np.nan, np.nan, np.nan, np.nan])
    assert a.uns["sigma_anchor_provenance"]["mode"] == "named_regions"
