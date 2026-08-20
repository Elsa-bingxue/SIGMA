import anndata as ad
import numpy as np
import pytest

from sigma_spatial.assessment import assess_dataset, prepare_anchors


def example():
    obj = ad.AnnData(np.ones((6, 3)))
    obj.obsm["spatial"] = np.c_[np.arange(6), np.zeros(6)]
    obj.obs["region"] = ["tumor", "tumor", "stroma", "stroma", "mixed", None]
    return obj


@pytest.mark.parametrize("evidence,workflow,policy", [
    ({"same_section_pathology": True}, "direct_pathology", "direct_binary"),
    ({"matched_st": True}, "multiomics_inferred", "high_confidence_three_state"),
    ({"adjacent_section_pathology": True}, "transferred_pathology", "cluster_assisted_transfer"),
    ({"region_defined": True}, "region_defined_disease", "named_regions"),
])
def test_assessment_recommends_from_declared_evidence(evidence, workflow, policy):
    obj = example()
    result = assess_dataset(obj, annotation_key="region", evidence=evidence)
    assert result.recommended_workflow == workflow
    assert result.recommended_anchor_policy == policy
    assert result.ready


def test_assessment_does_not_guess_missing_evidence():
    with pytest.raises(ValueError, match="evidence"):
        assess_dataset(example(), annotation_key="region", evidence={})


def test_prepare_anchors_keeps_unlisted_and_missing_values_unknown():
    obj = example()
    anchors = prepare_anchors(
        obj, annotation_key="region", positive_labels=["tumor"],
        negative_labels=["stroma"], unknown_labels=["mixed", "None", "nan"],
        workflow="direct_pathology", anchor_policy="direct_binary",
    )
    np.testing.assert_equal(anchors, [1, 1, 0, 0, np.nan, np.nan])
    assert obj.uns["sigma_workflow_provenance"]["workflow"] == "direct_pathology"


def test_prepare_anchors_copy_does_not_modify_input():
    obj = example()
    result = prepare_anchors(
        obj, annotation_key="region", positive_labels=["tumor"],
        negative_labels=["stroma"], workflow="direct_pathology",
        anchor_policy="direct_binary", copy=True,
    )
    assert "sigma_anchor" not in obj.obs
    assert "sigma_anchor" in result.obs
