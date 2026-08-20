import anndata as ad
import numpy as np
import pytest

from sigma_spatial.workflows import (
    DATASET_WORKFLOWS, VALID_WORKFLOWS, WorkflowConfig,
    get_workflow_config, record_workflow_provenance,
)


def test_registry_contains_all_four_evidence_source_workflows():
    assert {config.workflow for config in DATASET_WORKFLOWS.values()} == set(VALID_WORKFLOWS)


def test_dataset_evidence_sources_are_not_conflated():
    assert get_workflow_config("ccRCC_Y27T").workflow == "multiomics_inferred"
    assert get_workflow_config("ccRCC_Y27T").matched_st is True
    assert get_workflow_config("HCC_P1").workflow == "transferred_pathology"
    assert get_workflow_config("HCC_P1").matched_st is False
    assert get_workflow_config("HPD_A1").workflow == "region_defined_disease"


def test_workflow_provenance_is_metadata_only():
    a = ad.AnnData(np.arange(6, dtype=float).reshape(3, 2))
    before = a.X.copy()
    payload = record_workflow_provenance(a, "HCC_P1")
    np.testing.assert_array_equal(a.X, before)
    assert payload["anchor_policy"] == "cluster_assisted_transfer"
    assert a.uns["sigma_workflow_provenance"]["annotation_source"].startswith("adjacent-section")


def test_invalid_workflow_name_is_rejected():
    with pytest.raises(ValueError, match="workflow must be"):
        WorkflowConfig("x", "unknown", "x", "x", "x", "x", "x")
