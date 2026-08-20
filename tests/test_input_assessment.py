import anndata as ad
import numpy as np

from sigma_spatial import inspect_input


def _object(types):
    obj = ad.AnnData(np.ones((8, len(types))))
    obj.var["feature_type"] = types
    obj.obsm["spatial"] = np.column_stack([np.arange(8), np.zeros(8)])
    return obj


def test_sm_only_routes_to_weak_anchor_without_claiming_st():
    report = inspect_input(_object(["SM", "SM"]))
    assert report.modalities == ("SM",)
    assert report.recommended_entry_point == "run_sigma_weak_anchor"


def test_joint_routes_to_main_sigma():
    report = inspect_input(_object(["ST", "SM"]))
    assert report.modalities == ("SM", "ST")
    assert report.recommended_entry_point == "run_sigma"


def test_st_only_does_not_claim_metabolite_programs():
    report = inspect_input(_object(["ST", "ST"]))
    assert report.recommended_entry_point == "st_interface_validation"
    assert report.supports_metabolic_sigma is False


def test_ambiguous_matrix_is_not_guessed():
    obj = ad.AnnData(np.ones((8, 2)))
    obj.obsm["spatial"] = np.ones((8, 2))
    report = inspect_input(obj)
    assert report.recommended_entry_point == "declare_modality"
