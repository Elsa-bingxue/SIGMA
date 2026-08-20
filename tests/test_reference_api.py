import pytest

from sigma_spatial import get_reference_preset, list_reference_datasets


def test_reference_registry_contains_all_current_datasets():
    assert len(list_reference_datasets()) == 11
    assert get_reference_preset("HBC515").leading_program == 3
    assert get_reference_preset("HCC_P1").near_quantile == .25
    assert get_reference_preset("HPD_A1").ranking_col == "r2_logI"


def test_unknown_reference_dataset_is_rejected():
    with pytest.raises(KeyError, match="No reference preset"):
        get_reference_preset("unknown")
