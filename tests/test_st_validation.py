import anndata as ad
import numpy as np
import pytest

from sigma_spatial.st_validation import (
    program_signature_association, record_st_validation_provenance,
    side_near_far_signature_statistics,
)


def test_program_signature_association_detects_aligned_signal():
    x = np.linspace(-2, 2, 80)
    table = program_signature_association({0: x, 1: -x}, {"hypoxia_score": x})
    assert table.loc[table.program.eq(0), "correlation"].iloc[0] > .99
    assert table.loc[table.program.eq(1), "correlation"].iloc[0] < -.99


def test_signature_near_far_retains_biological_side_names():
    distance = np.r_[-100:-0, 1:101]
    scores = {"dopamine_score": -np.abs(distance)}
    table = side_near_far_signature_statistics(
        scores, distance,
        side_names={"negative": "CI", "positive": "Cd"},
        min_group_size=20,
    )
    assert set(table.side_name) == {"CI", "Cd"}
    assert table.near_minus_far.gt(0).all()


def test_provenance_rejects_circular_signature_validation():
    obj = ad.AnnData(np.ones((4, 2)))
    with pytest.raises(ValueError, match="overlap"):
        record_st_validation_provenance(
            obj, annotation_signatures=["tumor_core"],
            validation_signatures=["tumor_core", "hypoxia"],
        )
    payload = record_st_validation_provenance(
        obj, annotation_signatures=["tumor_core"],
        validation_signatures=["hypoxia", "invasion"],
    )
    assert payload["independent_validation"] is True
