"""SIGMA: spatial interface analysis for spatial metabolomics and transcriptomics."""
from .preprocessing import get_msi_matrix, msi_to_embedding
from .boundary import (
    build_boundary_from_binary_mask,
    build_boundary_from_field,
    signed_distance_from_boundary_points,
)
from .io import load_spatial_metabolomics
from .validation import validate_sigma_input, validate_sigma_output

__version__ = "0.1.1"
__all__ = [
    "SIGMA", "run_sigma", "run_sigma_weak_anchor", "load_spatial_metabolomics",
    "validate_sigma_input", "validate_sigma_output",
    "get_msi_matrix", "msi_to_embedding",
    "build_boundary_from_binary_mask", "build_boundary_from_field",
    "signed_distance_from_boundary_points",
]


def __getattr__(name):
    """Keep light-weight utilities importable without the optional torch stack."""
    if name == "SIGMA":
        from .pipeline import SIGMA
        return SIGMA
    if name == "run_sigma":
        from .api import run_sigma
        return run_sigma
    if name == "run_sigma_weak_anchor":
        from .weak_anchor import run_sigma_weak_anchor
        return run_sigma_weak_anchor
    raise AttributeError(name)
