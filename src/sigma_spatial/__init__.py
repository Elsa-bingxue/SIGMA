"""SIGMA: spatial interface analysis for spatial metabolomics and transcriptomics."""
from .preprocessing import get_msi_matrix, msi_to_embedding
from .boundary import (
    build_boundary_from_binary_mask,
    build_boundary_from_field,
    signed_distance_from_boundary_points,
)

__version__ = "0.1.0a1"
__all__ = [
    "SIGMA", "get_msi_matrix", "msi_to_embedding",
    "build_boundary_from_binary_mask", "build_boundary_from_field",
    "signed_distance_from_boundary_points",
]


def __getattr__(name):
    """Keep light-weight utilities importable without the optional torch stack."""
    if name == "SIGMA":
        from .pipeline import SIGMA
        return SIGMA
    raise AttributeError(name)
