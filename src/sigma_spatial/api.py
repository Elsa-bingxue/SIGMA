"""High-level public entry points for SIGMA."""

from __future__ import annotations

from .pipeline import SIGMA
from .validation import validate_sigma_input, validate_sigma_output


def run_sigma(
    adata,
    *,
    anchor_key="sigma_anchor",
    representation_key="X_harmony",
    spatial_key="spatial",
    copy=True,
    n_components=64,
    n_neighbors=15,
    random_state=0,
    device=None,
    **fit_kwargs,
):
    """Run the validated end-to-end SIGMA workflow on an AnnData object.

    Anchors must be encoded as 1 (tumor), 0 (non-tumor), and NaN (unknown).
    The current validated model also requires an auxiliary representation such
    as matched RNA PCs in ``adata.obsm[representation_key]``.
    """
    validate_sigma_input(
        adata, anchor_key=anchor_key, representation_key=representation_key, spatial_key=spatial_key
    )
    result = adata.copy() if copy else adata
    if spatial_key != "spatial":
        result.obsm["spatial"] = result.obsm[spatial_key].copy()
    temporary_key = "_sigma_public_anchor"
    anchor = result.obs[anchor_key].to_numpy()
    labels = ["Tumor" if value == 1 else "Stroma" if value == 0 else None for value in anchor]
    result.obs[temporary_key] = labels
    estimator = SIGMA(n_components=n_components, k=n_neighbors, seed=random_state, device=device)
    estimator.fit(
        result,
        annotation_key=temporary_key,
        tumor_label="Tumor",
        stroma_label="Stroma",
        rna_key=representation_key,
        **fit_kwargs,
    )
    del result.obs[temporary_key]
    validate_sigma_output(result)
    result.uns.setdefault("sigma", {})["public_api"] = {
        "anchor_key": anchor_key,
        "representation_key": representation_key,
        "spatial_key": spatial_key,
    }
    return result
