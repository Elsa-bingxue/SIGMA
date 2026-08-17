"""SM-only weak-anchor SIGMA preserved from the HCC P1/P4 notebooks."""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.sparse import issparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import RobustScaler
from torch_geometric.data import Data

from .boundary import build_boundary_from_field, signed_distance_from_boundary_points
from .graph import gaussian_knn_graph, gaussian_label_smoothing, smooth_embedding
from .model import SpectralHighPassBlock
from .preprocessing import get_msi_matrix
from .utils import set_seed
from .validation import validate_sigma_output


def weak_anchor_msi_embedding(matrix, n_components=64, random_state=0, already_log=True):
    """HCC-reference MSI embedding, including its ``already_log`` behavior."""
    limit = min(matrix.shape[0] - 1, matrix.shape[1] - 1)
    if limit < 2:
        raise ValueError("Weak-anchor MSI input requires at least three spots and three features.")
    k = min(int(n_components), limit)
    if issparse(matrix):
        values = matrix.copy()
        if not already_log:
            values.data = np.log1p(values.data)
    else:
        values = np.asarray(matrix, dtype=np.float32)
        if not already_log:
            values = np.log1p(values)
        values = RobustScaler(
            with_centering=True, with_scaling=True, quantile_range=(10, 90)
        ).fit_transform(values).astype(np.float32)
    embedding = TruncatedSVD(n_components=k, random_state=random_state).fit_transform(values).astype(np.float32)
    return (embedding - embedding.mean(0)) / (embedding.std(0) + 1e-6)


class _WeakAnchorSpectralResidualNet(nn.Module):
    """Exact layer construction order used by both HCC reference notebooks."""

    def __init__(self, in_dim, hidden=128, out_dim=64, z_dim=32):
        super().__init__()
        self.hp1 = SpectralHighPassBlock(in_dim, hidden, dropout=0.1)
        self.hp2 = SpectralHighPassBlock(hidden, out_dim, dropout=0.1)
        self.proj = nn.Linear(out_dim, z_dim)
        self.rec_head = nn.Linear(out_dim, out_dim)
        self.cls = nn.Linear(out_dim, 1)

    def forward(self, x, edge_index, edge_weight=None):
        h = self.hp2(self.hp1(x, edge_index, edge_weight), edge_index, edge_weight)
        return {
            "r": self.rec_head(h),
            "z_hat": self.proj(h),
            "h": h,
            "logit": self.cls(h).squeeze(1),
        }


def _validate_weak_anchor_input(adata, anchor_key, spatial_key):
    if spatial_key not in adata.obsm:
        raise KeyError(f"adata.obsm[{spatial_key!r}] is required")
    xy = np.asarray(adata.obsm[spatial_key])
    if xy.ndim != 2 or xy.shape[0] != adata.n_obs or xy.shape[1] < 2 or not np.all(np.isfinite(xy[:, :2])):
        raise ValueError("Spatial coordinates must be finite with shape (n_obs, >=2).")
    if anchor_key not in adata.obs:
        raise KeyError(f"adata.obs[{anchor_key!r}] is required")
    anchors = adata.obs[anchor_key].to_numpy()
    finite = ~adata.obs[anchor_key].isna().to_numpy()
    values = set(anchors[finite].tolist())
    if not {0, 1} <= values:
        raise ValueError("Weak anchors must contain 1 (tumor) and 0 (non-tumor); use NaN for unknown spots.")


def run_sigma_weak_anchor(
    adata,
    *,
    anchor_key="sigma_anchor",
    spatial_key="spatial",
    copy=True,
    already_log=True,
    n_components=64,
    n_neighbors=15,
    random_state=0,
    device=None,
    epochs=1000,
    learning_rate=1e-3,
    lambda_supervised=0.05,
    boundary_level=0.5,
    boundary_neighbors=10,
    verbose=False,
):
    """Run the frozen HCC SM-only weak-anchor workflow.

    This function deliberately has no RNA loss. Defaults reproduce the executed
    P1/P4 notebook branch and must be changed explicitly by the caller.
    """
    _validate_weak_anchor_input(adata, anchor_key, spatial_key)
    set_seed(random_state)
    result = adata.copy() if copy else adata
    xy = np.asarray(result.obsm[spatial_key], dtype=np.float32)
    matrix = get_msi_matrix(result)
    embedding = weak_anchor_msi_embedding(
        matrix, n_components=n_components, random_state=random_state, already_log=already_log
    )
    edge_index, edge_weight, sigma = gaussian_knn_graph(xy, k=n_neighbors, symmetrize=False)
    anchors = result.obs[anchor_key].to_numpy()
    labeled = ~result.obs[anchor_key].isna().to_numpy()
    y = np.full(result.n_obs, -1, dtype=np.int64)
    y[labeled] = anchors[labeled].astype(np.int64)
    gaussian_anchor = gaussian_label_smoothing(
        edge_index, edge_weight, y, labeled, n_iter=50, alpha=0.85
    )
    gaussian_embedding = smooth_embedding(edge_index, edge_weight, embedding, n_iter=10, alpha=0.9)
    residual_target = (embedding - gaussian_embedding).astype(np.float32)

    target_device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
    data = Data(
        x=torch.tensor(gaussian_embedding, dtype=torch.float32),
        edge_index=torch.tensor(edge_index, dtype=torch.long),
        edge_attr=torch.tensor(edge_weight, dtype=torch.float32),
    ).to(target_device)
    residual_tensor = torch.tensor(residual_target, dtype=torch.float32, device=target_device)
    labeled_index = torch.tensor(np.where(labeled)[0], dtype=torch.long, device=target_device)
    label_tensor = torch.tensor(y[labeled].astype(np.float32), dtype=torch.float32, device=target_device)
    model = _WeakAnchorSpectralResidualNet(
        embedding.shape[1], hidden=128, out_dim=embedding.shape[1], z_dim=32
    ).to(target_device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=1e-4)
    positive = int((y[labeled] == 1).sum()); negative = int((y[labeled] == 0).sum())
    positive_weight = torch.tensor([negative / max(positive, 1)], device=target_device)
    loss_history = []
    for epoch in range(1, int(epochs) + 1):
        model.train(); optimizer.zero_grad()
        output = model(data.x, data.edge_index, edge_weight=data.edge_attr)
        reconstruction_loss = F.smooth_l1_loss(output["r"], residual_tensor)
        supervised_loss = F.binary_cross_entropy_with_logits(
            output["logit"][labeled_index], label_tensor, pos_weight=positive_weight
        )
        loss = reconstruction_loss + float(lambda_supervised) * supervised_loss
        loss.backward(); optimizer.step(); loss_history.append(float(loss.detach().cpu()))
        if verbose and (epoch == 1 or epoch % 100 == 0 or epoch == int(epochs)):
            print(
                f"epoch={epoch:04d} loss={loss_history[-1]:.6f} "
                f"rec={float(reconstruction_loss.detach().cpu()):.6f} "
                f"sup={float(supervised_loss.detach().cpu()):.6f}",
                flush=True,
            )

    model.eval()
    with torch.no_grad():
        output = model(data.x, data.edge_index, edge_weight=data.edge_attr)
        residual = output["r"].cpu().numpy()
        raw_probability = torch.sigmoid(output["logit"]).cpu().numpy()
    probability = (raw_probability - np.nanmin(raw_probability)) / (
        np.nanmax(raw_probability) - np.nanmin(raw_probability) + 1e-8
    )
    corrected = gaussian_embedding + residual
    hp_score = np.linalg.norm(residual, axis=1)
    hp_score = (hp_score - hp_score.min()) / (hp_score.max() - hp_score.min() + 1e-8)
    boundary, inside = build_boundary_from_field(
        xy, probability, level=boundary_level, k_nn=boundary_neighbors
    )
    if boundary.sum() < 2:
        raise ValueError("Boundary construction produced fewer than two boundary spots.")
    signed_distance = signed_distance_from_boundary_points(xy, boundary, inside)

    result.obsm["X_sigma_msi"] = embedding
    result.obsm["X_sigma_gauss"] = gaussian_embedding
    result.obsm["X_sigma_residual"] = residual
    result.obsm["X_sigma_corrected"] = corrected
    result.obs["sigma_gaussian_anchor"] = gaussian_anchor
    result.obs["sigma_region_probability_raw"] = raw_probability
    result.obs["sigma_region_probability"] = probability
    result.obs["sigma_inside"] = inside
    result.obs["sigma_boundary"] = boundary
    result.obs["sigma_d_signed"] = signed_distance
    # Reference HCC aliases retained for direct comparison with frozen notebooks.
    result.obsm["X_msi_embed"] = embedding
    result.obsm["E_gauss"] = gaussian_embedding
    result.obsm["R_spectral"] = residual
    result.obsm["E_hat_gauss_spectral"] = corrected
    result.obs["f_gauss"] = gaussian_anchor
    result.obs["hp_score_spectral"] = hp_score
    result.obs["p_tumor_gauss_spectral_raw"] = raw_probability
    result.obs["p_tumor_gauss_spectral"] = probability
    result.obs["is_boundary"] = boundary.astype(int)
    result.obs["inside_tumor"] = inside.astype(int)
    result.obs["d_signed"] = signed_distance
    result.obs["d_abs"] = np.abs(signed_distance)
    result.uns["sigma"] = {
        "mode": "sm_only_weak_anchor_hcc_reference",
        "sigma": sigma,
        "k": int(n_neighbors),
        "seed": int(random_state),
        "already_log": bool(already_log),
        "epochs": int(epochs),
        "lambda_supervised": float(lambda_supervised),
        "rna_loss": False,
        "loss_history": loss_history,
    }
    validate_sigma_output(result)
    return result
