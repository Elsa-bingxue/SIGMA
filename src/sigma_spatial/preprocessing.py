import numpy as np
from scipy.sparse import issparse
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import RobustScaler


def get_msi_matrix(adata, layer=None):
    """Return an n_spots x n_features MSI matrix from AnnData."""
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer {layer!r} not found in adata.layers")
        return adata.layers[layer]
    if "msi" in adata.uns:
        x = adata.uns["msi"]
        if hasattr(x, "shape") and x.shape[0] == adata.n_obs:
            return x
    if "raw" in adata.layers:
        x = adata.layers["raw"]
        if hasattr(x, "shape") and x.shape[0] == adata.n_obs:
            return x
    return adata.X


def msi_to_embedding(x, n_components=64, random_state=0):
    """log1p + robust scaling (dense only) + TruncatedSVD + z-score."""
    max_components = min(x.shape) - 1
    k = min(int(n_components), max_components)
    if k < 1:
        raise ValueError("MSI matrix is too small for SVD.")
    if issparse(x):
        xlog = x.copy()
        xlog.data = np.log1p(xlog.data)
        e = TruncatedSVD(n_components=k, random_state=random_state).fit_transform(xlog)
    else:
        xlog = np.log1p(np.asarray(x, dtype=np.float32))
        xs = RobustScaler(quantile_range=(10, 90)).fit_transform(xlog)
        e = TruncatedSVD(n_components=k, random_state=random_state).fit_transform(xs)
    e = np.asarray(e, dtype=np.float32)
    return (e - e.mean(0)) / (e.std(0) + 1e-6)
