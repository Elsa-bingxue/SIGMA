import scanpy as sc
from sigma_spatial import run_sigma

adata = sc.read_h5ad("BC_515_Section_1.h5ad")
adata.obs["sigma_anchor"] = adata.obs["annotation"].map({"Tumor": 1, "Stroma": 0})
result = run_sigma(
    adata,
    anchor_key="sigma_anchor",
    representation_key="X_harmony",
    random_state=0,
    epochs=1000,
)
result.write_h5ad("BC_515_SIGMA.h5ad")
