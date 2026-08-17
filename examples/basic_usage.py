import scanpy as sc
from sigma_spatial import SIGMA

adata = sc.read_h5ad("BC_515_Section_1.h5ad")
model = SIGMA(n_components=64, k=15, seed=0)
model.fit(adata, annotation_key="annotation", rna_key="X_harmony", epochs=1000)
adata.write_h5ad("BC_515_SIGMA.h5ad")
