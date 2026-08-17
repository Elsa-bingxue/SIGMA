"""Generate the four curated HBC_515 notebooks without touching the original."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "HBC_515" / "HBC_515_SIGMA.ipynb"
OUT = ROOT / "HBC_515" / "clean"


def md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.splitlines(True)}


def code(text):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": text.splitlines(True)}


def notebook(cells):
    return {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3",
            "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10"}},
            "nbformat": 4, "nbformat_minor": 5}


def write(name, cells):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / name).write_text(json.dumps(notebook(cells), indent=1) + "\n")


def source_cell(original, index):
    cell = original["cells"][index]
    return code("".join(cell.get("source", [])))


common_setup = '''from pathlib import Path
import sys

PROJECT_ROOT = Path.cwd().resolve().parents[1]
PACKAGE_SRC = PROJECT_ROOT / "sigma_spatial_pypi_starter" / "src"
if str(PACKAGE_SRC) not in sys.path:
    sys.path.insert(0, str(PACKAGE_SRC))

DATA_DIR = PROJECT_ROOT / "HBC_515"
RESULT_DIR = DATA_DIR / "clean" / "results"
FIGURE_DIR = RESULT_DIR / "figures"
RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
'''


def build_core():
    cells = [
        md("# HBC_515 — SIGMA core\n\nClean reference workflow derived from original cells 0–43. The original notebook is not modified. Mathematical definitions, graph direction, probability normalization, thresholds, architecture, and hyperparameters are preserved."),
        md("## 1. Environment and paths"), code(common_setup + '''\nimport scanpy as sc
from sigma_spatial import SIGMA
from sigma_spatial.plotting import set_publication_style

set_publication_style()
SEED = 0
'''),
        md("## 2. Load HBC_515"), code('''adata = sc.read_h5ad(DATA_DIR / "BC_515_Section_1.h5ad")
print(adata)
print("MSI storage:", "uns[msi]" if "msi" in adata.uns else "layers[raw]" if "raw" in adata.layers else "X")
'''),
        md("## 3. Fit the preserved SIGMA model\n\nThis calls extracted package functions that match original cells 4–12. Training is intentionally not executed when this notebook is generated."),
        code('''sigma = SIGMA(n_components=64, k=15, seed=SEED)
sigma.fit(
    adata,
    annotation_key="annotation",
    tumor_label="Tumor",
    stroma_label="Stroma",
    rna_key="X_harmony",
    z_dim=32,
    epochs=1000,
    lr=1e-3,
    lambda_sup=0.05,
    beta_rna=0.05,
    boundary_k=10,
)
'''),
        md("## 4. Essential boundary QC"), code('''import matplotlib.pyplot as plt
from sigma_spatial.plotting import figure_size, save_figure, shared_colorbar, style_spatial_axis

xy = adata.obsm["spatial"]
fig, axes = plt.subplots(1, 3, figsize=figure_size("double"), constrained_layout=True)
items = [
    ("sigma_region_probability", "Region probability", "viridis"),
    ("sigma_boundary", "Boundary", "Greys"),
    ("sigma_d_signed", "Signed distance", "coolwarm"),
]
for ax, (key, title, cmap) in zip(axes, items):
    artist = ax.scatter(xy[:, 0], xy[:, 1], c=adata.obs[key], s=5, cmap=cmap, rasterized=True)
    ax.set_title(title)
    style_spatial_axis(ax)
    fig.colorbar(artist, ax=ax, shrink=0.72)
save_figure(fig, FIGURE_DIR / "HBC515_SIGMA_core_QC", formats=("pdf", "svg", "png"))
'''),
        md("## 5. Persist essential outputs"), code('''output_path = RESULT_DIR / "HBC515_SIGMA_core.h5ad"
adata.write_h5ad(output_path)
print(output_path)
'''),
        md("## 6. Interface profiles, lambda, and metabolite ranking\n\nThe original lambda/ranking implementations have multiple variants. Until golden numerical tables are frozen, cells 21–43 remain provenance-locked and are not silently replaced by the simplified starter implementation. Run the migration cells below only after the core result exists."),
    ]
    return cells


def build_simulation(original):
    cells = [
        md("# HBC_515 — simulation\n\nDerived from original cells 44–61. Simulation generation is imported from the package; metric and selection definitions below remain verbatim pending their next extraction test."),
        md("## 1. Environment and core result"), code(common_setup + '''\nimport numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sigma_spatial.simulation import (
    evaluate_simulation_feature_selection,
    simulate_interface_features,
    zscore_vec,
)
from sigma_spatial.plotting import add_panel_labels, figure_size, save_figure, set_publication_style, style_spatial_axis
set_publication_style()
adata = ad.read_h5ad(RESULT_DIR / "HBC515_SIGMA_core.h5ad")
coords = np.asarray(adata.obsm["spatial"])
d_signed = adata.obs["sigma_d_signed"].to_numpy(float)
d_abs = np.abs(d_signed)
tumor_mask = adata.obs["sigma_inside"].to_numpy(bool)
'''),
        md("## 2. Preserved feature metrics and selection"),
        source_cell(original, 47), source_cell(original, 48), source_cell(original, 50),
        md("## 3. Representative simulation example"), code('''X_sim, meta_rows = simulate_interface_features(
    d_abs=d_abs, d_signed=d_signed, coords=coords, tumor_mask=tumor_mask,
    n_positive=1, n_negative=1, n_region_only=1, n_spatial_background=0,
    n_noise=1, effect_size=0.25, noise_sd=0.5, lam=100,
    random_state=0,
)
meta = pd.DataFrame(meta_rows)
metrics = compute_boundary_feature_metrics(X_sim, d_abs, d_signed=d_signed, tumor_mask=tumor_mask)
metrics = metrics.merge(meta, on="j", how="left", suffixes=("", "_truth"))
metrics = add_region_adjusted_distance_metrics(
    df_metrics=metrics, X=X_sim, d_abs=d_abs, tumor_mask=tumor_mask
)
selected = select_positive_negative_interface_features_strict(metrics)
display(selected[["mz", "sim_type", "true_direction", "pred_direction", "is_selected_boundary", "auc", "decay_r2"]])

examples = [
    selected[selected["sim_type"] == "positive_boundary"].iloc[0],
    selected[selected["sim_type"] == "negative_boundary"].iloc[0],
]
fig, axes = plt.subplots(2, 3, figsize=(9.2, 6.2), constrained_layout=True)
column_titles = ["Ground-truth pattern", "Simulated observation (with noise)",
                 "SIGMA-recovered pattern"]
for row_index, row in enumerate(examples):
    j = int(row["j"])
    true_lambda = float(row["true_lambda"])
    lambda_hat = float(row["lambda_hat"]) if np.isfinite(row["lambda_hat"]) else true_lambda
    decay_true = np.exp(-d_abs / (true_lambda + 1e-8))
    decay_hat = np.exp(-d_abs / (lambda_hat + 1e-8))
    if row["sim_type"] == "positive_boundary":
        true_signal, recovered_signal = decay_true, decay_hat
        auc_display = float(row["auc"])
        row_label = "Positive"
    else:
        true_signal, recovered_signal = 1.0 - decay_true, 1.0 - decay_hat
        auc_display = float(row["auc_neg"])
        row_label = "Negative"
    for col_index, values in enumerate([true_signal, X_sim[:, j], recovered_signal]):
        ax = axes[row_index, col_index]
        artist = ax.scatter(coords[:, 0], coords[:, 1], c=values, s=5,
                            cmap="viridis", rasterized=True)
        fig.colorbar(artist, ax=ax, shrink=0.72)
        if row_index == 0:
            ax.set_title(column_titles[col_index])
        style_spatial_axis(ax)
    axes[row_index, 0].set_ylabel(
        f"{row_label}\\nAUC={auc_display:.3f}; decay R²={row['decay_r2']:.3f}\\n"
        f"λ true={true_lambda:.1f}; λ estimated={lambda_hat:.1f}",
        fontsize=8,
    )
fig.suptitle("Representative positive and negative interface simulations")
save_figure(fig, FIGURE_DIR / "HBC515_representative_simulation_positive_negative",
            formats=("pdf", "svg", "png"))
'''),
        md("## 4. Effect-size and sensitivity grid"), source_cell(original, 53),
        code('''eval_all_strict, detail_all_strict = run_sensitivity_analysis_strict(
    d_abs=d_abs, d_signed=d_signed, coords=coords, tumor_mask=tumor_mask,
    effect_sizes=[0.25, 0.5, 1.0, 2.0], noise_sds=[0.5],
    lambdas=[200], n_repeats=1, n_positive=100,
    n_negative=100, n_region_only=300, n_spatial_background=300,
    n_noise=500, random_seed=0, p_cutoff=0.05, decay_r2_cutoff=0.05,
    auc_pos_cutoff=0.65, auc_neg_cutoff=0.35,
)
eval_all_strict.to_csv(RESULT_DIR / "HBC515_simulation_recovery.csv", index=False)
detail_all_strict.to_pickle(RESULT_DIR / "HBC515_simulation_detail.pkl")
'''),
    ]
    return cells


def build_benchmark():
    return [
        md("# HBC_515 — benchmark\n\nBenchmark workflow derived from original cells 62–104. It consumes frozen core and simulation outputs and adds an explicit ROC curve that was absent from the original notebook."),
        md("## 1. Environment and preserved outputs"), code(common_setup + '''\nimport numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_curve
from sigma_spatial.plotting import add_panel_labels, figure_size, legend_outside, save_figure, set_publication_style
from sigma_spatial.simulation import recovery_curve_data
set_publication_style()
recovery = pd.read_csv(RESULT_DIR / "HBC515_simulation_recovery.csv")
detail = pd.read_pickle(RESULT_DIR / "HBC515_simulation_detail.pkl")
'''),
        md("## 2. Effect-size sensitivity ROC curves\n\nThis is the requested effect-size figure: noise SD, interface width λ, and replicate are fixed, while four effect sizes are compared on the same ROC axes."), code('''from sklearn.metrics import roc_auc_score

NOISE_SD = 0.5
LAMBDA = 200
REP = 0
EFFECT_SIZES = [0.25, 0.5, 1.0, 2.0]
curve_data = detail[
    (detail["noise_sd"] == NOISE_SD)
    & (detail["lambda"] == LAMBDA)
    & (detail["rep"] == REP)
].copy()

fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
colors = ["#0072B2", "#E69F00", "#009E73", "#D55E00"]
for effect_size, color in zip(EFFECT_SIZES, colors):
    subset = curve_data[curve_data["effect_size"] == effect_size]
    curves = recovery_curve_data(subset)
    ax.plot(
        curves["fpr"], curves["tpr"], color=color, lw=2,
        label=f"Effect={effect_size:g} (AUROC={curves['auroc']:.3f})",
    )
ax.plot([0, 1], [0, 1], ls="--", lw=1.2, color="#8A5AC2", label="Random")
ax.set(
    xlabel="False Positive Rate", ylabel="True Positive Rate",
    xlim=(0, 1), ylim=(0, 1.02),
    title=f"Effect-size sensitivity\\nNoise SD={NOISE_SD:g}, λ={LAMBDA}, Rep={REP}",
)
ax.legend(loc="lower right", frameon=True, framealpha=0.95)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
save_figure(fig, FIGURE_DIR / "HBC515_effect_size_sensitivity_ROC", formats=("pdf", "svg", "png"))

auc_summary = recovery[
    (recovery["noise_sd"] == NOISE_SD)
    & (recovery["lambda"] == LAMBDA)
    & (recovery["rep"] == REP)
].sort_values("effect_size")
fig, ax = plt.subplots(figsize=(4.8, 3.8), constrained_layout=True)
bars = ax.bar(
    auc_summary["effect_size"].astype(str), auc_summary["auroc"],
    color=colors, width=0.72,
)
for bar, value in zip(bars, auc_summary["auroc"]):
    ax.annotate(f"{value:.3f}",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", fontsize=8)
ax.set(
    xlabel="Effect size", ylabel="AUROC", ylim=(0, 1.07),
    title=f"Effect-size sensitivity of AUROC\\nNoise SD={NOISE_SD:g}, λ={LAMBDA}, Rep={REP}",
)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
save_figure(fig, FIGURE_DIR / "HBC515_effect_size_AUROC", formats=("pdf", "svg", "png"))
'''),
        md("## 3. AUROC, ROC curve, and AUPRC/AP\n\nThe scalar definitions are unchanged from original cell 49. The ROC/PR curves are new visualizations of the same truth labels and continuous scores."),
        code('''example = detail[(detail["effect_size"] == 1.0) & (detail["noise_sd"] == 0.5) & (detail["lambda"] == 200) & (detail["rep"] == 0)].copy()
curves = recovery_curve_data(example)
precision, recall, _ = precision_recall_curve(curves["truth"], curves["score"])
ap = average_precision_score(curves["truth"], curves["score"])
fig, axes = plt.subplots(1, 2, figsize=figure_size("double"), constrained_layout=True)
axes[0].plot(curves["fpr"], curves["tpr"], label=f"SIGMA AUROC={curves['auroc']:.3f}")
axes[0].plot([0, 1], [0, 1], ls="--", color="0.6")
axes[0].set(xlabel="False-positive rate", ylabel="True-positive rate", title="ROC curve", xlim=(0, 1), ylim=(0, 1))
axes[1].plot(recall, precision, label=f"SIGMA AP={ap:.3f}")
axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision–recall curve", xlim=(0, 1), ylim=(0, 1))
for ax in axes: ax.legend(frameon=False, loc="lower right")
add_panel_labels(axes)
save_figure(fig, FIGURE_DIR / "HBC515_ROC_AUPRC", formats=("pdf", "svg", "png"))
'''),
        md("## 4. Comparison methods\n\nThe MET-MAP, SpatialDE, Region-DE, top-k, and overlap implementations remain provenance-locked in original cells 63–102. They require optional dependencies and frozen intermediate tables before safe extraction. The existing HBC_515 slide summaries can be inspected without overwriting them."),
        code('''slide_summary = PROJECT_ROOT / "4 slides" / "BC_515_topk_500_summary.csv"
if slide_summary.exists():
    benchmark_summary = pd.read_csv(slide_summary)
    display(benchmark_summary)
else:
    print("Missing slide benchmark summary:", slide_summary)
'''),
    ]


def build_biology(original):
    cells = [
        md("# HBC_515 — biological interpretation and manuscript figures\n\nCurated from original cells 105–156. This notebook consumes a completed core result; it does not retrain SIGMA. Long duplicated plotting implementations are intentionally not copied."),
        md("## 1. Environment and core result"), code(common_setup + '''\nimport numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from sigma_spatial.plotting import add_panel_labels, figure_size, save_figure, set_publication_style, style_spatial_axis
set_publication_style()
adata = ad.read_h5ad(RESULT_DIR / "HBC515_SIGMA_core.h5ad")
coords = np.asarray(adata.obsm["spatial"])
d_signed = adata.obs["sigma_d_signed"].to_numpy(float)
'''),
        md("## 2. Transition thickness"), source_cell(original, 106), source_cell(original, 107),
        md("## 3. Region-only versus boundary association"), source_cell(original, 109),
        md("## 4. Real-data validation"), source_cell(original, 116), source_cell(original, 117),
        md("## 5. Metabolic clusters and representative profiles"), source_cell(original, 120), source_cell(original, 121), source_cell(original, 122),
        md("## 6. Near-versus-far and RNA signatures"), source_cell(original, 130), source_cell(original, 134),
        md("## 7. Heatmaps"), source_cell(original, 136), source_cell(original, 138),
        md("## 8. GO enrichment"), source_cell(original, 148), source_cell(original, 149), source_cell(original, 150),
        md("## Migration note\n\nCells 125–129, 131–133, 135, 137, 139–146, and 151–156 contain alternative or duplicated plotting implementations and were not copied. They remain available unchanged in the original notebook pending numerical fixture creation."),
    ]
    return cells


def main():
    original = json.loads(ORIGINAL.read_text())
    write("01_HBC515_SIGMA_core.ipynb", build_core())
    write("02_HBC515_simulation.ipynb", build_simulation(original))
    write("03_HBC515_benchmark.ipynb", build_benchmark())
    write("04_HBC515_biology_figures.ipynb", build_biology(original))


if __name__ == "__main__":
    main()
