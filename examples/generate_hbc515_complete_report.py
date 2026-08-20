"""Generate all downstream HBC515 validation figures without refitting SIGMA."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from sigma_spatial import (
    anisotropy_report,
    plot_anisotropy_polar,
    plot_program_report,
    plot_stroma_near_far_validation,
    program_report_data,
    score_st_signatures,
    stroma_near_far_statistics,
)
from sigma_spatial.plotting import set_publication_style


ROOT = Path(__file__).resolve().parents[2]
RESULT = ROOT / "HBC_515/clean/results/pypi_validation_0.1.1/HBC515_sigma_omics_0.1.1.h5ad"
ASSIGNMENTS = ROOT / "HBC_515/clean/results/BC515_metabolic_program_assignments.csv"
OUTPUT = ROOT / "HBC_515/clean/results/pypi_validation_0.1.1/figures/full_report"
TARGET_MZ = 817.5345203639481
BOUNDARY_PROGRAM = 3  # Frozen biological interpretation; no re-selection here.
GENE_SETS = {
    "epithelial_score": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CD24", "CLDN7", "PRSS8", "ELF3"],
    "luminal_secretory_score": ["MUC1", "AGR2", "FOXA1", "ELF3", "XBP1", "GOLM1", "SCGB2A2", "SLC44A4"],
    "stress_adaptation_score": ["ATF3", "JUNB", "DUSP4", "GADD45G", "DDIT3", "HSPA1A", "HSPA1B", "HMOX1", "SOD2", "NQO1", "UCP2", "HMGCS2", "ALDH3B2", "ALOX15B"],
    "plasma_cell_score": ["IGKC", "IGHG1", "IGHG3", "IGHM", "JCHAIN", "MZB1", "TXNDC5", "XBP1"],
    "immune_leukocyte_score": ["PTPRC", "LST1", "LYZ", "TYROBP", "CD3D", "CD3E", "MS4A1", "CD79A"],
    "ecm_matrix_score": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "POSTN", "SPARC", "MMP2", "TIMP1"],
    "caf_fibroblast_score": ["ACTA2", "TAGLN", "FAP", "THY1", "PDGFRA", "PDGFRB", "CXCL12", "COL1A1", "COL1A2", "FN1"],
    "emt_migration_score": ["VIM", "FN1", "ITGA5", "ITGB1", "MMP2", "MMP14", "SERPINE1", "TGFBI", "SPARC", "SNAI2"],
    "proliferation_score": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "MCM6", "UBE2C", "BIRC5"],
}


def main():
    set_publication_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    result = ad.read_h5ad(RESULT)
    assignments = pd.read_csv(ASSIGNMENTS)

    report = program_report_data(result, assignments)
    representatives = plot_program_report(
        result, report, OUTPUT, boundary_cluster=BOUNDARY_PROGRAM, prefix="HBC515"
    )
    representatives.to_csv(OUTPUT / "HBC515_boundary_program_representatives.csv", index=False)

    mz_names = np.asarray(result.uns["mz_features"], float)
    feature_index = int(np.argmin(np.abs(mz_names - TARGET_MZ)))
    anisotropy = anisotropy_report(result, feature_index)
    anisotropy["table"].to_csv(OUTPUT / "HBC515_mz817_5345_sector_lambda.csv", index=False)
    plot_anisotropy_polar(anisotropy, OUTPUT, prefix="HBC515_mz817_5345")

    signature_scores = score_st_signatures(result, GENE_SETS, random_state=0)
    validation, selected_program = stroma_near_far_statistics(
        result, report["scores"], signature_scores=signature_scores,
        spot_spacing=21.0, near_width=5, far_start=8,
    )
    validation.to_csv(OUTPUT / "HBC515_stroma_side_near_far_statistics.csv", index=False)
    plot_stroma_near_far_validation(
        validation, OUTPUT, leading_program=selected_program, prefix="HBC515"
    )
    if selected_program != BOUNDARY_PROGRAM:
        raise AssertionError(
            f"Expected frozen HBC515 boundary program {BOUNDARY_PROGRAM}, got {selected_program}"
        )
    print(f"Saved complete report to {OUTPUT}")
    print(f"Anisotropy CV: {anisotropy['anisotropy_cv']:.6f}")
    print(f"Automatically selected boundary-associated program: {selected_program}")


if __name__ == "__main__":
    main()
