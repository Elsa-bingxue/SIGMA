"""Generate the complete downstream HBC525 regression report."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from generate_hbc515_complete_report import GENE_SETS
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
RESULT = ROOT / "HBC_525/clean/results/HBC525_SIGMA_core.h5ad"
ASSIGNMENTS = ROOT / "HBC_525/clean/results/BC525_metabolic_program_assignments_reference.csv"
OUTPUT = ROOT / "HBC_525/clean/results/pypi_validation_full_report"
EXPECTED_BOUNDARY_PROGRAM = 2  # Frozen BC525 biology-figure interpretation.


def main():
    set_publication_style()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    result = ad.read_h5ad(RESULT)
    assignments = pd.read_csv(ASSIGNMENTS)
    report = program_report_data(result, assignments)

    signature_scores = score_st_signatures(result, GENE_SETS, random_state=0)
    validation, selected_program = stroma_near_far_statistics(
        result, report["scores"], signature_scores=signature_scores,
        spot_spacing=21.0, near_width=5, far_start=8,
    )
    validation.to_csv(OUTPUT / "HBC525_stroma_side_near_far_statistics.csv", index=False)
    plot_stroma_near_far_validation(
        validation, OUTPUT, leading_program=selected_program, prefix="HBC525"
    )
    if selected_program != EXPECTED_BOUNDARY_PROGRAM:
        raise AssertionError(
            f"Expected frozen HBC525 interface program {EXPECTED_BOUNDARY_PROGRAM}, "
            f"got {selected_program}"
        )

    representatives = plot_program_report(
        result, report, OUTPUT, boundary_cluster=selected_program, prefix="HBC525"
    )
    representatives.to_csv(OUTPUT / "HBC525_boundary_program_representatives.csv", index=False)

    # Preserve the leading ranked BC525 feature as the anisotropy example.
    target = assignments.sort_values("r2_logI", ascending=False).iloc[0]
    anisotropy = anisotropy_report(result, int(target["j"]), n_sectors=8, min_points=200)
    anisotropy["table"].to_csv(OUTPUT / "HBC525_representative_sector_lambda.csv", index=False)
    plot_anisotropy_polar(anisotropy, OUTPUT, prefix="HBC525_representative")

    summary = pd.DataFrame([{
        "dataset": "HBC525", "selected_program": selected_program,
        "expected_program": EXPECTED_BOUNDARY_PROGRAM,
        "representative_mz": float(target["mz"]),
        "anisotropy_cv": anisotropy["anisotropy_cv"],
        "n_signatures": len(signature_scores),
    }])
    summary.to_csv(OUTPUT / "HBC525_complete_report_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Saved complete report to {OUTPUT}")


if __name__ == "__main__":
    main()
