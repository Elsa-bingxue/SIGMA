# Reference data

SIGMA does not redistribute the reference datasets through PyPI. Obtain data
from the cited study or repository, comply with the source licence and convert
the required matrices to the documented [AnnData input format](input_format.md).

## Dataset provenance

| SIGMA samples | Disease | Data used in this study | Annotation setting | Source |
|---|---|---|---|---|
| `BC_515`, `BC_525` | Breast cancer | Matched spatial metabolomics and transcriptomics from the same tissue, with pathological annotations | Direct pathological anchors | Godfrey *et al.*, *Angewandte Chemie International Edition* (2025), [doi:10.1002/anie.202502028](https://doi.org/10.1002/anie.202502028) |
| `LC_091`, `LC_276` | Lung adenocarcinoma | Matched spatial metabolomics and transcriptomics from the same tissue, with pathological annotations | Direct pathological anchors | Godfrey *et al.*, *Angewandte Chemie International Edition* (2025), [doi:10.1002/anie.202502028](https://doi.org/10.1002/anie.202502028) |
| `ccRCC_Y27T` | Clear cell renal cell carcinoma | SpatialMETA-aligned joint SM–ST AnnData (`Y_27T`) | Multi-omics-inferred anchors | [Processed ccRCC data](https://zenodo.org/records/14986870); original study: Hu *et al.*, *Nature Genetics* (2024) |
| `GBM_248T` | Glioblastoma | SpatialMETA-aligned joint SM–ST AnnData (`248_T`) | Multi-omics-inferred anchors | [Processed GBM data](https://doi.org/10.5061/dryad.h70rxwdmj); original study: Ravi *et al.*, *Cancer Cell* (2022) |
| `HCC_P1`, `HCC_P4` | Hepatocellular carcinoma, MVI− and MVI+ | Spatial metabolomics with weak anchors derived from SM clustering and annotations transferred from adjacent histological sections | SM-only transferred-pathology weak anchors | Luo *et al.*, *PLoS Medicine* (2026), “Spatial transcriptomic-metabolic features of tumor foci and tumor capsule in microvascular invasion with hepatocellular carcinoma” |
| `HPD_A1`, `HPD_B1`, `HPD_C1` (`Cd1`–`Cd3` in the manuscript) | Parkinson's disease | Spatial metabolomics, spatial transcriptomics, H&E and disease-region annotations | Disease-region transition | Vicari *et al.*, *Nature Biotechnology* (2024), “Spatial multimodal analysis of transcriptomes and metabolomes in tissues” |

## SpatialMETA inputs

The ccRCC and GBM entries are processed, spatially registered joint SM–ST
inputs rather than unaltered files from the original studies. SpatialMETA
aligns the two modalities to a shared spatial resolution. Additional processed
outputs are available from [Zenodo](https://zenodo.org/records/12528191), and
the alignment code is available in the
[SpatialMETA repository](https://github.com/WanluLiuLab/SpatialMETA).

In these analyses, spatial metabolomics remains the target modality. Matched
transcriptomics supports annotation inference and downstream biological
validation where available.

## Availability notes

- A citation is shown instead of a download button when a stable public
  dataset URL has not yet been verified in this repository.
- Sample identifiers in SIGMA may differ from names used by the source study.
- The package does not include manuscript data, frozen results or large AnnData
  objects.
