# Data sources

The section notebooks include small saved results for viewing. Raw MSI/RNA
matrices and full fitted models must be obtained separately. PyPI distributes
the software, not the study datasets.

## Studies and downloads

| Data | Sections used by SIGMA | Original study | Data access |
|---|---|---|---|
| Breast cancer and lung adenocarcinoma | BC515 S1, BC515 S2, BC525, LC091, LC276 | Godfrey et al., *Angewandte Chemie International Edition* (2025). [Integrating Ambient Ionization Mass Spectrometry Imaging and Spatial Transcriptomics on the Same Cancer Tissues to Identify RNA–Metabolite Correlations](https://doi.org/10.1002/anie.202502028) | [Harvard Dataverse: GZFCWC](https://doi.org/10.7910/DVN/GZFCWC). BC515 S2 was obtained as `BC_2_515` from the [MAGPIE data archive](https://zenodo.org/records/17789448). |
| Hepatocellular carcinoma | HCC P1, HCC P4 | Luo et al., *PLOS Medicine* (2026). [Spatial transcriptomic-metabolic features of tumor foci and tumor capsule in microvascular invasion with hepatocellular carcinoma](https://doi.org/10.1371/journal.pmed.1004703) | MSI: [OMIX006935](https://ngdc.cncb.ac.cn/omix/release/OMIX006935). ST: [HRA011344](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011344). [Study code](https://github.com/Jingwe-Zhao/HCC_Spatial_Multiomics). |
| Clear cell renal cell carcinoma | Y27T (`Y_27T`) | Hu et al., *Nature Genetics* (2024). [Multi-omic profiling of clear cell renal cell carcinoma identifies metabolic reprogramming associated with disease progression](https://doi.org/10.1038/s41588-024-01662-5) | [SpatialMETA processed joint SM–ST inputs](https://zenodo.org/records/14986870). |
| Colorectal cancer | B_2, C_2 | Nunes et al., *Nature Methods* (2024). [Integration of mass cytometry and mass spectrometry imaging for spatially resolved single-cell metabolic profiling](https://doi.org/10.1038/s41592-024-02392-6) | [Figshare: 23599044, version 1](https://doi.org/10.6084/m9.figshare.23599044.v1). [Registration and integration code](https://github.com/deMirandaLab/MALDI-MSI-IMC). |
| IDH-mutant glioma | SMINT P1; Glioma574 processing object | Kriel et al., *Cell Reports Methods* (2026). [An integrative spatial multi-omic workflow for unified analysis of tumor tissue](https://doi.org/10.1016/j.crmeth.2026.101306) | SMINT: [BioImage Archive S-BIAD1426](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1426), including `processed_data/unified_spe.rds`. Glioma574: [SpaMTP processed data](https://zenodo.org/records/18282678), `brain_tumour_integrated.rds`. |
| Glioblastoma | 248T (`248_T`; notebook ID `GBM`) | Ravi et al., *Cancer Cell* (2022). [Spatially resolved multi-omics deciphers bidirectional tumor-host interdependence in glioblastoma](https://doi.org/10.1016/j.ccell.2022.05.009) | [Original study data on Dryad](https://doi.org/10.5061/dryad.h70rxwdmj). The SIGMA input was spatially aligned through SpatialMETA; see its [processed-data archive](https://zenodo.org/records/12528191) and [code](https://github.com/WanluLiuLab/SpatialMETA). |

For data preparation, also cite [SpatialMETA](https://doi.org/10.1038/s41467-025-63915-z),
[MAGPIE](https://doi.org/10.1038/s41467-025-68003-w), or
[SpaMTP](https://doi.org/10.1038/s41592-026-03140-8) where the corresponding
processed input is used. Follow each archive's access conditions and licence.
Download links identify source archives; they are not automatic SIGMA imports.

## Sample and annotation notes

- **Breast and lung:** paired MSI/RNA with pathological annotations. BC515 S1
  and S2 are sections from the same case. Breast notebook previews use the
  common pathology-reference distance from Figure 3.
- **HCC:** the displayed fits use SM with weak anchors derived from clustering
  and pathology transferred from adjacent sections. The source study's ST
  archive does not mean same-section RNA was used in these fits.
- **Renal and GBM:** the inputs are aligned SM–ST objects, not unmodified raw
  data. Multi-omics-informed anchors differ from independent pathological labels.
- **CRC:** B_2 and C_2 are tissue regions. Labels come from IMC cell phenotypes;
  RNA was not used. Public channel indices have not been reliably matched to
  chemical identities in this analysis.
- **IDH-mutant glioma:** SMINT P1 is reported as grade 3 astrocytoma. Glioma574
  is a separate processing object; its patient-level independence from P1 has
  not been verified. Neither is labelled GBM. Leading-edge tissue is not
  equivalent to healthy tissue.

RNA used to fit a model can provide cellular context but is not independent
validation of that fit. Programme numbers are local to each analysis.

## Files included here

See [data/README.md](../data/README.md) for the saved-result file formats and
[notebooks/README.md](../notebooks/README.md) for viewing and recalculation.
The [section manifest](../notebooks/sections.json) records sample names,
distance definitions, sources and checksums. To train from source data, prepare
the matrices described in the [input guide](input_format.md).
