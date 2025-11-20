# mCRC Manuscript Scripts

This repository contains all R scripts and analysis pipelines used for the mCRC (metastatic Colorectal Cancer) manuscript. Scripts are organized by analysis category and figure generation.

## Folder Structure

### `A_Data_preprocessing/`
Data preprocessing pipelines for single-cell and spatial transcriptomics data.

- **`A_snMultiome/`**: Single-nucleus multiome (RNA + ATAC) processing
  - Individual object processing (CellRanger, Seurat pipeline)
  - Integrated object processing (SCVI integration)
  - Cell type annotation (epithelial, mesenchymal/endothelial, myeloid, T/NK cells)

- **`B_Xenium/`**: Xenium spatial transcriptomics processing
  - Individual object processing
  - Integrated object processing (sketch integration)
  - Cell type annotation
  - Banksy clustering
  - Morphological annotation (tissue domain classification)

- **`C_CellChat/`**: Cell-cell communication analysis

- **`D_AUCell/`**: Gene set enrichment analysis using AUCell
  - Epithelial gene sets (snRNA)
  - TMR region analysis (Visium, Xenium)
  - Fibroblast gene sets (Xenium)

### `B_Data_collection_visualization/`
Scripts for generating data collection overview figures.

- **`F1A_case_level_comparisons.R`**: Figure 1A - case-level comparisons and sample statistics
- **`SF1A_data_collections.R`**: Supplementary Figure 1A - data collection heatmaps
- **`SF3_Xenium_Tissue_Domain_Spatial_Plot.R`**: Supplementary Figure 3 - Xenium spatial tissue domain plots

### `C_Somatic_mutation_analysis/`
Somatic mutation analysis and visualization.

- **`F5A_KRAS_mut_freq_barplot.R`**: Figure 5A - KRAS mutation frequency barplot
- **`SF2B_oncoplot.R`**: Supplementary Figure 2B - oncoplot visualization

### `D_Bulk_RNAseq/`
Bulk RNA-seq analysis.

- **`SF3C_bulkRNA_CMS.R`**: Supplementary Figure 3C - Consensus Molecular Subtype (CMS) classification

### `E_Global_tissue_domain_analysis/`
Global tissue domain analysis scripts.

### `F_Putative_proinvasinve_tumor/`
Analysis of putative pro-invasive tumor features.

### `G_TMR_PSI_interaction/`
Tumor Margin Region (TMR) and Peri-Stromal Interface (PSI) interaction analysis.

- **`F4C_Boxplot_tumor_subtype_prop_by_histology_regions.R`**: Figure 4C - tumor subtype proportions
- **`F4D_Heatmap_selected_geneset_TMR_regions.R`**: Figure 4D - gene set heatmap in TMR regions
- **`S6A_Correlation_heatmap_TMR_by_histological_layers.R`**: Supplementary Figure 6A - TMR correlation heatmap
- **`S6B_Correlation_heatmap_PSI_fibroblast_by_histological_layers.R`**: Supplementary Figure 6B - PSI fibroblast correlation heatmap

### `H_Organ_specific_adaption/`
Organ-specific adaptation analysis.

### `I_Liver_mCRC_growth_pattern/`
Liver metastasis growth pattern analysis.

### `J_Angiogenesis_3D_model/`
Angiogenesis 3D modeling scripts.

## Naming Convention

Scripts are named with figure identifiers:
- **`F#X_`**: Main figure panels (e.g., `F1A_`, `F4C_`, `F5A_`)
- **`SF#X_`** or **`S#X_`**: Supplementary figure panels (e.g., `SF1A_`, `SF2B_`, `S6A_`)

## Requirements

- **R Version**: R 4.0 or higher
- **Key R Packages**: Seurat, tidyverse, ComplexHeatmap, AUCell, qs, future
- **System**: Linux (tested on RHEL 7), sufficient RAM for large Seurat objects (30GB+)
- **Data Access**: Scripts require access to specific server file paths and data directories

## Usage

Each folder may contain its own README with specific instructions. For figure generation scripts, set the working directory to the script location and source the R file:

```r
setwd('/path/to/mCRC_Manuscript_Script/[folder_name]')
source('script_name.R')
```

## Notes

- Scripts read from Google Sheets (read-only, no authentication required)
- Output files are typically saved to the current working directory
- Many scripts depend on custom functions located in `/diskmnt/Users2/epeng/Projects/mCRC/scripts/`
- Data files are stored in `/diskmnt/Projects/MetNet_analysis/Colorectal/data/`

---

**Manuscript**: mCRC (metastatic Colorectal Cancer)  
**Status**: Pre-submission / Under Review

