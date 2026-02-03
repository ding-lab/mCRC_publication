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
  - CellChat object creation

- **`D_AUCell/`**: Gene set enrichment analysis using AUCell
  - Epithelial gene sets (snRNA)
  - TMR region analysis (Visium, Xenium)
  - Fibroblast gene sets (Xenium)
  - Pseudobulk boundary analysis (Visium)

### `B_Data_collection_visualization/`
Scripts for generating data collection overview figures.

- **`F1A_case_level_comparisons.R`**: Figure 1A - case-level comparisons and sample statistics
- **`S1A_data_collections.R`**: Supplementary Figure 1A - data collection heatmaps
- **`S3_Xenium_Tissue_Domain_Spatial_Plot.R`**: Supplementary Figure 3 - Xenium spatial tissue domain plots

### `C_Somatic_mutation_analysis/`
Somatic mutation analysis and visualization.

- **`F5A_KRAS_mut_freq_barplot.R`**: Figure 5A - KRAS mutation frequency barplot
- **`S2B_oncoplot.R`**: Supplementary Figure 2B - oncoplot visualization

### `D_Bulk_RNAseq/`
Bulk RNA-seq analysis.

- **`S3C_bulkRNA_CMS.R`**: Supplementary Figure 3C - Consensus Molecular Subtype (CMS) classification

### `E_Global_tissue_domain_analysis/`
Global tissue domain analysis scripts.

- **`F2A_NB_celltype_heatmap.R`**: Figure 2A - neighborhood cell type heatmap

### `F_Putative_proinvasinve_tumor/`
Analysis of putative pro-invasive tumor features.

- **`F3E_S5F_snATAC_cell_type_correlation.R`**: Figure 3E and Supplementary Figure 5F - snATAC cell type correlation
- **`F3J_S5H_TCGA_survival.R`**: Figure 3J and Supplementary Figure 5H - TCGA survival analysis

### `G_TMR_PSI_interaction/`
Tumor Margin Region (TMR) and Peri-Stromal Interface (PSI) interaction analysis.

- **`F4C_S6BC_Boxplot_tumor_subtype_prop_by_histology_regions.R`**: Figure 4C and Supplementary Figure 6B/C - tumor subtype proportions
- **`F4D_Heatmap_selected_geneset_TMR_regions.R`**: Figure 4D - gene set heatmap in TMR regions
- **`S6A_Correlation_heatmap_TMR_by_histological_layers.R`**: Supplementary Figure 6A - TMR correlation heatmap
- **`S6D_fibroblast_subtypes_by_GTR_Organ.R`**: Supplementary Figure 6D - fibroblast subtypes by GTR and organ
- **`S6E_Barplot_Fibroblast_by_TMR.R`**: Supplementary Figure 6E - fibroblast barplot by TMR
- **`S6F_Correlation_heatmap_PSI_fibroblast_by_histological_layers.R`**: Supplementary Figure 6F - PSI fibroblast correlation heatmap
- **`S6G_snRNAseq_Fibroblast_Tumor_subtype_Association.R`**: Supplementary Figure 6G - snRNAseq fibroblast-tumor subtype association
- **`S6H_Visium_Tumor_Boundary.R`**: Supplementary Figure 6H - Visium tumor boundary analysis
- **`S6I_TMR_boundary_Tumor_fibroblast_subtype_association.R`**: Supplementary Figure 6I - TMR boundary tumor-fibroblast subtype association

### `H_Organ_specific_adaption/`
Organ-specific adaptation analysis.

- **`F5E_Liver_Lung_CCI_differences_tumor_cell_targeted.R`**: Figure 5E - liver vs lung cell-cell interaction differences
- **`F5F_S7B_HGF_EGF_SLIT_example_spatial_plot.R`**: Figure 5F and Supplementary Figure 7B - HGF/EGF/SLIT spatial plots
- **`F5G_S7CDG_snRNA_selected_CCI_receptors_ligands_dotplot_violinplot.R`**: Figure 5G and Supplementary Figure 7C/D/G - snRNA CCI receptors/ligands visualization
- **`F5G_S7F_CellChat_selected_tumor_recepters_boxplot.R`**: Figure 5G and Supplementary Figure 7F - CellChat tumor receptors boxplot
- **`S7D_Xenium_selected_CCI_receptors_ligands_dotplot.R`**: Supplementary Figure 7D - Xenium CCI receptors/ligands dotplot
- **`S7E_CellChat_selected_tumor_recepters_dotplot.R`**: Supplementary Figure 7E - CellChat tumor receptors dotplot
- **`S7H_Xenium_Macrophage_HBEGF_dotplot.R`**: Supplementary Figure 7H - Xenium macrophage HBEGF dotplot
- **`S7I_Xenium_Macrophage_subtype_Liver_vs_Lung_GTR_proportion.R`**: Supplementary Figure 7I - Xenium macrophage subtype liver vs lung GTR proportion

### `I_Liver_mCRC_growth_pattern/`
Liver metastasis growth pattern analysis.

### `J_Angiogenesis_3D_model/`
Angiogenesis 3D modeling scripts.

## Naming Convention

Scripts are named with figure identifiers:
- **`F#X_`**: Main figure panels (e.g., `F1A_`, `F2A_`, `F3E_`, `F4C_`, `F5A_`)
- **`S#X_`**: Supplementary figure panels (e.g., `S1A_`, `S2B_`, `S3C_`, `S6A_`, `S7D_`)
- Scripts generating multiple figures use format: `F#X_S#Y_` (e.g., `F3E_S5F_`, `F4C_S6BC_`)

## Requirements

### Software and Environments

- **R Version**: R 4.4.3 (via `seurat5_env` conda environment)
- **Python Versions**: 
  - Python 3.13.3 (for R analysis environment)
  - Python 3.13.7 (for 3D analysis)
  - Python 3.10.19 (for morphological annotation)
- **Key R Packages**: Seurat v5.3.0, tidyverse, ComplexHeatmap, AUCell, qs, future, CellChat
- **System**: Linux (tested on RHEL 7), sufficient RAM for large Seurat objects (30GB+)
- **Data Access**: Scripts require access to specific server file paths and data directories

### Environment Setup

**Recommended**: Use the conda environment files for reproducible setup:

```bash
# For R analysis (most scripts)
conda env create -f envs/seurat5_env.yml
conda activate seurat5_env

# For morphological annotation
conda env create -f envs/morph_env.yml
conda activate morph_env
# Install morph from GitHub
pip install git+https://github.com/ding-lab/morph.git

# For 3D reconstruction
conda env create -f envs/3d-analysis_env.yml
conda activate 3d-analysis
```

For detailed package documentation, see:
- `PACKAGE_DOCUMENTATION.md` - **Start here**: Guide to all package documentation files
- `Software_packages_list.md` - Complete package documentation with descriptions
- `Software_packages_concise.md` - Concise version for methods section
- `envs/README.md` - Environment setup instructions

## Bioinformatics Tools

Many bioinformatics tools were used in the course of this work. All tools written and/or published by the authors are freely available at our public GitHub repository (https://github.com/ding-lab/), including:

- **somaticwrapper** - Variant calling pipeline: https://github.com/ding-lab/somaticwrapper
- **10Xmapping** - Code for mutation mapping from bulk to single cells: https://github.com/ding-lab/10Xmapping
- **ffpefiltering** - Pipeline for FFPE WES filtering: https://github.com/ding-lab/ffpefiltering
- **pecgs-pipeline** - Bulk RNA-seq alignment and transcript counting: https://github.com/ding-lab/pecgs-pipeline
- **Morph** - Spatial transcriptomics toolset for tumor boundary detection: https://github.com/ding-lab/morph

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
