# mCRC Manuscript Scripts

This repository contains all R scripts and analysis pipelines used for the manuscript: **"Spatial Zonation of Tumor Cell States, Stromal-Immune Networks, and Growth Patterns in Metastatic Colorectal Cancer"**. Scripts are organized by analysis category and figure generation.

## Folder Structure

- **`A_Data_preprocessing/`** - Data preprocessing pipelines
  - `A_snMultiome/` - Single-nucleus multiome (RNA + ATAC) processing
  - `B_Xenium/` - Xenium spatial transcriptomics processing
  - `C_CellChat/` - Cell-cell communication analysis
  - `D_AUCell/` - Gene set enrichment analysis
- **`B_Data_collection_visualization/`** - Data collection overview figures
- **`C_Somatic_mutation_analysis/`** - Mutation analysis and visualization
- **`D_Bulk_RNAseq/`** - Bulk RNA-seq analysis
- **`E_NB_and_Functional_Zone/`** - Neighborhood and functional zone analysis
- **`F_Tumor_subtypes/`** - Tumor subtype analysis
- **`G_TMR_PSI_interaction/`** - Tumor Margin Region and Peri-Stromal Interface analysis
- **`H_Organ_specific_adaption/`** - Organ-specific adaptation analysis

## Requirements

- **R**: 4.4.3 (via `seurat5_env` conda environment)
- **Python**: 3.13.3 (R analysis), 3.13.7 (3D analysis), 3.10.19 (morphological annotation)
- **Key packages**: Seurat v5.3.0, tidyverse, ComplexHeatmap, AUCell, CellChat
- **System**: Linux (tested on RHEL 7), 30GB+ RAM recommended

### Quick Start

```bash
# R analysis (most scripts)
conda env create -f envs/seurat5_env.yml
conda activate seurat5_env

# Morphological annotation
conda env create -f envs/morph_env.yml
conda activate morph_env
pip install git+https://github.com/ding-lab/morph.git

# 3D reconstruction
conda env create -f envs/3d-analysis_env.yml
conda activate 3d-analysis
```

See `envs/README.md` for build times and detailed setup instructions. For package documentation, see `PACKAGE_DOCUMENTATION.md`.

## Bioinformatics Tools

Many bioinformatics tools were used in the course of this work. All tools written and/or published by the authors are freely available at our public GitHub repository (https://github.com/ding-lab/), including:

- **somaticwrapper** - Variant calling pipeline: https://github.com/ding-lab/somaticwrapper
- **10Xmapping** - Code for mutation mapping from bulk to single cells: https://github.com/ding-lab/10Xmapping
- **ffpefiltering** - Pipeline for FFPE WES filtering: https://github.com/ding-lab/ffpefiltering
- **pecgs-pipeline** - Bulk RNA-seq alignment and transcript counting: https://github.com/ding-lab/pecgs-pipeline
- **Morph** - Spatial transcriptomics toolset for tumor boundary detection: https://github.com/ding-lab/morph

## Usage

Scripts are named with figure identifiers (`F#X_` for main figures, `S#X_` for supplementary). Set working directory to script location and source:

```r
setwd('/path/to/mCRC_Manuscript_Script/[folder_name]')
source('script_name.R')
```

**Note**: Scripts require access to specific server file paths and data directories. Many depend on custom functions in `/diskmnt/Users2/epeng/Projects/mCRC/scripts/`.

---

**Manuscript**: Spatial Zonation of Tumor Cell States, Stromal-Immune Networks, and Growth Patterns in Metastatic Colorectal Cancer  
**Status**: Under Review
