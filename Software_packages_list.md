# Software and R/Python Packages Used for Single-Cell and Xenium Analysis

This document lists all software packages and libraries used for single-cell RNA-seq, single-nucleus multiome (snRNA-seq + snATAC-seq), and Xenium spatial transcriptomics analyses in this project, along with their versions.

## Core Single-Cell Analysis Packages

### R Packages

#### Primary Single-Cell Analysis
- **Seurat** v5.3.0 - Primary framework for single-cell and spatial transcriptomics analysis
- **Signac** v1.14.0 - Single-cell chromatin accessibility analysis (for snATAC-seq)
- **qs** v0.27.3 - Fast serialization for R objects (used for storing large Seurat objects)
- **Matrix** v1.7-3 - Sparse and dense matrix classes and methods
- **future** v1.49.0 - Parallel processing framework for R

#### Integration and Clustering
- **harmony** v1.2.3 - Integration of single-cell datasets
- **reticulate** v1.42.0 - R interface to Python (for SCVI integration)

#### Cell Type Annotation and Analysis
- **scCustomize** v3.0.1 - Custom plotting functions for Seurat objects
- **MAST** v1.32.0 - Model-based Analysis of Single-cell Transcriptomics (differential expression)

### Python Packages (for Xenium and Integration)

**Note**: Python packages are distributed across multiple conda environments. See `envs/README.md` for environment-specific package lists.

#### Spatial Transcriptomics
- **scanpy** - Single-cell analysis in Python (used for Xenium data processing)
  - v1.11.2 in `seurat5_env` (via pip)
  - v1.9.5 in `banksy_env` (for Banksy clustering)
  - v1.11.4 in `3d-analysis_env` (for 3D reconstruction)
- **anndata** v0.11.4 - Annotated data objects for single-cell genomics (in `seurat5_env` via pip)
- **squidpy** - Spatial single-cell omics analysis
  - v1.6.5 in `banksy_env` (for Banksy clustering)
  - v1.2.2 in `3d-analysis_env` (for 3D reconstruction)
- **spatialdata** v0.4.0 - Spatial omics data structures (in `banksy_env`)
- **spatialdata_io** v0.2.0 - I/O for spatial omics data (including Xenium) (in `banksy_env`)

#### Clustering and Integration
- **banksy** - Spatial clustering algorithm (Banksy clustering for Xenium). Used via `banksy_env` environment
- **harmonypy** v0.0.10 - Harmony integration in Python (in `seurat5_env` via pip, also in `banksy_env`)
- **secuer** v1.1 - Additional spatial analysis tools (in `banksy_env` via pip)

#### Dimensionality Reduction
- **umap-learn** - Uniform Manifold Approximation and Projection
  - v0.5.7 in `seurat5_env` (via pip)
  - v0.5.4 in `banksy_env`
  - v0.5.9.post2 in `3d-analysis_env`
- **sklearn** (scikit-learn) v1.7.0 - Machine learning, including PCA (in `seurat5_env` via pip)

#### Image Processing (for Morphological Annotation)
- **Morph** - Spatial transcriptomics toolset for tumor boundary detection and morphological operations. **Installation**: Must be installed from GitHub after setting up `morph_env`:
  ```bash
  pip install git+https://github.com/ding-lab/morph.git
  ```
  See: https://github.com/ding-lab/morph
- **skimage** (scikit-image) v0.25.2 - Image processing
- **scipy** v1.16.2 - Scientific computing

#### Data Manipulation and Visualization
- **numpy** v2.3.0 - Numerical computing
- **pandas** v2.2.3 - Data manipulation and analysis
- **matplotlib** v3.10.3 - Plotting library

## Gene Set Enrichment and Pathway Analysis

- **AUCell** v1.28.0 - Gene set enrichment scoring using AUC (Area Under the Curve)
- **GSVA** v2.0.7 - Gene Set Variation Analysis
- **GSEABase** v1.68.0 - Base classes and methods for Gene Set Enrichment Analysis
- **gProfileR** v0.7.0 - Functional enrichment analysis

## Chromatin Accessibility Analysis (snATAC-seq)

- **chromVAR** v1.28.0 - Chromatin variation analysis
- **TFBSTools** v1.44.0 - Transcription factor binding site analysis
- **JASPAR2020** v0.99.10 - Transcription factor binding site database
- **ChIPseeker** v1.42.1 - ChIP peak annotation and visualization
- **motifmatchr** v1.28.0 - Motif matching in genomic regions
- **SummarizedExperiment** v1.36.0 - Container for matrix-like genomic data
- **BiocParallel** v1.40.0 - Parallel evaluation for Bioconductor
- **BSgenome.Hsapiens.UCSC.hg38** v1.4.5 - Human reference genome (hg38)
- **EnsDb.Hsapiens.v100** v0.0.1 - Ensembl database for human (v100)
- **EnsDb.Hsapiens.v86** v2.99.0 - Ensembl database for human (v86)
- **ensembldb** v2.30.0 - Ensembl database interface
- **GenomeInfoDb** v1.42.0 - Genome information database
- **GenomicRanges** v1.58.0 - Representation and manipulation of genomic intervals

## Cell-Cell Communication

- **CellChat** v2.2.0 - Analysis of cell-cell communication from single-cell data

## Visualization

- **ComplexHeatmap** v2.22.0 - Advanced heatmap visualization
- **circlize** v0.4.16 - Circular visualization
- **ggplot2** v3.5.2 - Grammar of graphics plotting
- **ggpubr** v0.6.0 - Publication-ready plots based on ggplot2
- **ggrepel** v0.9.6 - Text and label geoms for ggplot2
- **ggalluvial** v0.12.5 - Alluvial plots
- **ggrastr** v1.0.2 - Rasterization for ggplot2
- **patchwork** v1.3.0 - Composing plots
- **cowplot** v1.1.3 - Publication-ready theme for ggplot2
- **viridis** v0.6.5 - Color scales for visualization
- **RColorBrewer** v1.1-3 - Color palettes
- **grid** v4.4.3 - Grid graphics system (base R)
- **gridExtra** v2.3 - Additional grid graphics functions
- **EnhancedVolcano** v1.24.0 - Enhanced volcano plots

## Data Manipulation and Analysis

- **tidyverse** v2.0.0 - Collection of R packages for data science
- **dplyr** v1.1.4 - Data manipulation
- **tidyr** v1.3.1 - Tidy data
- **readr** v2.1.5 - Read rectangular data
- **purrr** v1.0.4 - Functional programming tools
- **data.table** v1.17.4 - Fast data manipulation
- **forcats** v1.0.0 - Tools for working with categorical variables
- **tibble** v3.2.1 - Modern data frames
- **magrittr** v2.0.3 - Forward pipe operator

## Statistical Analysis

- **survival** v3.8-3 - Survival analysis
- **survminer** v0.5.0 - Survival analysis visualization
- **rstatix** v0.7.2 - Pipe-friendly framework for basic statistical tests
- **clinfun** v1.1.5 - Clinical trial design and analysis
- **broom** v1.0.8 - Convert statistical objects to tidy tibbles
- **car** v3.1-3 - Companion to Applied Regression
- **emmeans** v1.11.2-8 - Estimated marginal means
- **scales** v1.4.0 - Scale functions for visualization and formatting
- **cutpointr** v1.2.1 - Determine and evaluate optimal cutpoints

## Differential Expression and Bulk RNA-seq

- **DESeq2** v1.46.0 - Differential gene expression analysis
- **edgeR** v4.4.2 - Empirical analysis of digital gene expression data
- **limma** v3.62.2 - Linear models for microarray and RNA-seq data
- **CMScaller** v2.0.1 - Consensus Molecular Subtype (CMS) classification for colorectal cancer

## TCGA and Public Data Analysis

- **TCGAbiolinks** v2.34.1 - Download and analyze TCGA data
- **biomaRt** v2.62.1 - Interface to BioMart databases
- **AnnotationDbi** v1.68.0 - Annotation database interface
- **org.Hs.eg.db** v3.20.0 - Human genome-wide annotation database

## Machine Learning and Classification

- **glmnet** v4.1-8 - Lasso and elastic-net regularized generalized linear models
- **caret** v6.0-94 - Classification and regression training
- **FactoMineR** v2.12 - Multivariate exploratory data analysis
- **factoextra** v1.0.7 - Extract and visualize results of multivariate data analyses

## Additional Utilities

- **googlesheets4** v1.1.1 - Read Google Sheets from R
- **optparse** v1.7.5 - Command-line option parser
- **argparse** v2.2.5 - Command-line argument parsing (Python)
- **logging** - Logging facility (Python, standard library)
- **pickle** - Object serialization (Python, standard library)
- **csv** - CSV file handling (Python, standard library)
- **os** - Operating system interface (Python, standard library)
- **time** - Time-related functions (Python, standard library)

## Additional Bioconductor Packages

- **DelayedMatrixStats** v1.28.1 - Delayed matrix operations
- **matrixStats** v1.5.0 - Matrix statistics
- **genefilter** v1.88.0 - Methods for filtering genes from microarray experiments

## Package Availability

All packages listed above are available in the conda environments specified in `envs/`. Some packages may be optional dependencies or used only in specific analysis workflows. For exact package versions and availability, refer to the conda environment YAML files in `envs/`.

## Version Information

- **R Version**: R 4.4.3 (see `envs/seurat5_env.yml`)
- **Seurat**: v5.3.0
- **Python**: 
  - v3.13.3 in `seurat5_env` (for Xenium Banksy clustering and spatial analysis)
  - v3.13.7 in `3d-analysis` (for 3D reconstruction)
  - v3.10.19 in `morph_env` (for morphological annotation)

**Note**: Exact package versions are specified in the conda environment YAML files in `envs/`. For reproducible analysis, always use the conda environments rather than installing packages individually.

## System Requirements

- **Operating System**: Linux (tested on RHEL 7)
- **RAM**: 30GB+ recommended for large Seurat objects
- **Storage**: Sufficient space for large spatial transcriptomics datasets

## Installation Notes

**Recommended**: Use the conda environment files for reproducible setup:

```bash
# For R analysis (Seurat, Bioconductor, etc.)
conda env create -f envs/seurat5_env.yml
conda activate seurat5_env

# For morphological annotation
conda env create -f envs/morph_env.yml
conda activate morph_env
# Install morph from GitHub
pip install git+https://github.com/ding-lab/morph.git

# For 3D reconstruction and spatial analysis
conda env create -f envs/3d-analysis_env.yml
conda activate 3d-analysis
```

**Alternative manual installation** (not recommended for reproducibility):

R packages can be installed from CRAN or Bioconductor:
```r
# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Seurat", "Signac", "AUCell", "ComplexHeatmap", ...))

# CRAN packages
install.packages(c("tidyverse", "qs", "future", ...))
```

Python packages can be installed via pip or conda:
```bash
pip install scanpy anndata squidpy spatialdata spatialdata-io banksy harmonypy
```

See `envs/README.md` for more details on environment setup.

## Package Documentation Structure

This repository contains multiple documentation files for packages and environments. **See `PACKAGE_DOCUMENTATION.md` for a complete guide to all documentation files.**

1. **Conda Environment Files** (`envs/` directory):
   - `seurat5_env.yml` - Complete R analysis environment (Seurat v5, Bioconductor, Python packages)
   - `morph_env.yml` - Morphological annotation environment
   - `3d-analysis_env.yml` - 3D reconstruction and spatial analysis environment
   - These YAML files are the **source of truth** for reproducible environment setup

2. **Version Tables** (machine-readable):
   - `packages_versions_table.md` - R packages with versions (markdown table)
   - `python_packages_versions_table.md` - Python packages with versions (markdown table)

3. **Documentation Files**:
   - `Software_packages_list.md` (this file) - Detailed documentation with descriptions
   - `Software_packages_concise.md` - Concise version for methods section

## Complete Package List with Versions

For complete package lists with exact versions:
- **Complete environments**: See `envs/` directory for conda environment YAML files (source of truth)
- **Quick reference**: See `Software_packages_concise.md` for a concise list organized by category

To recreate the analysis environment, use the conda environment files in `envs/`:
```bash
conda env create -f envs/seurat5_env.yml
conda activate seurat5_env
```
