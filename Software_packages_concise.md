# Software Packages - Concise Version for Methods

## R Packages (with versions)

**Single-cell and spatial transcriptomics analysis:**
- Seurat v5.3.0, Signac v1.14.0, qs v0.27.3, Matrix v1.7-3, future v1.49.0, harmony v1.2.3, reticulate v1.42.0, scCustomize v3.0.1, MAST v1.32.0

**Gene set enrichment:**
- AUCell v1.28.0, GSVA v2.0.7, GSEABase v1.68.0

**Chromatin accessibility (snATAC-seq):**
- chromVAR v1.28.0, TFBSTools v1.44.0, JASPAR2020 v0.99.10, ChIPseeker v1.42.1, motifmatchr v1.28.0, SummarizedExperiment v1.36.0, BiocParallel v1.40.0, BSgenome.Hsapiens.UCSC.hg38 v1.4.5, EnsDb.Hsapiens.v100 v0.0.1, EnsDb.Hsapiens.v86 v2.99.0, ensembldb v2.30.0, GenomeInfoDb v1.42.0, GenomicRanges v1.58.0

**Cell-cell communication:**
- CellChat v2.2.0

**Visualization:**
- ComplexHeatmap v2.22.0, circlize v0.4.16, ggplot2 v3.5.2, ggpubr v0.6.0, ggrepel v0.9.6, patchwork v1.3.0, cowplot v1.1.3, viridis v0.6.5, RColorBrewer v1.1-3, EnhancedVolcano v1.24.0

**Data manipulation:**
- tidyverse v2.0.0 (includes dplyr v1.1.4, tidyr v1.3.1, readr v2.1.5, purrr v1.0.4, forcats v1.0.0), data.table v1.17.4

**Statistical analysis:**
- survival v3.8-3, survminer v0.5.0, rstatix v0.7.2, clinfun v1.1.5, broom v1.0.8, car v3.1-3, emmeans v1.11.2-8, scales v1.4.0

**Differential expression:**
- DESeq2 v1.46.0, edgeR v4.4.2, limma v3.62.2

**TCGA analysis:**
- TCGAbiolinks v2.34.1, biomaRt v2.62.1, AnnotationDbi v1.68.0, org.Hs.eg.db v3.20.0

**Machine learning:**
- glmnet v4.1-8, caret v6.0-94, FactoMineR v2.12, factoextra v1.0.7

## Python Packages

**Spatial transcriptomics and single-cell analysis:**
- scanpy, anndata, squidpy, spatialdata, spatialdata_io

**Clustering and integration:**
- banksy, harmonypy, secuer

**Dimensionality reduction:**
- umap, sklearn (scikit-learn)

**Image processing:**
- Morph (install from https://github.com/ding-lab/morph), skimage (scikit-image), scipy

**Data manipulation:**
- numpy, pandas, matplotlib

## Software Versions

- **R**: 4.4.3 (via `seurat5_env` conda environment)
- **Seurat**: v5.3.0
- **Python**: 
  - 3.13.3 (in `seurat5_env` for Xenium Banksy clustering and spatial analysis)
  - 3.13.7 (in `3d-analysis` for 3D reconstruction)
  - 3.10.19 (in `morph_env` for morphological annotation)

**Note**: Exact package versions are specified in conda environment YAML files in `envs/`. For complete documentation, see `Software_packages_list.md`.
