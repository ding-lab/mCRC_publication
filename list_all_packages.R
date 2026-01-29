#!/usr/bin/env Rscript
# Script to extract and list all R packages used in this project

library(dplyr)

# Get all R files in the project
r_files <- list.files(
  path = getwd(),
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract all library() and require() calls
all_packages <- character()

for (file in r_files) {
  lines <- readLines(file, warn = FALSE)
  # Extract library() calls
  lib_lines <- grep("^library\\(|^require\\(", lines, value = TRUE)
  for (line in lib_lines) {
    # Extract package name
    pkg <- gsub("^library\\(|^require\\(", "", line)
    pkg <- gsub("\\).*$", "", pkg)
    pkg <- gsub('"', "", pkg)
    pkg <- gsub("'", "", pkg)
    pkg <- gsub(",.*$", "", pkg)  # Remove any additional arguments
    pkg <- trimws(pkg)
    if (nchar(pkg) > 0) {
      all_packages <- c(all_packages, pkg)
    }
  }
}

# Get unique packages and sort
unique_packages <- sort(unique(all_packages))

# Categorize packages
single_cell_packages <- c("Seurat", "Signac", "qs", "Matrix", "future", 
                          "harmony", "SeuratWrappers", "reticulate", 
                          "scCustomize", "MAST", "CellChat", "AUCell")

visualization_packages <- c("ComplexHeatmap", "circlize", "ggplot2", "ggpubr", 
                           "ggrepel", "ggalluvial", "patchwork", "cowplot", 
                           "viridis", "RColorBrewer", "grid", "gridExtra", 
                           "EnhancedVolcano")

bioc_packages <- c("AUCell", "GSVA", "GSEABase", "SummarizedExperiment", 
                   "BiocParallel", "BSgenome.Hsapiens.UCSC.hg38", 
                   "EnsDb.Hsapiens.v100", "ensembldb", "GenomeInfoDb", 
                   "GenomicRanges", "chromVAR", "TFBSTools", "JASPAR2020", 
                   "ChIPseeker", "motifmatchr", "AnnotationDbi", "org.Hs.eg.db",
                   "TCGAbiolinks", "biomaRt", "DESeq2", "edgeR", "limma", 
                   "CMScaller", "DelayedMatrixStats", "matrixStats")

tidyverse_packages <- c("tidyverse", "dplyr", "tidyr", "readr", "purrr", 
                       "forcats", "broom", "scales")

statistical_packages <- c("survival", "survminer", "rstatix", "clinfun", 
                         "car", "emmeans", "cutpointr")

ml_packages <- c("glmnet", "caret", "FactoMineR", "factoextra")

# Print results
cat("=== ALL R PACKAGES USED IN PROJECT ===\n\n")
cat("Total unique packages:", length(unique_packages), "\n\n")

cat("=== SINGLE-CELL & SPATIAL ANALYSIS ===\n")
sc_pkgs <- unique_packages[unique_packages %in% single_cell_packages]
cat(paste(sc_pkgs, collapse = ", "), "\n\n")

cat("=== VISUALIZATION ===\n")
viz_pkgs <- unique_packages[unique_packages %in% visualization_packages]
cat(paste(viz_pkgs, collapse = ", "), "\n\n")

cat("=== BIOCONDUCTOR PACKAGES ===\n")
bioc_pkgs <- unique_packages[unique_packages %in% bioc_packages]
cat(paste(bioc_pkgs, collapse = ", "), "\n\n")

cat("=== TIDYVERSE & DATA MANIPULATION ===\n")
tidy_pkgs <- unique_packages[unique_packages %in% tidyverse_packages]
cat(paste(tidy_pkgs, collapse = ", "), "\n\n")

cat("=== STATISTICAL ANALYSIS ===\n")
stat_pkgs <- unique_packages[unique_packages %in% statistical_packages]
cat(paste(stat_pkgs, collapse = ", "), "\n\n")

cat("=== MACHINE LEARNING ===\n")
ml_pkgs <- unique_packages[unique_packages %in% ml_packages]
cat(paste(ml_pkgs, collapse = ", "), "\n\n")

cat("=== OTHER PACKAGES ===\n")
other_pkgs <- unique_packages[!unique_packages %in% c(single_cell_packages, 
                                                      visualization_packages, 
                                                      bioc_packages, 
                                                      tidyverse_packages, 
                                                      statistical_packages, 
                                                      ml_packages)]
cat(paste(other_pkgs, collapse = ", "), "\n\n")

cat("=== COMPLETE ALPHABETICAL LIST ===\n")
cat(paste(unique_packages, collapse = "\n"), "\n")

# Save to file
write.table(data.frame(Package = unique_packages), 
            file = "packages_list.txt", 
            row.names = FALSE, 
            quote = FALSE)

cat("\n\nPackage list also saved to: packages_list.txt\n")
