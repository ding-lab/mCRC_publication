# seurat_scvi_env
# Script Author: Evan P.
# Date: August 10, 2024
#
# Usage:
# -----
# To run this script from the command line, use the following command:
#
#   Rscript run_SCVI_integration_v2.R <seurat_object_path> <output_folder_path> <file_prefix> <cell_type> [group_by_column]
#
# <seurat_object_path>: Full path to the Seurat object RDS file (e.g., "/path/to/seurat_object.rds")
# <output_folder_path>: Directory where the output files (e.g., integrated object, plots) will be saved (e.g., "/path/to/output/folder")
# <file_prefix>: Prefix to be added to the names of the output files (e.g., "mcrc_v4_")
# <cell_type>: Cell type label to be used in the cluster naming.
# [group_by_column]: Optional column name to group the UMAP plot by.
#
# Example:
# -----
#   Rscript run_SCVI_integration_v2.R /diskmnt/Projects/Users/Evan.p/data/merged_mcrc.rds /diskmnt/Projects/Users/Evan.p/output mcrc_v4_ MyCellType orig.ident
#
# Ensure that the necessary R packages and the conda environment are set up correctly before running the script.


# Load necessary libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(reticulate)
library(SeuratWrappers)
library(optparse)

# Increase memory limit
options(future.globals.maxSize = 1e9)

# Set conda environment
use_condaenv("seurat_scvi_env", required = TRUE)

# Define function to perform scVI integration
perform_scvi_integration <- function(seurat_obj_path, output_folder, file_prefix, cell_type, batch_column = 'orig.ident', group_by_column = NULL) {
  
  print('Load Seurat object')
  # Load the Seurat object from the provided path
  seurat_obj <- readRDS(seurat_obj_path)
  
  # Ensure output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  print('Remove samples with only one cell')
  # Identify and remove samples (orig.ident) with only one cell
  seurat_obj$batch <- seurat_obj[[batch_column]]
  cell_counts <- table(seurat_obj$batch)
  valid_samples <- names(cell_counts[cell_counts > 1])
  seurat_obj <- subset(seurat_obj, cells = WhichCells(seurat_obj, expression = batch %in% valid_samples))
  print(paste0(valid_samples, ' with cell number > 1.'))
  
  print(paste('Split Layers', batch_column))
  # Convert RNA assay to Assay5
  #seurat_obj[["RNA5"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay5")
  DefaultAssay(seurat_obj) <- "RNA"
  Idents(object = seurat_obj) <- batch_column
  
  # Split assay by identity
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$batch)
  
  print('RNA normalization')
  # Standard Seurat preprocessing
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 20)
  
  print('Run scVI integration')
  # Perform scVI integration
  start.time <- Sys.time()
  integrated_obj <- IntegrateLayers(
    object = seurat_obj, method = scVIIntegration,
    new.reduction = "integrated.scvi",
    conda_env = "/data/epeng/miniforge3/envs/seurat_scvi_env/",
    verbose = TRUE
  )
  end.time <- Sys.time()
  
  # Calculate time taken for integration
  time.taken <- round(end.time - start.time, 2)
  print(paste("Time taken for scVI integration:", time.taken))

  # Run UMAP and clustering
  print('Run UMAP')
  integrated_obj <- RunUMAP(integrated_obj, reduction = "integrated.scvi", dims = 1:20, reduction.name = paste0(cell_type, "_umap.scvi"))
  integrated_obj <- FindNeighbors(integrated_obj, reduction = "integrated.scvi", dims = 1:20)
  
  print('Run clustering')
  # Perform clustering at resolutions from 0.1
  resolutions <- seq(0.1, 2, by = 0.1)
  for (res in resolutions) {
    cluster_name <- paste0(cell_type, "_clusters_", res)
    integrated_obj <- FindClusters(integrated_obj, resolution = res, cluster.name = cluster_name, graph = 'RNA_snn')
    
    # Save a DimPlot of UMAP for each resolution
    dimplot_path <- file.path(output_folder, paste0(file_prefix, "_umap_", cluster_name, ".png"))
    p <- DimPlot(integrated_obj, reduction = paste0(cell_type, "_umap.scvi"), group.by = cluster_name, label = TRUE)
    ggsave(filename = dimplot_path, plot = p)
  }

  # If a group_by_column is provided, generate a UMAP grouped by this column
  if (!is.null(group_by_column)) {
    print(paste('Generating UMAP grouped by:', group_by_column))
    group_by_dimplot_path <- file.path(output_folder, paste0(file_prefix, "_umap_grouped_by_", group_by_column, ".png"))
    p_group_by <- DimPlot(integrated_obj, reduction = paste0(cell_type, "_umap.scvi"), group.by = group_by_column, label = TRUE)
    ggsave(filename = group_by_dimplot_path, plot = p_group_by)
  }

  DefaultAssay(integrated_obj) <- 'RNA'
  #integrated_obj[["RNA5"]] <- NULL
  saveRDS(integrated_obj, file = file.path(output_folder, paste0(file_prefix, '.rds')))
}

# Define command-line options
option_list <- list(
  make_option(c("-s", "--seurat_obj_path"), type = "character", default = NULL, help = "Path to Seurat object RDS file", metavar = "character"),
  make_option(c("-o", "--output_folder"), type = "character", default = NULL, help = "Output folder path", metavar = "character"),
  make_option(c("-p", "--file_prefix"), type = "character", default = NULL, help = "Prefix for output files", metavar = "character"),
  make_option(c("-c", "--cell_type"), type = "character", default = NULL, help = "Cell type to process", metavar = "character"),
  make_option(c("-g", "--group_by_column"), type = "character", default = NULL, help = "Column name to group UMAP by", metavar = "character"),
  make_option(c("-b", "--batch"), type = "character", default = "orig.ident", help = "Column name for batch assignment (default: orig.ident)", metavar = "character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure required arguments are provided
if (is.null(opt$seurat_obj_path) || is.null(opt$output_folder) || is.null(opt$file_prefix) || is.null(opt$cell_type)) {
  print_help(opt_parser)
  stop("Please provide the required arguments.")
}

# Run the integration
perform_scvi_integration(opt$seurat_obj_path, opt$output_folder, opt$file_prefix, opt$cell_type, opt$batch, opt$group_by_column)