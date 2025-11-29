#!/usr/bin/env Rscript

library(argparse)
library(Seurat)
library(CellChat)
library(Matrix)
library(dplyr)
library(tibble)
library(data.table)
library(qs)
library(future)

# Argument parsing
parser <- ArgumentParser()
parser$add_argument("--seurat", required=TRUE, help="Path to the Seurat object (.qs or .rds)")
parser$add_argument("--annotation", required=TRUE, help="Column in meta.data for cell type annotation")
parser$add_argument("--subset_annotation", required=NULL, help="Column in meta.data for cell type annotation that will be used for subset")
parser$add_argument("--cell_types", nargs='*', default=NULL, help="Cell types to keep (optional, space-separated)")
parser$add_argument("--filterDB", type="logical", default=TRUE, help="Whether to filter DB (TRUE/FALSE)")
parser$add_argument("--signal", default="Secreted Signaling", help="Signal type in DB (e.g. 'Secreted Signaling')")
parser$add_argument("--output_dir", required=TRUE, help="Output directory")
parser$add_argument("--output_file", required=TRUE, help="Output filename (without extension)")

args <- parser$parse_args()

# ---- Utility Functions ----

process_data <- function(seurat_object, annotation_col) {
    fov_coordinates <- list()
    for (i in seq_along(Images(seurat_object))) {
        fov_name <- Images(seurat_object)[i]
        fov_data <- as.data.frame(GetTissueCoordinates(seurat_object[[fov_name]][["centroids"]]))
        fov_data$x <- fov_data$x + 10000 * (i - 1)
        fov_data$fov <- fov_name
        fov_coordinates[[fov_name]] <- fov_data
    }
    coordinates <- do.call(rbind, fov_coordinates)
    rownames(coordinates) <- coordinates$cell
    coordinates <- coordinates[, c("x", "y")]
    meta <- seurat_object@meta.data %>% select(labels = all_of(annotation_col))
    exp_matrix <- as.matrix(GetAssayData(seurat_object, slot = "data", assay = "Xenium"))
    exp_sparse <- as(exp_matrix, "dgCMatrix")
    spatial_factors <- data.frame(ratio = 1.0, tol = 10.0)
    return(list(exp_sparse = exp_sparse, meta = meta, coordinates = coordinates, spatial_factors = spatial_factors))
}

process_cellchat <- function(cellchat, workers = 30, future_max_size = 32000 * 1024^2, 
                             signal = "Secreted Signaling", filterDB = TRUE) {
    options(future.globals.maxSize = future_max_size)
    on.exit(future::plan("sequential"), add = TRUE)
    CellChatDB <- CellChatDB.human
    if (isTRUE(filterDB)) {
        CellChatDB.use <- subsetDB(CellChatDB, search = signal, key = "annotation")
    } else {
        CellChatDB.use <- subsetDB(CellChatDB)
    }
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    future::plan("multisession", workers = workers)
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
    cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = FALSE)
    cellchat <- computeCommunProb(
        cellchat, type = "truncatedMean", trim = 0.1,
        distance.use = TRUE, interaction.range = 250, scale.distance = 6.8,
        contact.dependent = TRUE, contact.range = 10
    )
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    return(cellchat)
}

# ---- Load Seurat object ----
if (grepl("\\.qs$", args$seurat)) {
    seurat_obj <- qread(args$seurat)
    seurat_obj
} else if (grepl("\\.rds$", args$seurat)) {
    seurat_obj <- readRDS(args$seurat)
    seurat_obj
} else {
    stop("Unsupported file format for Seurat object.")
}

# ---- Subset by cell types if specified ----
if (!is.null(args$cell_types) && length(args$cell_types) > 0) {
    meta_col <- args$subset_annotation
    if (!(meta_col %in% colnames(seurat_obj@meta.data))) {
        stop(paste("subset_annotation column", meta_col, "not found in seurat_obj@meta.data!"))
    }
    # Print unique values for sanity check
    message("Unique values in subset_annotation: ", paste(unique(seurat_obj@meta.data[[meta_col]]), collapse=", "))
    # Identify cells to keep
    keep_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[meta_col]] %in% args$cell_types]
    message("Number of cells after filtering: ", length(keep_cells))
    # Subset using cell names
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
    seurat_obj
}


# ---- Process data and run CellChat ----
slots <- process_data(seurat_obj, annotation_col = args$annotation)
cellchat <- createCellChat(object = slots$exp_sparse, 
                          meta = slots$meta, 
                          group.by = "labels",
                          datatype = "spatial", 
                          coordinates = as.matrix(slots$coordinates), 
                          spatial.factors = slots$spatial_factors)

cellchat <- process_cellchat(cellchat, filterDB = args$filterDB, signal = args$signal)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# ---- Output ----
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
qsave(cellchat, file = file.path(args$output_dir, paste0(args$output_file, ".qs")))

cat("CellChat analysis completed and saved to", file.path(args$output_dir, paste0(args$output_file, ".qs")), "\n")

