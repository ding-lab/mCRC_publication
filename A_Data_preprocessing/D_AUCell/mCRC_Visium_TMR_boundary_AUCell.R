#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(Seurat)
    library(AUCell)
    library(future)
    library(BiocParallel)
    library(qs)
    library(data.table)
    library(tidyverse)
    library(googlesheets4)
})

# Global variables
geneset_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets"
output_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Visium/TMR_genesets/"
file_prefix <- "mCRC_Visium_TMR_boundary"

gs4_deauth()

# Load gene sets
Harada_fib_geneset <- readRDS(file.path(geneset_dir, 'Harada_colon_fibroblast_genesets.rds'))
WU_fib_geneset <- readRDS(file.path(geneset_dir, 'WU_CRC_fibroblast_genesets.rds'))
fib_geneset <- c(Harada_fib_geneset, WU_fib_geneset)

# Function to calculate AUC scores
calculate_auc_scores <- function(gene_set_name, gene_sets, data_matrix) {
    if (is.null(gene_sets[[gene_set_name]])) {
        stop(paste("Gene set", gene_set_name, "is NULL or not found."))
    }
    
    auc_result <- AUCell_run(
        exprMat = data_matrix,
        geneSets = gene_sets[[gene_set_name]],
        BPPARAM = BiocParallel::MulticoreParam(workers = 30)
    )
    
    auc_scores <- auc_result@assays@data$AUC
    auc_vector <- as.vector(auc_scores)
    names(auc_vector) <- colnames(auc_scores)
    return(auc_vector)
}

# Load sample collection info
collection.df <- read_sheet(
    'https://docs.google.com/spreadsheets/d/105cqGuqY_WAeQdvHY-MsI8fBco6lz61OL09_MgZaplE/edit?gid=0#gid=0',
    sheet = "Sample_level"
) %>%
    filter(QCstatus != 'FAILED', cancer_type == 'CRC', !is.na(AnnotationFinal)) %>%
    select(LibraryName, sample_type_tumor, SeuratObject, AnnotationFinal, MorphOutput_FINAL)

head(collection.df)    

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop through samples
for (i in seq_len(nrow(collection.df))) {
    sample_name <- collection.df$LibraryName[i]
    cat(sprintf("\n[%d/%d] Processing sample: %s\n", i, nrow(collection.df), sample_name))
    
    message("Loading annotation file")
    annot_file <- collection.df$MorphOutput_FINAL[i]
    annot <- fread(annot_file, check.names = TRUE) %>%
             as.data.frame() %>%         
             dplyr::rename(TME_boundary = TME.boundary)
    rownames(annot) <- annot$V1
    
    message("Loading Seurat object")
    seurat_obj <- readRDS(collection.df$SeuratObject[i])
    
    message("Adding metadata")
    seurat_obj <- AddMetaData(seurat_obj, annot, col.name = "TME_boundary")
    
    message("Subsetting to annotated spots")
    seurat_obj <- subset(seurat_obj, TME_boundary == 'TME boundary')
    message(sprintf("%d spots with annotation", ncol(seurat_obj)))
    
    message("Getting RNA matrix")
    rna_matrix <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
    
    message("Intersecting genesets")
    Visium_module_geneset <- lapply(fib_geneset, function(gs) intersect(gs, rownames(rna_matrix)))
    
    auc_scores <- list()
    for (gene_set_name in names(Visium_module_geneset)) {
        message(sprintf("  ↳ Calculating AUC for %s", gene_set_name))
        auc_scores[[gene_set_name]] <- calculate_auc_scores(gene_set_name,
                                                            Visium_module_geneset,
                                                            rna_matrix)
    }
    
    auc_matrix <- do.call(cbind, auc_scores)
    auc_matrix_df <- as.data.frame(auc_matrix)
    auc_matrix_df$Spot_ID <- rownames(auc_matrix_df)
    
    output_file <- file.path(output_dir,
        sprintf("%s_%s__TME_boundary_fib_genesets_auc_score.csv", file_prefix, sample_name))
    fwrite(auc_matrix_df, output_file, row.names = FALSE)
    
    message("Done.\n")
    
}

cat("\n All samples processed.\n")