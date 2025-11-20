#!/usr/bin/env Rscript

library(Seurat)
library(AUCell)
library(future)
library(BiocParallel)
library(qs)
library(data.table)
library(tidyverse)
library(googlesheets4)

# Global variables
geneset_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets"
output_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Visium/TMR_genesets/"
file_prefix <- "mCRC_Visium_TMR"

gs4_deauth()

# Load gene sets
module_geneset <- readRDS(file.path(geneset_dir, 'gene_modules.rds'))
hallmark_hypoxia_geneset <- readRDS(file.path(geneset_dir, 'Hallmark_hpoxia_genesets.rds'))
mCRC_selected_genesets <- readRDS(file.path(geneset_dir, 'mCRC_selected_genesets.rds'))
washu_snrna_geneset <- readRDS(file.path(geneset_dir, 'WASHU_snRNA_CRC_epithelial_genesets.rds'))
Sushie_CRC_curated_genesets = readRDS(file.path(geneset_dir, 'Sushie_CRC_curated_genesets.rds'))
module_geneset = c(module_geneset, hallmark_hypoxia_geneset, mCRC_selected_genesets, washu_snrna_geneset, Sushie_CRC_curated_genesets)


# Clean up geneset names
names(module_geneset) <- gsub(" ", "_", names(module_geneset))
names(module_geneset) <- gsub("-", "_", names(module_geneset))

# Function to calculate AUC scores
calculate_auc_scores <- function(gene_set_name, gene_sets, data_matrix) {
    if (is.null(gene_sets[[gene_set_name]])) { 
        stop(paste("Gene set", gene_set_name, "is NULL or not found in gene_sets"))
    }
    
    auc_result <- AUCell_run(data_matrix, gene_sets[[gene_set_name]], 
                            BPPARAM = BiocParallel::MulticoreParam(30))
    auc_scores <- auc_result@assays@data$AUC
    auc_vector <- as.vector(auc_scores)
    names(auc_vector) <- colnames(auc_scores)
    return(auc_vector)
}

# Load sample collection information
collection.df <- read_sheet('https://docs.google.com/spreadsheets/d/105cqGuqY_WAeQdvHY-MsI8fBco6lz61OL09_MgZaplE/edit?gid=0#gid=0', 
                          sheet = "Sample_level")
collection.df <- collection.df %>% 
    dplyr::filter(QCstatus != 'FAILED') %>% 
    dplyr::filter(cancer_type == 'CRC') %>% 
    dplyr::filter(!is.na(AnnotationFinal)) %>% 
    select(LibraryName, sample_type_tumor, SeuratObject, AnnotationFinal)

# Create directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Process each sample
for(i in 1:nrow(collection.df)) {
    sample_name <- collection.df$LibraryName[i]
    cat(sprintf("\nProcessing sample %d of %d: %s\n", i, nrow(collection.df), sample_name))
    
    # Check if annotation file exists
    if (!file.exists(collection.df$AnnotationFinal[i])) {
        warning(sprintf("Annotation file not found for %s: %s\nSkipping this sample.", 
                       sample_name, collection.df$AnnotationFinal[i]))
        next
    }
    
    tryCatch({
        # Load annotation and Seurat object
        annot <- fread(collection.df$AnnotationFinal[i]) %>% as.data.frame()
        colnames(annot)[2] <- "TMR"
        rownames(annot) <- annot$Barcode
        seurat_obj <- readRDS(collection.df$SeuratObject[i])
        
        # Add TMR annotation
        seurat_obj <- AddMetaData(seurat_obj, annot, col.name = "TMR")
        
        # Subset Seurat object to only include spots with TMR annotations
        spots_with_tmr <- colnames(seurat_obj)[!(seurat_obj$TMR %in% c("", "NA"))]
        seurat_obj <- subset(seurat_obj, cells = spots_with_tmr)
        
        cat(sprintf("Number of spots with TMR annotations: %d\n", ncol(seurat_obj)))
        
        # Get RNA matrix for spots with TMR annotations
        rna_matrix <- seurat_obj[["SCT"]]$data
        
        # Get gene intersection between genesets and Visium genes
        Visium_module_geneset <- lapply(module_geneset, function(gene_set) {
            intersect(gene_set, rownames(rna_matrix))
        })
        
        # Calculate AUC scores
        auc_scores <- list()
        for(gene_set_name in names(Visium_module_geneset)) {
            cat(sprintf("Calculating AUC scores for %s\n", gene_set_name))
            auc_scores[[gene_set_name]] <- calculate_auc_scores(gene_set_name,
                                                              Visium_module_geneset,
                                                              rna_matrix)
        }
        
        # Convert to matrix
        auc_matrix <- do.call(cbind, auc_scores)
        auc_matrix_df <- as.data.frame(auc_matrix)
        
        # Save AUC matrix
        output_file <- file.path(output_dir, 
                               paste0(file_prefix, "_", sample_name, "_genesets_auc_score.csv"))
        fwrite(auc_matrix_df, output_file, row.names=TRUE)
        
        cat(sprintf("Completed processing %s\n", sample_name))
        
    }, error = function(e) {
        warning(sprintf("Error processing sample %s: %s", sample_name, e$message))
    })
}

cat("\nProcessing completed for all samples\n") 