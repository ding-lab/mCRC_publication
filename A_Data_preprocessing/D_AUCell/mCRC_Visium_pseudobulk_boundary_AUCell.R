#!/usr/bin/env Rscript

cat("=== Starting Pseudobulk Boundary AUCell Score Calculation Script ===\n")
cat("Loading required libraries...\n")

library(Seurat)
library(AUCell)
library(future)
library(BiocParallel)
library(data.table)
library(tidyverse)

cat("Libraries loaded successfully\n\n")

# Global variables
cat("Setting up global variables...\n")
geneset_path <- file.path(getwd(), "genesets")
tracking_file <- file.path(getwd(), "/CRC_Visium_tracking.csv")
output_dir <- getwd()
output_file <- file.path(output_dir, "mCRC_Visium_pseudobulk_boundary_genesets_auc_score.csv")

cat(sprintf("Geneset path: %s\n", geneset_path))
cat(sprintf("Tracking file: %s\n", tracking_file))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Output file: %s\n\n", output_file))

# Load gene sets
cat("Step 1: Loading gene sets...\n")
cat(sprintf("  Loading CRC genesets from: %s\n", 
           file.path(geneset_path, "WASHU_snRNA_CRC_epithelial_genesets.rds")))
CRC_genesets <- readRDS(file.path(geneset_path, "WASHU_snRNA_CRC_epithelial_genesets.rds"))
cat(sprintf("  CRC genesets loaded: %d gene sets\n", length(CRC_genesets)))

cat(sprintf("  Loading FIB genesets from: %s\n", 
           file.path(geneset_path, "WASHU_snRNA_CRC_fibroblast_genesets.rds")))
FIB_genesets <- readRDS(file.path(geneset_path, "WASHU_snRNA_CRC_fibroblast_genesets.rds"))
cat(sprintf("  FIB genesets loaded: %d gene sets\n", length(FIB_genesets)))

# Filter CRC_genesets to only include specified gene sets
cat("Step 2: Filtering CRC genesets...\n")
CRC_selected <- c('Canonical_CRC_Intestine', 
                  'Canonical_CRC_Intestine_Proliferation',
                  'Canonical_CRC_Stem',
                  'Canonical_CRC_Stem_Proliferation',
                  'Non_Canonical_CRC_1',
                  'Non_Canonical_CRC_2')

cat(sprintf("  Selected CRC gene sets: %s\n", paste(CRC_selected, collapse=", ")))
CRC_genesets <- CRC_genesets[CRC_selected]
cat(sprintf("  Filtered CRC genesets: %d gene sets\n", length(CRC_genesets)))

# Combine gene sets
cat("Step 3: Combining gene sets...\n")
module_geneset <- c(CRC_genesets, FIB_genesets)
cat(sprintf("  Total combined gene sets: %d\n", length(module_geneset)))

# Clean up geneset names
cat("Step 4: Cleaning geneset names...\n")
names(module_geneset) <- gsub(" ", "_", names(module_geneset))
names(module_geneset) <- gsub("-", "_", names(module_geneset))
cat(sprintf("  Gene set names cleaned\n\n"))

# Function to calculate pseudobulk (mean expression across spots)
create_pseudobulk <- function(expression_matrix, spot_groups) {
    # Remove NA groups
    valid_spots <- !is.na(spot_groups) & spot_groups != "" & spot_groups != "NA"
    expression_matrix <- expression_matrix[, valid_spots, drop=FALSE]
    spot_groups <- spot_groups[valid_spots]
    
    unique_groups <- unique(spot_groups)
    cat(sprintf("    Creating pseudobulk for %d groups: %s\n", 
               length(unique_groups), paste(unique_groups, collapse=", ")))
    
    if (length(unique_groups) == 0) {
        stop("No valid spot groups found")
    }
    
    # Calculate mean expression for each group
    pseudobulk_list <- lapply(unique_groups, function(group) {
        group_spots <- spot_groups == group
        n_spots <- sum(group_spots)
        if (n_spots == 0) {
            return(NULL)
        }
        group_expr <- expression_matrix[, group_spots, drop=FALSE]
        # Calculate mean across spots (handles both single and multiple spots)
        if (n_spots == 1) {
            as.vector(group_expr)
        } else {
            rowMeans(as.matrix(group_expr))
        }
    })
    
    names(pseudobulk_list) <- unique_groups
    pseudobulk_list <- pseudobulk_list[!sapply(pseudobulk_list, is.null)]
    
    if (length(pseudobulk_list) == 0) {
        stop("No valid pseudobulk groups created")
    }
    
    # Convert to matrix (genes x groups)
    pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
    rownames(pseudobulk_matrix) <- rownames(expression_matrix)
    
    cat(sprintf("    Pseudobulk matrix created: %d genes x %d groups\n", 
               nrow(pseudobulk_matrix), ncol(pseudobulk_matrix)))
    return(pseudobulk_matrix)
}

# Function to calculate AUC scores on pseudobulk
calculate_auc_scores_pseudobulk <- function(gene_set_name, gene_sets, data_matrix) {
    if (is.null(gene_sets[[gene_set_name]])) { 
        stop(paste("Gene set", gene_set_name, "is NULL or not found in gene_sets"))
    }
    
    cat(sprintf("    [AUCell] Starting calculation for %s...\n", gene_set_name))
    cat(sprintf("    [AUCell] Gene set has %d genes\n", length(gene_sets[[gene_set_name]])))
    cat(sprintf("    [AUCell] Running AUCell with 30 cores...\n"))
    
    auc_result <- AUCell_run(data_matrix, gene_sets[[gene_set_name]], 
                            BPPARAM = BiocParallel::MulticoreParam(30))
    
    cat(sprintf("    [AUCell] AUCell calculation completed\n"))
    auc_scores <- auc_result@assays@data$AUC
    auc_vector <- as.vector(auc_scores)
    names(auc_vector) <- colnames(auc_scores)
    
    cat(sprintf("    [AUCell] Extracted AUC scores for %d pseudobulk samples\n", length(auc_vector)))
    return(auc_vector)
}

# Load sample collection information from CSV
cat("Step 5: Loading sample tracking information...\n")
cat(sprintf("  Reading CSV file: %s\n", tracking_file))
collection.df <- read_csv(tracking_file, show_col_types = FALSE)
cat(sprintf("  Loaded %d samples from tracking file\n", nrow(collection.df)))

# Select required columns
cat("Step 6: Selecting required columns...\n")
collection.df <- collection.df %>% 
    select(LibraryName, SeuratObject, MorphOutput_FINAL)
cat(sprintf("  Selected columns: LibraryName, SeuratObject, MorphOutput_FINAL\n"))

# Create directory if it doesn't exist
cat("Step 7: Creating output directory...\n")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  Output directory ready: %s\n\n", output_dir))

# Process each sample
cat("Step 8: Starting sample processing loop...\n")
cat(sprintf("Total samples to process: %d\n\n", nrow(collection.df)))

# Store results for all samples
all_results <- list()

for(i in 1:nrow(collection.df)) {
    sample_name <- collection.df$LibraryName[i]
    cat(sprintf("\n========================================\n"))
    cat(sprintf("Processing sample %d of %d: %s\n", i, nrow(collection.df), sample_name))
    cat(sprintf("========================================\n"))
    
    # Check if Seurat object exists
    cat(sprintf("[Sample %s] Step 8.1: Checking Seurat object path...\n", sample_name))
    if (is.na(collection.df$SeuratObject[i])) {
        cat(sprintf("[Sample %s] WARNING: SeuratObject path is NA\n", sample_name))
        warning(sprintf("Seurat object not found for %s: %s\nSkipping this sample.", 
                       sample_name, collection.df$SeuratObject[i]))
        next
    }
    
    seurat_path <- collection.df$SeuratObject[i]
    if (!file.exists(seurat_path)) {
        cat(sprintf("[Sample %s] WARNING: Seurat object file does not exist: %s\n", 
                   sample_name, seurat_path))
        warning(sprintf("Seurat object not found for %s: %s\nSkipping this sample.", 
                       sample_name, seurat_path))
        next
    }
    
    # Check if annotation file exists
    if (is.na(collection.df$MorphOutput_FINAL[i]) || 
        !file.exists(collection.df$MorphOutput_FINAL[i])) {
        cat(sprintf("[Sample %s] WARNING: Annotation file not found, skipping\n", sample_name))
        next
    }
    
    tryCatch({
        # Load Seurat object
        cat(sprintf("[Sample %s] Step 8.2: Loading Seurat object...\n", sample_name))
        cat(sprintf("[Sample %s]   File path: %s\n", sample_name, seurat_path))
        seurat_obj <- readRDS(seurat_path)
        cat(sprintf("[Sample %s]   Seurat object loaded successfully\n", sample_name))
        
        # Load annotation file
        cat(sprintf("[Sample %s] Step 8.3: Loading annotation file...\n", sample_name))
        annot <- fread(collection.df$MorphOutput_FINAL[i]) %>% as.data.frame()
        cat(sprintf("[Sample %s]   Annotation file loaded: %d rows, %d columns\n", 
                   sample_name, nrow(annot), ncol(annot)))
        cat(sprintf("[Sample %s]   Column names: %s\n", 
                   sample_name, paste(colnames(annot), collapse=", ")))
        
        # Check for required columns
        if (!"Barcode" %in% colnames(annot)) {
            cat(sprintf("[Sample %s]   ERROR: 'Barcode' column not found\n", sample_name))
            next
        }
        
        if (!"Tumor boundary" %in% colnames(annot) && !"TME boundary" %in% colnames(annot)) {
            cat(sprintf("[Sample %s]   WARNING: Neither 'Tumor boundary' nor 'TME boundary' columns found\n", sample_name))
            next
        }
        
        # Get RNA matrix
        cat(sprintf("[Sample %s] Step 8.4: Extracting expression matrix...\n", sample_name))
        if ("SCT" %in% names(seurat_obj@assays)) {
            cat(sprintf("[Sample %s]   Using SCT assay\n", sample_name))
            rna_matrix <- seurat_obj[["SCT"]]$data
        } else if ("RNA" %in% names(seurat_obj@assays)) {
            cat(sprintf("[Sample %s]   SCT assay not found, using RNA assay\n", sample_name))
            rna_matrix <- seurat_obj[["RNA"]]$data
        } else {
            stop(sprintf("Neither SCT nor RNA assay found. Available assays: %s", 
                        paste(names(seurat_obj@assays), collapse=", ")))
        }
        
        cat(sprintf("[Sample %s]   RNA matrix dimensions: %d genes x %d spots\n", 
                   sample_name, nrow(rna_matrix), ncol(rna_matrix)))
        
        # Match barcodes between annotation and Seurat object
        cat(sprintf("[Sample %s] Step 8.5: Matching barcodes...\n", sample_name))
        common_barcodes <- intersect(annot$Barcode, colnames(rna_matrix))
        cat(sprintf("[Sample %s]   Common barcodes: %d\n", sample_name, length(common_barcodes)))
        
        if (length(common_barcodes) == 0) {
            cat(sprintf("[Sample %s]   ERROR: No common barcodes found\n", sample_name))
            next
        }
        
        # Subset annotation and expression matrix
        annot_matched <- annot[match(common_barcodes, annot$Barcode), ]
        rna_matrix_matched <- rna_matrix[, common_barcodes, drop=FALSE]
        
        # Create boundary groups
        cat(sprintf("[Sample %s] Step 8.6: Creating boundary groups...\n", sample_name))
        boundary_groups <- rep(NA_character_, length(common_barcodes))
        
        if ("Tumor boundary" %in% colnames(annot_matched)) {
            tumor_boundary_spots <- !is.na(annot_matched[["Tumor boundary"]]) & 
                                   annot_matched[["Tumor boundary"]] != "" &
                                   annot_matched[["Tumor boundary"]] == "Tumor boundary"
            boundary_groups[tumor_boundary_spots] <- "Tumor_boundary"
            cat(sprintf("[Sample %s]   Tumor boundary spots: %d\n", 
                       sample_name, sum(tumor_boundary_spots)))
        }
        
        if ("TME boundary" %in% colnames(annot_matched)) {
            tme_boundary_spots <- !is.na(annot_matched[["TME boundary"]]) & 
                                 annot_matched[["TME boundary"]] != "" &
                                 annot_matched[["TME boundary"]] == "TME boundary"
            boundary_groups[tme_boundary_spots] <- "TME_boundary"
            cat(sprintf("[Sample %s]   TME boundary spots: %d\n", 
                       sample_name, sum(tme_boundary_spots)))
        }
        
        # Check if we have any boundary spots
        valid_groups <- boundary_groups[!is.na(boundary_groups)]
        if (length(valid_groups) == 0) {
            cat(sprintf("[Sample %s]   WARNING: No boundary spots found, skipping\n", sample_name))
            next
        }
        
        cat(sprintf("[Sample %s]   Total boundary spots: %d\n", 
                   sample_name, length(valid_groups)))
        
        # Create pseudobulk
        cat(sprintf("[Sample %s] Step 8.7: Creating pseudobulk...\n", sample_name))
        pseudobulk_matrix <- create_pseudobulk(rna_matrix_matched, boundary_groups)
        
        if (ncol(pseudobulk_matrix) == 0) {
            cat(sprintf("[Sample %s]   WARNING: No pseudobulk groups created, skipping\n", sample_name))
            next
        }
        
        # Get gene intersection between genesets and expression genes
        cat(sprintf("[Sample %s] Step 8.8: Finding gene intersections...\n", sample_name))
        Visium_module_geneset <- lapply(module_geneset, function(gene_set) {
            intersect(gene_set, rownames(pseudobulk_matrix))
        })
        
        # Filter out empty gene sets
        Visium_module_geneset <- Visium_module_geneset[sapply(Visium_module_geneset, length) > 0]
        cat(sprintf("[Sample %s]   Number of gene sets with overlapping genes: %d\n", 
                   sample_name, length(Visium_module_geneset)))
        
        if (length(Visium_module_geneset) == 0) {
            cat(sprintf("[Sample %s]   WARNING: No gene sets have overlapping genes\n", sample_name))
            next
        }
        
        # Calculate AUC scores for each pseudobulk group
        cat(sprintf("[Sample %s] Step 8.9: Calculating AUC scores for pseudobulk...\n", sample_name))
        auc_results_sample <- list()
        
        for(gene_set_name in names(Visium_module_geneset)) {
            cat(sprintf("[Sample %s]   Processing gene set: %s (%d genes)\n", 
                       sample_name, gene_set_name, length(Visium_module_geneset[[gene_set_name]])))
            auc_scores <- calculate_auc_scores_pseudobulk(gene_set_name,
                                                          Visium_module_geneset,
                                                          pseudobulk_matrix)
            auc_results_sample[[gene_set_name]] <- auc_scores
        }
        
        # Convert to data frame - each pseudobulk group becomes a row
        cat(sprintf("[Sample %s] Step 8.10: Formatting results...\n", sample_name))
        
        # Get all pseudobulk group names (columns of pseudobulk_matrix)
        pseudobulk_groups <- colnames(pseudobulk_matrix)
        
        # Get all gene set names (use the full module_geneset to ensure consistency)
        all_gene_set_names <- names(module_geneset)
        
        # For each pseudobulk group, extract scores across all gene sets
        for(group_name in pseudobulk_groups) {
            row_name <- paste0(sample_name, "_", group_name)
            
            # Extract AUC scores for this group across all gene sets
            # Ensure we have all gene sets, even if they weren't calculated for this sample
            row_data <- sapply(all_gene_set_names, function(gs_name) {
                if (gs_name %in% names(auc_results_sample)) {
                    auc_scores <- auc_results_sample[[gs_name]]
                    if (group_name %in% names(auc_scores)) {
                        return(auc_scores[group_name])
                    }
                }
                return(NA)
            })
            
            # Create data frame row with all gene sets as columns
            row_df <- data.frame(t(row_data), row.names=row_name, check.names=FALSE)
            colnames(row_df) <- all_gene_set_names
            all_results[[row_name]] <- row_df
            
            cat(sprintf("[Sample %s]   Added row: %s with %d gene sets\n", 
                       sample_name, row_name, length(all_gene_set_names)))
        }
        
        cat(sprintf("[Sample %s] ✓ Completed processing %s\n", sample_name, sample_name))
        
    }, error = function(e) {
        cat(sprintf("\n[Sample %s] ✗ ERROR occurred!\n", sample_name), file=stderr())
        cat(sprintf("[Sample %s] Error message: %s\n", sample_name, e$message), file=stderr())
        warning(sprintf("Error processing sample %s: %s", sample_name, e$message))
    })
}

# Combine all results
cat("\n========================================\n")
cat("Step 9: Combining all results...\n")
if (length(all_results) == 0) {
    cat("WARNING: No results to combine!\n")
} else {
    cat(sprintf("Combining %d result rows...\n", length(all_results)))
    
    # Check column names consistency
    all_colnames <- lapply(all_results, colnames)
    unique_colnames <- unique(all_colnames)
    if (length(unique_colnames) > 1) {
        cat("WARNING: Different column names found across results\n")
        cat("Ensuring all data frames have the same columns...\n")
        
        # Get all unique column names
        all_unique_cols <- unique(unlist(all_colnames))
        
        # Ensure all data frames have the same columns
        all_results <- lapply(all_results, function(df) {
            missing_cols <- setdiff(all_unique_cols, colnames(df))
            if (length(missing_cols) > 0) {
                for(col in missing_cols) {
                    df[[col]] <- NA
                }
            }
            # Reorder columns to match
            df <- df[, all_unique_cols, drop=FALSE]
            return(df)
        })
    }
    
    final_df <- do.call(rbind, all_results)
    
    cat(sprintf("Final table dimensions: %d rows x %d columns\n", 
               nrow(final_df), ncol(final_df)))
    cat(sprintf("Row names: %s\n", paste(rownames(final_df), collapse=", ")))
    cat(sprintf("Column names: %s\n", paste(colnames(final_df), collapse=", ")))
    
    # Save results
    cat("Step 10: Saving results...\n")
    cat(sprintf("Saving to: %s\n", output_file))
    fwrite(final_df, output_file, row.names=TRUE)
    
    if (file.exists(output_file)) {
        file_size <- file.info(output_file)$size
        cat(sprintf("File saved successfully! Size: %d bytes\n", file_size))
    } else {
        cat("WARNING: File was not created!\n")
    }
}

cat("\n========================================\n")
cat("Processing completed for all samples\n")
cat("========================================\n")

