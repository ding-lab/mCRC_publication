#!/usr/bin/env Rscript
#conda activate seurat5_env

library(Seurat)
library(AUCell)
library(future)
library(BiocParallel)
library(qs)
library(data.table)

# Global variables
geneset_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets'
seurat_obj_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds'
# cell_type_csv <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Xenium/Cell_type/All_cell_type/mCRC_N18_Xenium_banky_celltype_metadata_20250214.csv'
output_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/snRNAseq" 
file_prefix <- "mCRC_snRNAseq_epithelial"

# Functions
calculate_auc_scores <- function(gene_set_name, gene_sets, data_matrix) {

	if (is.null(gene_sets[[gene_set_name]])) { 
	stop(paste("Gene set", gene_set_name, "is NULL or not found in gene_sets"))
	}

	auc_result <- AUCell_run(data_matrix, gene_sets[[gene_set_name]], BPPARAM = BiocParallel::MulticoreParam(30))
	auc_scores <- auc_result@assays@data$AUC
	auc_vector <- as.vector(auc_scores)
	names(auc_vector) <- colnames(auc_scores)
	return(auc_vector)
}


# Load geneset
colon_geneset = readRDS(file.path(geneset_dir, 'colon_epithelium_genesets.rds'))
module_geneset = readRDS(file.path(geneset_dir, 'gene_modules.rds'))
hallmark_hypoxia_geneset = readRDS(file.path(geneset_dir, 'Hallmark_hpoxia_genesets.rds'))
mCRC_selected_genesets = readRDS(file.path(geneset_dir, 'mCRC_selected_genesets.rds'))
Sushie_CRC_curated_genesets = readRDS(file.path(geneset_dir, 'Sushie_CRC_curated_genesets.rds'))
module_geneset = c(module_geneset, hallmark_hypoxia_geneset, mCRC_selected_genesets, colon_geneset, Sushie_CRC_curated_genesets)

names(module_geneset) <- gsub(" ", "_", names(module_geneset))
names(module_geneset) <- gsub("-", "_", names(module_geneset))
names(module_geneset)

# Load seurat obj
seurat_obj = readRDS(seurat_obj_path)

# # Load cell types
# cell_tbl <- read.csv(cell_type_csv, row.names = 1)
# seurat_obj <- AddMetaData(seurat_obj, cell_tbl)
# seurat_obj <- subset(seurat_obj, broad_cell_type == 'Tumor')

# Get RNA matrix
tumor_rna_matrix <- seurat_obj[['RNA']]$'counts'

# # Get gene intersection between genesets and Xenium genes
# Xenium_5K_module_geneset <- lapply(module_geneset, function(gene_set) {
#   intersect(gene_set, rownames(tumor_rna_matrix))
# })


# Get AUC score matrix
auc_scores <- list() 
for(gene_set_name in names(module_geneset)) {
	auc_scores[[gene_set_name]] <- calculate_auc_scores(gene_set_name, module_geneset, tumor_rna_matrix) 
}

auc_matrix <- do.call(cbind, auc_scores)
colnames(auc_matrix) <- names(auc_scores)

# Save AUC matrix
auc_matrix_df <- as.data.frame(auc_matrix)
output_file <- file.path(output_dir, paste0(file_prefix, "_genesets_auc_score.csv")) 
fwrite(auc_matrix_df, output_file, row.names=TRUE)

