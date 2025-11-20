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
seurat_obj_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_5K/subset_objects/mCRC_Xenium_N26_stroma_harmonized_mesenchymal_layersjoined.rds'
#cell_type_csv <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Xenium/Cell_type/All_cell_type/mCRC_N21_Xenium_banky_celltype_metadata_20250326.csv'
output_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Xenium/Mesenchymal_genesets/" 
file_prefix <- "mCRC_Xenium5K_N26_mesenchymal"

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
Harada_fib_geneset = readRDS(file.path(geneset_dir, 'Harada_colon_fibroblast_genesets.rds'))
Washu_fib_geneset <- readRDS(file.path(geneset_dir, 'WU_CRC_fibroblast_genesets.rds'))
fib_geneset = c(Harada_fib_geneset, Washu_fib_geneset)



# Load seurat obj
seurat_obj = readRDS(seurat_obj_path)

# Load cell types
#cell_tbl <- read.csv(cell_type_csv, row.names = 1)
#seurat_obj <- AddMetaData(seurat_obj, cell_tbl)
#seurat_obj <- subset(seurat_obj, broad_cell_type == 'Tumor')

# Get Xenium matrix
mesenchymal_xenium_matrix <- seurat_obj[['Xenium']]$'counts'

# Get gene intersection between genesets and Xenium genes
Xenium_5K_module_geneset <- lapply(fib_geneset, function(gene_set) {
  intersect(gene_set, rownames(mesenchymal_xenium_matrix))
})


# Get AUC score matrix
auc_scores <- list() 
for(gene_set_name in names(Xenium_5K_module_geneset)) {
	auc_scores[[gene_set_name]] <- calculate_auc_scores(gene_set_name, Xenium_5K_module_geneset, mesenchymal_xenium_matrix) 
}

auc_matrix <- do.call(cbind, auc_scores)
colnames(auc_matrix) <- names(auc_scores)

# Save AUC matrix
auc_matrix_df <- as.data.frame(auc_matrix)
output_file <- file.path(output_dir, paste0(file_prefix, "_genesets_auc_score.csv")) 
fwrite(auc_matrix_df, output_file, row.names=TRUE)

