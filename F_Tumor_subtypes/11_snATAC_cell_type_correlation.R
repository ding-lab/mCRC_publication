library(Seurat)
library(Signac)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

set.seed(1234)

out_dir = getwd()

## Run differential accessible ATAC peaks
# rds_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/ATAC/tumor/"
# tumor_obj <- readRDS(file.path(rds_path, "chromvar", 
#                                "41_merged_normalized_mCRC_Combo_clean_tumor.chromvar.rds"))

# Idents(tumor_obj) <- "cell_type_xenium"
# DefaultAssay(tumor_obj) <- "ATAC_merged"

# metadata_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/"
# metadata <- read.csv(file.path(metadata_path, 
#                                "mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv"))

# rownames(metadata) <- metadata$X
# metadata$X <- NULL

# tumor_obj <- AddMetaData(tumor_obj, metadata)

# dap <- FindAllMarkers(
#   object = tumor_obj,
#   logfc.threshold = 0.2,
#   min.pct = 0.2,
#   pseudocount.use = 1,
#   only.pos = FALSE,
#   test.use = 'LR',
#   latent.vars = 'nCount_ATAC_merged'
# ) %>%  mutate(peak=rownames(.))

# Load dap data, so don't need to run it again
dap_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_3/supp/snATAC_analysis'
dap = read.csv(file.path(dap_dir, 'mCRC_tumor_dap2.csv'))


interesting.peaks.tb <- dap %>% 
    dplyr::filter(abs(avg_log2FC) > 0.25 & p_val_adj < 0.05) 
interesting.peaks <- interesting.peaks.tb %>% pull(peak)
interesting.peaks <- unique(interesting.peaks)
length(interesting.peaks)

## Get average peaks by cell types
# avg_obj <- AverageExpression(
#   tumor_obj,
#   assays = "ATAC_merged",
#   group.by = "cell_type_xenium",
#   return.seurat = TRUE,
#   layer = "counts"
# )

# # Add QC metrics for regression
# avg_obj$nCount_ATAC_merged <- Matrix::colSums(
#   LayerData(avg_obj[["ATAC_merged"]], layer = "counts")
# )

# avg_obj$nFeature_ATAC_merged <- Matrix::colSums(
#   LayerData(avg_obj[["ATAC_merged"]], layer = "counts") > 0
# )

# avg_obj <- ScaleData(
#   avg_obj,
#   assay = "ATAC_merged",
#   vars.to.regress = "nCount_ATAC_merged",
#   verbose = FALSE
# )

## Load avg_obj, so don't need to run it again
data_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects'
avg_obj = readRDS(file.path(data_dir, "41_merged_normalized_mCRC_Combo_clean_tumor.average.scaled.rds"))


available_peaks <- rownames(avg_obj[["ATAC_merged"]])

interesting.peaks <- interesting.peaks[interesting.peaks %in% available_peaks]

cat("Number of selected peaks present in ATAC assay:", length(interesting.peaks), "\n")

if (length(interesting.peaks) < 2) {
  stop("Not enough peaks found in ATAC assay for correlation analysis.")
}

mat_scaled_all <- GetAssayData(
  avg_obj,
  assay = "ATAC_merged",
  slot  = "scale.data"
)

mat_scaled_sub <- mat_scaled_all[interesting.peaks, , drop = FALSE]
mat_scaled_sub <- as.matrix(mat_scaled_sub)
cor_mat <- cor(mat_scaled_sub, method = "pearson")

col_fun_scaled <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

p1 <- Heatmap(
  cor_mat,
  name = "Correlation",
  col = col_fun_scaled,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(legend_height = unit(6, "cm"))
)

pdf(file.path(out_dir, "ATAC_celltype_scaled_correlation_heatmap_significant_daps.pdf"), width = 6, height = 5)
draw(p1)
dev.off()

# -----------------------------
# PCA on correlation matrix
# -----------------------------

# cor_scaled is your original correlation matrix (cell type × cell type)

# Convert correlation to a PCA-ready matrix
# PCA expects a covariance-like matrix; correlation is acceptable
pca <- prcomp(cor_mat, center = TRUE, scale. = FALSE)

# Pull PCA coordinates
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  CellType = rownames(pca$x)
)

# Variance explained
var_expl <- pca$sdev^2 / sum(pca$sdev^2)
pc1_lab <- paste0("PC1 (", round(var_expl[1] * 100, 1), "%)")
pc2_lab <- paste0("PC2 (", round(var_expl[2] * 100, 1), "%)")

# Highlight Canonical-CRC-Stem
pca_df$highlight <- ifelse(pca_df$CellType == "Canonical-CRC-Stem", "yes", "no")

library(ggplot2)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = highlight), size = 5) +
  geom_text(aes(label = CellType), vjust = -0.8, size = 4) +
  scale_color_manual(values = c("yes" = "red", "no" = "grey40")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA of ATAC Cell-Type Correlation Matrix",
    x = pc1_lab,
    y = pc2_lab
  ) +
  theme(legend.position = "none")

pdf(file.path(out_dir, "ATAC_celltype_scaled_pca_significant_daps.pdf"), width = 5, height = 5)
p_pca
dev.off()