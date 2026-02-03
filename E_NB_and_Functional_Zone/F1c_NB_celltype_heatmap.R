library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(qs)
library(circlize)
library(RColorBrewer)
library(googlesheets4)
library(ggpubr)

out_dir <- getwd()

#------------------------------------------------------------
# Load metadata
#------------------------------------------------------------
metadata_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv"

meta_df <- readr::read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"

# Make rownames = barcode (ensure it's plain data.frame)
meta_df <- meta_df |>
  as.data.frame() |>
  column_to_rownames("barcode")

#------------------------------------------------------------
# Filter to desired BANKSY clusters & broad cell types,
# and assign neighborhood labels
#------------------------------------------------------------
meta_df4 <- meta_df |>
  dplyr::filter(
    banksy_leiden_0.7 %in% 0:24,
    Broad_cell_type2 %in% c(
      "Tumor", "Intestinal_epithelium", 
      "Cholangiocyte", "Hepatocyte", "Breast_duct",
      "Aveolar_epithelium", "Necrosis",
      "Myeloid cells", "T_NK_cell", "B_cell", "Plasma_cell",
      "Fibroblast", "Endothelial_cell", "Pericyte", "Smooth_muscle"
    )
  ) |>
  dplyr::mutate(
    neighborhoods = dplyr::case_when(
      banksy_leiden_0.7 %in% c(5, 6)                    ~ "NB_tumor_core",
      banksy_leiden_0.7 %in% c(2, 17)                  ~ "NB_tumor_stroma_interface",
      banksy_leiden_0.7 %in% c(0, 11, 13, 16, 18, 19,
                               20, 24)                 ~ "NB_tumor_body",
      banksy_leiden_0.7 == 3                           ~ "NB_perivascular_stroma",
      banksy_leiden_0.7 == 9                           ~ "NB_desmoplastic_stroma",
      banksy_leiden_0.7 == 15                          ~ "NB_smooth_muscle_stroma",
      banksy_leiden_0.7 == 4                           ~ "NB_myelocyte_enriched_stroma",
      banksy_leiden_0.7 == 1                           ~ "NB_lymphocyte_enriched_stroma",
      banksy_leiden_0.7 == 7                           ~ "NB_plasma_cell_enriched_stroma",
      banksy_leiden_0.7 == 22                          ~ "NB_germinal_center",
      banksy_leiden_0.7 == 10                          ~ "NB_liver",
      banksy_leiden_0.7 == 21                          ~ "NB_ductal_epithelium",
      banksy_leiden_0.7 == 14                          ~ "NB_colon_epithelium",
      banksy_leiden_0.7 %in% c(8, 12, 23)              ~ "NB_necrosis",
      TRUE                                             ~ NA_character_
    )
  ) |>
  dplyr::filter(!is.na(neighborhoods))   # drop any unassigned clusters just in case

#------------------------------------------------------------
# Prepare 2 x 2 contingency tables per NB × cell type
#------------------------------------------------------------
neighborhoods <- sort(unique(meta_df4$neighborhoods))
cell_types   <- sort(unique(meta_df4$Broad_cell_type2))

odds_matrix <- matrix(
  NA_real_,
  nrow = length(neighborhoods),
  ncol = length(cell_types),
  dimnames = list(neighborhoods, cell_types)
)

pval_matrix <- matrix(
  NA_real_,
  nrow = length(neighborhoods),
  ncol = length(cell_types),
  dimnames = list(neighborhoods, cell_types)
)

for (nb in neighborhoods) {
  for (ct in cell_types) {
    
    # a = in nb AND of type ct
    a <- sum(meta_df4$neighborhoods == nb & meta_df4$Broad_cell_type2 == ct)
    # b = in nb AND NOT ct
    b <- sum(meta_df4$neighborhoods == nb & meta_df4$Broad_cell_type2 != ct)
    # c = NOT nb AND ct
    c <- sum(meta_df4$neighborhoods != nb & meta_df4$Broad_cell_type2 == ct)
    # d = NOT nb AND NOT ct
    d <- sum(meta_df4$neighborhoods != nb & meta_df4$Broad_cell_type2 != ct)
    
    # Safety: skip if degenerate table
    if ((a + b + c + d) == 0) {
      odds_matrix[nb, ct] <- NA_real_
      pval_matrix[nb, ct] <- NA_real_
      next
    }
    
    tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    # Fisher odds ratio (enrichment)
    fisher <- fisher.test(tbl)
    odds_matrix[nb, ct] <- fisher$estimate
    
    # Hypergeometric p-value (enrichment of ct in nb)
    m <- a + c  # total ct cells
    n <- b + d  # total non-ct cells
    k <- a + b  # total cells in nb
    q <- ifelse(a > 0, a - 1, 0)
    
    pval_matrix[nb, ct] <- phyper(q, m, n, k, lower.tail = FALSE)
  }
}

#------------------------------------------------------------
# Custom row/column order (keep only those present)
#------------------------------------------------------------
row_order3 <- c(
  "NB_necrosis", "NB_tumor_core", "NB_tumor_body", "NB_tumor_stroma_interface", 
  "NB_desmoplastic_stroma", "NB_perivascular_stroma",
  "NB_myelocyte_enriched_stroma", "NB_lymphocyte_enriched_stroma",
  "NB_plasma_cell_enriched_stroma", 
  "NB_smooth_muscle_stroma", "NB_germinal_center", "NB_colon_epithelium",
  "NB_liver", "NB_ductal_epithelium"
)

col_order3 <- c(
  "Necrosis", "Tumor",
  "Fibroblast", "Endothelial_cell", "Pericyte", "Smooth_muscle",
  "Myeloid cells", "T_NK_cell", "B_cell", "Plasma_cell",
  "Intestinal_epithelium", "Hepatocyte", "Cholangiocyte",
  "Breast_duct", "Aveolar_epithelium"
)

row_order3 <- row_order3[row_order3 %in% rownames(odds_matrix)]
col_order3 <- col_order3[col_order3 %in% colnames(odds_matrix)]

#------------------------------------------------------------
# Log2 transform odds ratios, handle infinities
#------------------------------------------------------------
log2_odds_matrix <- log2(odds_matrix)

finite_vals <- log2_odds_matrix[is.finite(log2_odds_matrix)]
min_val <- min(finite_vals, na.rm = TRUE)
max_val <- max(finite_vals, na.rm = TRUE)

# Replace -Inf with min finite, +Inf with max finite
log2_odds_matrix[is.infinite(log2_odds_matrix) & log2_odds_matrix < 0] <- min_val
log2_odds_matrix[is.infinite(log2_odds_matrix) & log2_odds_matrix > 0] <- max_val

#------------------------------------------------------------
# FDR adjust p-values & mark significant cells
#------------------------------------------------------------
pval_vector <- as.vector(pval_matrix)
padj_vector <- p.adjust(pval_vector, method = "fdr")
padj_matrix <- matrix(
  padj_vector,
  nrow = nrow(pval_matrix),
  ncol = ncol(pval_matrix),
  dimnames = dimnames(pval_matrix)
)

star_matrix <- ifelse(padj_matrix < 0.05, "*", "")

#------------------------------------------------------------
# Heatmap
#------------------------------------------------------------
palette <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))

# You can adjust breaks to your actual range if you want
col_fun <- circlize::colorRamp2(
  breaks = c(-3, 0, 3),
  colors = palette[c(1, 6, 10)]
)

ht <- Heatmap(
  log2_odds_matrix[row_order3, col_order3, drop = FALSE],
  name = "log2OR",
  col = col_fun,
  na_col = "grey95",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(star_matrix[row_order3[i], col_order3[j]])) {
      grid.text(star_matrix[row_order3[i], col_order3[j]], x, y,
                gp = gpar(fontsize = 10))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_title = "Neighborhoods",
  column_title = "Broad cell types"
)

pdf(file.path(out_dir, "mCRC_N26_neighborhoods_cell_types_simplified.pdf"),
    width = 9, height = 8)
draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()

#------------------------------------------------------------
# NB UMAP
#------------------------------------------------------------
mCRC_subset <- qread(
  "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_with_banksy_umap.qs"
)

# Reproducible sampling; handle case where < 500k cells
set.seed(123)
n_cells <- ncol(mCRC_subset)
n_keep  <- min(500000, n_cells)
subset_cells <- sample(colnames(mCRC_subset), n_keep)

mCRC_500k <- subset(mCRC_subset, cells = subset_cells)

# Drop cells with missing neighborhood annotation
mCRC_500k <- subset(mCRC_500k, subset = !is.na(neighborhoods))

NB_colors <- c(
  "NB_colon_epithelium"          = "gold2",
  "NB_tumor_core"                = "hotpink4",
  "NB_tumor_body"                = "violetred",
  "NB_tumor_stroma_interface"    = "orange1",
  "NB_desmoplastic_stroma"       = "darkgreen",
  "NB_perivascular_stroma"       = "yellowgreen",
  "NB_smooth_muscle_stroma"      = "lightseagreen",
  "NB_myelocyte_enriched_stroma" = "#be0032",
  "NB_lymphocyte_enriched_stroma"= "#dcd300",
  "NB_plasma_cell_enriched_stroma"= "#0067a5",
  "NB_germinal_center"           = "yellow4",
  "NB_liver"                     = "#882d17",
  "NB_ductal_epithelium"         = "#2b3d26",
  "NB_necrosis"                  = "#604e97"
)

mCRC_500k$neighborhoods <- factor(
  mCRC_500k$neighborhoods,
  levels = names(NB_colors)
)

umap1 <- DimPlot(
  mCRC_500k,
  reduction = "banksy_umap",
  group.by  = "neighborhoods",
  cols      = NB_colors
)

pdf(file.path(out_dir, "mCRC_N26_500K_dimplot.pdf"),
    width = 6, height = 4)
print(umap1)
dev.off()
