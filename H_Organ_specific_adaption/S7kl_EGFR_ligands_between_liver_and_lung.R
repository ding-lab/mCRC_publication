library(Seurat)
library(tidyverse)
library(qs)
library(ggpubr)
library(future)
options(future.globals.maxSize = 30 * 1024^3)
out_dir = getwd()

PSI50_cells <- qread('/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_LiLu_PSI50_stromal.qs')
PSI50_cells

genes_present <- c("EGF", "HBEGF", "TGFA",  "EREG")

## ================================================================
##  Build per-cell expression + metadata table
## ================================================================
expr_mat <- GetAssayData(
  PSI50_cells,
  assay = "Xenium",
  layer = "data"
)[genes_present, , drop = FALSE]

expr_df <- expr_mat %>%
  t() %>%                     # cells x genes
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(
    PSI50_cells@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::select(cell, Organ, orig.ident),
    by = "cell"
  ) %>%
  # keep only Liver/Lung and non-missing Organ
  dplyr::filter(!is.na(Organ), Organ %in% c("Liver", "Lung"))

## ================================================================
##  Aggregate to sample-level expression (mean per sample per gene)
## ================================================================
sample_expr <- expr_df %>%
  tidyr::pivot_longer(
    cols = all_of(genes_present),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  dplyr::group_by(orig.ident, Organ, gene) %>%
  dplyr::summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    .groups = "drop"
  )

## ================================================================
##  Filter out unwanted samples
## ================================================================
remove_samples <- c(
  "HT472C1-M1FP1U1",   # fragmented
  "HT413C1-Th1K2A1U2"  # MSI-High
)

sample_expr_filt <- sample_expr %>%
  dplyr::filter(!orig.ident %in% remove_samples)

cat("Samples kept after filtering:\n")
print(unique(sample_expr_filt$orig.ident))

## ================================================================
##  Common color scheme
## ================================================================
organ_colors <- c(
  "Liver" = "brown",
  "Lung"  = "steelblue"
)

## ================================================================
##  Figure 1: Liver vs Lung comparison (sample-level)
##  - Each point = one sample
##  - Box shows distribution across samples by Organ
## ================================================================
p_fig1 <- ggplot(
  sample_expr_filt,
  aes(x = Organ, y = mean_expr, fill = Organ)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black") +
  geom_jitter(
    aes(color = Organ),
    width = 0.15,
    size = 2
  ) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = organ_colors) +
  scale_color_manual(values = organ_colors) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10)
  ) +
  ylab("Mean normalized expression per sample")


pdf(file=file.path(out_dir, 'S7K_EGF_level_between_liver_lung.pdf'), width=11, height=4) 
p_fig1
dev.off()

## ================================================================
##  Figure 2: Sample-wise heatmap
##  - Rows = samples
##  - Columns = EGF-family ligands
##  - Values = mean expression (Z-scored per gene)
## ================================================================

# Build matrix: genes x samples
heatmap_mat <- sample_expr_filt %>%
  dplyr::select(orig.ident, gene, mean_expr) %>%
  tidyr::pivot_wider(
    names_from = orig.ident, 
    values_from = mean_expr
  ) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Make sure gene rows order matches your gene list
heatmap_mat <- heatmap_mat[intersect(genes, rownames(heatmap_mat)), , drop = FALSE]

# Replace NAs with zero or choose another strategy
heatmap_mat[is.na(heatmap_mat)] <- 0

# Z-score by gene (across samples)
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

# ----- Organ annotation (top annotation per sample / per column) -----
sample_anno <- sample_expr_filt %>%
  dplyr::distinct(orig.ident, Organ) %>%
  tibble::column_to_rownames("orig.ident")

# match order
sample_anno <- sample_anno[colnames(heatmap_mat_scaled), , drop = FALSE]

ha_top <- HeatmapAnnotation(
  Organ = sample_anno$Organ,
  col = list(Organ = organ_colors),
  annotation_name_side = "left"
)

# ----- draw heatmap -----
ht_fig2 <- Heatmap(
  heatmap_mat_scaled,
  name = "Z-score",
  top_annotation = ha_top,
  cluster_columns = TRUE,     # clusters samples
  cluster_rows = TRUE,        # clusters genes
  show_row_names = TRUE,      
  show_column_names = TRUE,   
  column_names_rot = 90,      # readable sample IDs
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 7)
)

pdf(file=file.path(out_dir, 'S7L_EGF_level_between_liver_lung_heatmap.pdf'), width=11, height=4) 
ht_fig2
dev.off()