# Merged from former 00_00_setup through 06_geneset_auc_score.
# Loads epithelial object with cell typing already applied (NE_included); no FindSubCluster / re-annotation.

# ---- Libraries ----
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(viridisLite)
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(circlize)

support_r <- "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/support_scripts/jupyter_support_functions.R"
if (file.exists(support_r)) {
  source(support_r)
}

# ---- Paths ----
output_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Neuroendocrine_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumor_rds_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE_included.rds"
tumor_obj <- readRDS(tumor_rds_path)

# Tumor epithelial subtypes used for subset dotplot / AUC (same order as notebook)
tumor_epi_types <- c(
  "CMETS",
  "Neuroendocrine-like tumor",
  "APCDD1_CRC",
  "Canonical_CRC_Stem",
  "Canonical_CRC_Stem_Proliferation",
  "Canonical_CRC_Intestine_Proliferation",
  "Canonical_CRC_Intestine"
)

# =============================================================================
# UMAP — epithelial cell types
# =============================================================================
epi.cell.color <- c(
  "APCDD1_CRC" = "indianred",
  "CMETS" = "#be0032",
  "Neuroendocrine-like tumor" = "magenta4",
  "Canonical_CRC_Intestine_Proliferation" = "olivedrab4",
  "Canonical_CRC_Stem_Proliferation" = "darkolivegreen",
  "Canonical_CRC_Intestine" = "seagreen",
  "Canonical_CRC_Stem" = "darkseagreen",
  "Stem cells" = "yellow3",
  "Enterocytes" = "goldenrod2",
  "Transit-amplifying cells" = "goldenrod4",
  "Goblet cells" = "gold2",
  "Tuft cells" = "orange3",
  "Neuroendocrine cells" = "lightyellow3"
)

options(repr.plot.width = 10, repr.plot.height = 8)
p1 <- DimPlot(tumor_obj, reduction = "epithelial_umap.scvi", group.by = "cell_type_all3", label = FALSE) +
  scale_color_manual(values = epi.cell.color, name = "cell type")
p1

pdf(file.path(output_dir, "Umap_mCRC_epithelial_cell_type.pdf"), width = 8, height = 8)
print(p1)
dev.off()

# =============================================================================
# Normal epithelial — tissue / barplots
# =============================================================================
tumor_obj@meta.data <- tumor_obj@meta.data %>%
  mutate(
    Tissue_Type2 = case_when(
      Tissue_Type == "metastasis" ~ "metastasis",
      Tissue_Type == "primary" ~ "primary",
      orig.ident == "SP819H1-Mc1" ~ "CM819C1-normal",
      orig.ident == "SP369H1-Mc1" ~ "SP369H1-normal"
    )
  )

tumor_obj@meta.data <- tumor_obj@meta.data %>%
  mutate(
    Tissue_Type2 = factor(
      Tissue_Type2,
      levels = c("metastasis", "primary", "CM819C1-normal", "SP369H1-normal")
    )
  )

source_umap <- DimPlot(
  tumor_obj,
  group.by = "Tissue_Type2",
  reduction = "epithelial_umap.scvi",
  pt.size = 1.5,
  raster = TRUE,
  cols = c(
    "metastasis" = "grey75",
    "primary" = "palegreen3",
    "CM819C1-normal" = "red",
    "SP369H1-normal" = "blue"
  ),
  order = c("metastasis", "primary", "CM819C1-normal", "SP369H1-normal")
)

source_umap2 <- DimPlot(
  subset(tumor_obj, subset = orig.ident %in% c("SP819H1-Mc1", "SP369H1-Mc1")),
  group.by = "Tissue_Type2",
  reduction = "epithelial_umap.scvi",
  cols = c(
    "CM819C1-normal" = "red",
    "SP369H1-normal" = "blue"
  ),
  pt.size = 2,
  raster = TRUE
)

options(repr.plot.width = 10, repr.plot.height = 4)
source_umap + source_umap2

pdf(file.path(output_dir, "Umap_mCRC_sample_source.pdf"), width = 10, height = 4)
print(source_umap + source_umap2)
dev.off()

epi_df <- tumor_obj@meta.data
unique(epi_df$cell_type_all3)

plot2_df <- epi_df %>%
  dplyr::count(cell_type_all3, Tissue_Type) %>%
  group_by(cell_type_all3) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(plot2_df, aes(x = cell_type_all3, y = prop, fill = Tissue_Type)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Cell type",
    y = "Cell proportion",
    fill = "Tissue type",
    title = "Tissue-of-origin distribution within each epithelial cell type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(p2)

pdf(file.path(output_dir, "Barplot_mCRC_tissue_type_by_epithelial_subtypes.pdf"), width = 8, height = 8)
print(p2)
dev.off()

normal_epi_types <- c(
  "Neuroendocrine cells",
  "Goblet cells",
  "Transit-amplifying cells",
  "Enterocytes",
  "Stem cells",
  "Tuft cells"
)

epi_df <- epi_df %>%
  filter(cell_type_all3 %in% c(normal_epi_types, tumor_epi_types)) %>%
  mutate(epi_group = case_when(
    cell_type_all3 %in% normal_epi_types ~ "Normal_epithelial",
    cell_type_all3 %in% tumor_epi_types ~ "Tumor_epithelial"
  ))

plot3_df <- epi_df %>%
  dplyr::count(Tissue_Type, epi_group) %>%
  group_by(Tissue_Type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

plot3_df$Tissue_Type <- factor(plot3_df$Tissue_Type, levels = c("normal", "primary", "metastasis"))

p3 <- ggplot(plot3_df, aes(x = Tissue_Type, y = prop, fill = epi_group)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Tissue type",
    y = "Cell proportion",
    fill = "Group",
    title = "Normal versus tumor epithelial states across tissue types"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank())

options(repr.plot.width = 4, repr.plot.height = 4)
print(p3)

pdf(file.path(output_dir, "Barplot_mCRC_tissue_type_by_epithelial_type.pdf"), width = 8, height = 8)
print(p3)
dev.off()

# =============================================================================
# Dotplot — tumor epithelial subset (SCT)
# =============================================================================
tumor_obj2 <- subset(tumor_obj, subset = cell_type_all3 %in% tumor_epi_types)
tumor_obj2

tumor_obj2$cell_type_all3 <- factor(
  tumor_obj2$cell_type_all3,
  levels = rev(tumor_epi_types)
)

saveRDS(
  tumor_obj2,
  "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE.rds"
)

DefaultAssay(tumor_obj2) <- "SCT"

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

p4 <- DotPlot(
  tumor_obj2,
  group.by = "cell_type_all3",
  features = c(
    "NOX1", "SLC26A3", "MUC2", "MKI67", "TOP2A", "LGR5", "SMOC2",
    "APCDD1", "CHGB", "CLU", "ANXA1", "VEGFA", "SLC2A1", "TGFBI", "EMP1", "L1CAM"
  )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_fixed() +
  scale_color_gradientn(
    colors = rdylbu_colors,
    limits = c(-2.5, 2.5),
    breaks = color_breaks
  ) +
  scale_size_area(limits = c(0, 100), oob = scales::squish)

options(repr.plot.width = 10, repr.plot.height = 4)
p4

pdf(file.path(output_dir, "Dotplot_mCRC_tumor_cell_type_SCT.pdf"), width = 8, height = 3)
print(p4)
dev.off()

# =============================================================================
# Geneset AUC scores — heatmap of median AUC by subtype
# =============================================================================
auc_score_file <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_snRNAseq_epithelial_genesets_auc_score.csv"
auc_score <- read.table(auc_score_file, header = TRUE, sep = ",", row.names = 1)
head(auc_score, 3)

tumor_obj2 <- AddMetaData(tumor_obj2, auc_score)
head(tumor_obj2@meta.data, 3)

auc_cols <- c(
  "Hallmark_EMT",
  "Hallmark_TGFB",
  "Hallmark_KRAS_up",
  "HALLMARK_HYPOXIA",
  "Hallmark_angiogenesis",
  "Neuroendocrine",
  "Hallmark_WNT",
  "Stem_",
  "Cycling_TA_",
  "Epi._Secretory_All_",
  "Epi._Absorptive_All_"
)

median_df <- tumor_obj2@meta.data %>%
  select(cell_type_all3, all_of(auc_cols)) %>%
  group_by(cell_type_all3) %>%
  summarize(across(everything(), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
  tibble::column_to_rownames("cell_type_all3")

mat <- as.matrix(median_df)
mat_scaled <- scale(mat, center = TRUE, scale = TRUE)

row_order <- intersect(tumor_epi_types, rownames(mat_scaled))
mat_clamped <- pmin(pmax(mat_scaled, -1.5), 1.5)
mat_clamped <- mat_clamped[row_order, , drop = FALSE]

col_fun <- colorRamp2(
  c(-1.5, 0, 1.5),
  c("darkgreen", "white", "#be0032")
)

p5 <- Heatmap(
  mat_clamped,
  name = "scaled median AUC",
  col = col_fun,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Tumor subcluster",
  column_title = "Gene signature",
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal"
  )
)

options(repr.plot.width = 8, repr.plot.height = 4)
p5

pdf(file.path(output_dir, "Heatmap_mCRC_tumor_subtype_AUCell_score.pdf"), width = 8, height = 4)
print(p5)
dev.off()
