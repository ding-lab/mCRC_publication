#!/usr/bin/env Rscript

# InferCNV tumor-subclone composition: per-sample bars and CMETS vs non-CMETS
# stacked proportions. Subclones with <1% of cells in a library are dropped.
# Source: Revesion/inferCNV / snRNAseq_CMETS_inferCNV_subclusters.ipynb

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

output_dir <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/A_Data_preprocessing/F_inferCNV"
)
tumor_rds <- Sys.getenv(
  "MCRC_TUMOR_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE_included.rds"
)
infercnv_meta_path <- Sys.getenv(
  "MCRC_INFERCNV_CELL_METADATA",
  unset = "/diskmnt/Users2/epeng/Projects/mCRC/Revesion/infercnv_cell_metadata_20260420.csv"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

MIN_CLONE_PROP_PER_SAMPLE <- 0.01
TARGET_CELLTYPE <- "CMETS"
exclude_samples <- c("CM329C1-S1")

tumor_obj <- readRDS(tumor_rds)
infercnv_meta <- read.csv(infercnv_meta_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
rownames(infercnv_meta) <- infercnv_meta$cell
infercnv_meta_add <- infercnv_meta[, colnames(infercnv_meta) != "cell", drop = FALSE]
tumor_obj <- AddMetaData(tumor_obj, infercnv_meta_add)

md <- as.data.frame(tumor_obj@meta.data, stringsAsFactors = FALSE)
md <- md[, !duplicated(colnames(md)), drop = FALSE]
md <- if ("cell" %in% colnames(md)) as_tibble(md) else as_tibble(md, rownames = "cell")

md <- md %>%
  filter(!is.na(infercnv_tumor_clone), infercnv_tumor_clone != "") %>%
  mutate(
    orig.ident = as.character(orig.ident),
    CMETS_group = ifelse(cell_type_all == TARGET_CELLTYPE, "CMETS", "Non-CMETS")
  )

n_cells_per_sample <- md %>% count(orig.ident, name = "n_cells_in_sample")
infercnv_kept_clones <- md %>%
  count(orig.ident, infercnv_tumor_clone, name = "n_in_clone") %>%
  left_join(n_cells_per_sample, by = "orig.ident") %>%
  mutate(prop_in_sample = n_in_clone / n_cells_in_sample) %>%
  filter(prop_in_sample >= MIN_CLONE_PROP_PER_SAMPLE) %>%
  select(orig.ident, infercnv_tumor_clone)

md2 <- md %>% semi_join(infercnv_kept_clones, by = c("orig.ident", "infercnv_tumor_clone"))

clone_levels <- sort(unique(as.character(md2$infercnv_tumor_clone)))
clone_cols <- setNames(grDevices::hcl.colors(length(clone_levels), palette = "Dark 3"), clone_levels)

clone_bar_df <- md2 %>%
  count(orig.ident, infercnv_tumor_clone, name = "n_cells") %>%
  group_by(orig.ident) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()

sample_order_df <- clone_bar_df %>%
  group_by(orig.ident) %>%
  summarise(n_subclones = n_distinct(infercnv_tumor_clone), total_cells = sum(n_cells), .groups = "drop") %>%
  arrange(n_subclones, total_cells)

p_sample <- ggplot(
  clone_bar_df %>% mutate(orig.ident = factor(orig.ident, levels = sample_order_df$orig.ident)),
  aes(x = prop, y = orig.ident, fill = infercnv_tumor_clone)
) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = clone_cols, breaks = clone_levels, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Proportion of tumor clone", y = "Sample", fill = "InferCNV subclone") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

samples_with_subclones <- setdiff(
  sample_order_df$orig.ident[sample_order_df$n_subclones >= 2],
  exclude_samples
)
plot_df <- md2 %>%
  filter(orig.ident %in% samples_with_subclones) %>%
  count(orig.ident, CMETS_group, infercnv_tumor_clone, name = "n_cells") %>%
  group_by(orig.ident, CMETS_group) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    CMETS_group = factor(CMETS_group, levels = c("CMETS", "Non-CMETS")),
    infercnv_tumor_clone = factor(as.character(infercnv_tumor_clone), levels = clone_levels)
  )

p_cmets <- ggplot(plot_df, aes(x = CMETS_group, y = prop, fill = infercnv_tumor_clone)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  facet_wrap(~ orig.ident, nrow = 1) +
  scale_fill_manual(values = clone_cols, breaks = clone_levels, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Subclone proportion", fill = "InferCNV subclone") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "infercnv_tumor_clone_proportion_by_samples.pdf"), p_sample, width = 6, height = 6)
ggsave(file.path(output_dir, "infercnv_tumor_clone_proportion_by_CMETS_status.pdf"), p_cmets, width = 12, height = 5)
write.csv(clone_bar_df, file.path(output_dir, "infercnv_tumor_clone_by_sample.csv"), row.names = FALSE)
write.csv(plot_df, file.path(output_dir, "infercnv_tumor_clone_by_CMETS.csv"), row.names = FALSE)
message("Wrote InferCNV CMETS subclone bars to ", output_dir)
