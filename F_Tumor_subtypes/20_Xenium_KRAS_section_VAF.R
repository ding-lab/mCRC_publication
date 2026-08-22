#!/usr/bin/env Rscript

# Section-level KRAS VAF in CMETS vs non-CMETS tumor cells (Ext. Data 7e–f / 8f).
# One VAF per Xenium section using G12D or G12V probes, whichever is the
# sample-level call (ALT / (ALT + WT) >= 0.5).
# Source: Revesion/CMETS/CMETS_mutation_mapping.ipynb

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(ggplot2)
})

output_dir <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes"
)
mCRC_path <- Sys.getenv(
  "MCRC_XENIUM_SLIM_QS",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Spatial_driver/mCRC/1_xenium_slim_spatial_dirver_crc_5K.qs"
)
metadata_path <- Sys.getenv(
  "MCRC_XENIUM_METADATA",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mCRC <- qread(mCRC_path)
meta <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
common <- intersect(colnames(mCRC), rownames(meta))
mCRC <- subset(mCRC, cells = common)
mCRC <- AddMetaData(mCRC, meta[common, , drop = FALSE])

if ("Broad_cell_type1" %in% colnames(mCRC[[]])) {
  mCRC_tumor <- subset(mCRC, Broad_cell_type1 == "Tumor")
} else {
  mCRC_tumor <- mCRC
}

mCRC_tumor$CMETS_group <- case_when(
  mCRC_tumor$All_cell_type2 %in% c("Non-canonical", "Non_canonical") ~ "CMETS",
  mCRC_tumor$All_cell_type2 %in% c(
    "Stem-like", "Stem_like", "Proliferative-like", "Proliferative_like",
    "Intestine-like", "Intestine_like"
  ) ~ "Non-CMETS",
  TRUE ~ NA_character_
)

features <- c("KRAS-p.G12D-ALT:T", "KRAS-p.G12V-ALT:A", "KRAS-p.G12V-WT")
results_list <- list()

for (layer_name in Layers(mCRC_tumor)) {
  counts_mat <- GetAssayData(mCRC_tumor, assay = "Xenium", layer = layer_name)
  features_present <- features[features %in% rownames(counts_mat)]
  if (length(features_present) == 0) next
  counts <- counts_mat[features_present, , drop = FALSE]
  layer_cells <- colnames(counts_mat)

  meta_sub <- mCRC_tumor[[]] %>%
    rownames_to_column("cell_id") %>%
    filter(cell_id %in% layer_cells, CMETS_group %in% c("CMETS", "Non-CMETS"))
  common_cells <- intersect(layer_cells, meta_sub$cell_id)
  if (length(common_cells) == 0) next
  meta_sub <- meta_sub %>% filter(cell_id %in% common_cells)
  counts <- counts[, common_cells, drop = FALSE]

  for (grp in c("CMETS", "Non-CMETS")) {
    grp_cells <- meta_sub$cell_id[meta_sub$CMETS_group == grp]
    if (length(grp_cells) == 0) next
    grp_counts <- counts[, grp_cells, drop = FALSE]
    get_sum <- function(f) if (f %in% rownames(grp_counts)) sum(grp_counts[f, ]) else 0
    sum_g12d <- get_sum("KRAS-p.G12D-ALT:T")
    sum_g12v <- get_sum("KRAS-p.G12V-ALT:A")
    sum_wt <- get_sum("KRAS-p.G12V-WT")
    results_list[[paste(layer_name, grp, sep = "_")]] <- data.frame(
      Sample = sub("^counts\\.", "", layer_name),
      CMETS_group = grp,
      Total_Cells = length(grp_cells),
      total_KRAS_G12D_ALT = sum_g12d,
      total_KRAS_G12V_ALT = sum_g12v,
      total_KRAS_WT = sum_wt,
      VAF_KRAS_G12D = ifelse(sum_g12d + sum_wt > 0, sum_g12d / (sum_g12d + sum_wt), NA_real_),
      VAF_KRAS_G12V = ifelse(sum_g12v + sum_wt > 0, sum_g12v / (sum_g12v + sum_wt), NA_real_)
    )
  }
}

summary_table <- bind_rows(results_list)
sample_call <- summary_table %>%
  group_by(Sample) %>%
  summarise(
    sample_total_G12D_ALT = sum(total_KRAS_G12D_ALT, na.rm = TRUE),
    sample_total_G12V_ALT = sum(total_KRAS_G12V_ALT, na.rm = TRUE),
    sample_total_WT = sum(total_KRAS_WT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_VAF_G12D = ifelse(
      sample_total_G12D_ALT + sample_total_WT > 0,
      sample_total_G12D_ALT / (sample_total_G12D_ALT + sample_total_WT),
      NA_real_
    ),
    sample_VAF_G12V = ifelse(
      sample_total_G12V_ALT + sample_total_WT > 0,
      sample_total_G12V_ALT / (sample_total_G12V_ALT + sample_total_WT),
      NA_real_
    ),
    KRAS_call = case_when(
      !is.na(sample_VAF_G12D) & sample_VAF_G12D >= 0.5 &
        (is.na(sample_VAF_G12V) | sample_VAF_G12D >= sample_VAF_G12V) ~ "G12D",
      !is.na(sample_VAF_G12V) & sample_VAF_G12V >= 0.5 ~ "G12V",
      TRUE ~ NA_character_
    )
  )

plot_df <- summary_table %>%
  left_join(sample_call %>% select(Sample, KRAS_call), by = "Sample") %>%
  mutate(
    CMETS_group = factor(CMETS_group, levels = c("CMETS", "Non-CMETS")),
    VAF_plot = case_when(
      KRAS_call == "G12D" ~ VAF_KRAS_G12D,
      KRAS_call == "G12V" ~ VAF_KRAS_G12V,
      TRUE ~ NA_real_
    ),
    ALT_plot = case_when(
      KRAS_call == "G12D" ~ total_KRAS_G12D_ALT,
      KRAS_call == "G12V" ~ total_KRAS_G12V_ALT,
      TRUE ~ NA_real_
    ),
    label = paste0("(", ALT_plot, "/", total_KRAS_WT, ")"),
    facet_label = paste0(Sample, "\n", KRAS_call)
  ) %>%
  filter(!is.na(KRAS_call), !is.na(VAF_plot))

p <- ggplot(plot_df, aes(x = CMETS_group, y = VAF_plot, fill = CMETS_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(y = pmin(VAF_plot + 0.05, 1.08), label = label), size = 3) +
  facet_wrap(~ facet_label, nrow = 1) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(x = NULL, y = "KRAS VAF") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave(file.path(output_dir, "Xenium_CMETS_KRAS_section_VAF.pdf"), p, width = 14, height = 4)
write.csv(plot_df, file.path(output_dir, "Xenium_CMETS_KRAS_section_VAF.csv"), row.names = FALSE)
message("Wrote section-level KRAS VAF to ", output_dir)
