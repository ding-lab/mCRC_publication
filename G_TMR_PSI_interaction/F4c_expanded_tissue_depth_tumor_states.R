#!/usr/bin/env Rscript

# Fig. 4c / Extended Data Fig. 5p: tumor-state proportions in submucosa vs
# muscularis propria after adding 3 primary CRC samples (PE0519, HT1090, HT1085)
# to the original paired Xenium cases.
# Source: Revesion/CMETS/MP_depth/MP_depth_analysis.ipynb

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
})

output_dir <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/G_TMR_PSI_interaction"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

xenium_annotation_path <- Sys.getenv(
  "MCRC_XENIUM_METADATA",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv"
)
new_sample_path <- Sys.getenv(
  "MCRC_MP_DEPTH_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/CMETS/MP_depth"
)

xenium_annotation <- read.csv(xenium_annotation_path, header = TRUE)

PE0519C1_MPcells <- read.csv(file.path(new_sample_path, "PE0519C1-F5-S1Fp1-U2_MP_cells.csv"))
PE0519C1_tbl <- read.csv(file.path(new_sample_path, "PE0519C1_tumor_type2.csv")) %>%
  mutate(
    corrected_sample_id = "PE0519C1_Co",
    tmr_region = if_else(cell_id %in% PE0519C1_MPcells$Cell.ID, "muscle", "submucosa")
  ) %>%
  rename(All_cell_type2 = group)

HT1090C1_MPcells <- read.csv(file.path(new_sample_path, "HT1090C1-S1H1Fp1-U1_MP_cells_sushie.csv")) %>%
  filter(group == "muscle")
HT1090C1_tbl <- read.csv(file.path(new_sample_path, "HT1090C1_S1_tumor_type.csv")) %>%
  mutate(
    corrected_sample_id = "HT1090C1_Co",
    tmr_region = if_else(cell_id %in% HT1090C1_MPcells$cell_id, "muscle", "submucosa")
  ) %>%
  rename(All_cell_type2 = group)

HT1085C1H1_tbl <- read.csv(file.path(new_sample_path, "HT1085C1_S1H1_tumor_type.csv")) %>%
  mutate(corrected_sample_id = "HT1085C1_Co", tmr_region = "submucosa") %>%
  rename(All_cell_type2 = group)
HT1085C1H2_tbl <- read.csv(file.path(new_sample_path, "HT1085C1_S1H2_tumor_type.csv")) %>%
  mutate(corrected_sample_id = "HT1085C1_Co", tmr_region = "muscle") %>%
  rename(All_cell_type2 = group)

harmonize_state <- function(x) {
  case_when(
    x %in% c("Intestine_like", "Intestine-like") ~ "Intestine-like",
    x %in% c("Non_canonical", "Non-canonical") ~ "Non-canonical",
    x %in% c("Proliferative_like", "Proliferative-like") ~ "Proliferative-like",
    x %in% c("Stem_like", "Stem-like") ~ "Stem-like",
    TRUE ~ x
  )
}

additional_samples <- bind_rows(PE0519C1_tbl, HT1090C1_tbl, HT1085C1H1_tbl, HT1085C1H2_tbl) %>%
  mutate(All_cell_type2 = harmonize_state(All_cell_type2))

orig_tumor_tbl <- xenium_annotation %>%
  filter(
    Organ == "Colon",
    sample_id %in% c(
      "S13-18424-B5U1", "S13-31378-E7U1", "S15-1909-A5U1",
      "S16-38794-A4U2", "S16-44514-E7U2"
    ),
    Broad_cell_type1 == "Tumor",
    tmr_region %in% c("submucosa", "muscle")
  ) %>%
  select(corrected_sample_id, All_cell_type2, tmr_region)

tumor_subtype_prop <- bind_rows(
  orig_tumor_tbl,
  additional_samples %>% select(corrected_sample_id, All_cell_type2, tmr_region)
) %>%
  count(corrected_sample_id, tmr_region, All_cell_type2, name = "n_cells") %>%
  group_by(corrected_sample_id, tmr_region) %>%
  mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
  ungroup()

subtype_order <- c("Intestine-like", "Proliferative-like", "Stem-like", "Non-canonical")
tmr_col <- c(submucosa = "palegreen4", muscle = "firebrick1")

all_subtype_prop <- tumor_subtype_prop %>%
  filter(total_cells > 0) %>%
  mutate(
    Case_ID = as.character(corrected_sample_id),
    All_cell_type2 = factor(All_cell_type2, levels = subtype_order),
    tmr_region = factor(tmr_region, levels = c("submucosa", "muscle"))
  ) %>%
  group_by(Case_ID) %>%
  filter(n_distinct(tmr_region) == 2) %>%
  ungroup()

p_all <- ggplot(
  all_subtype_prop,
  aes(x = tmr_region, y = proportion, fill = tmr_region, color = tmr_region)
) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.55) +
  geom_line(aes(group = Case_ID), linewidth = 0.4, alpha = 0.7, color = "grey50") +
  geom_point(size = 2) +
  scale_fill_manual(values = tmr_col) +
  scale_color_manual(values = tmr_col) +
  facet_wrap(~ All_cell_type2, nrow = 1, scales = "free_y") +
  stat_compare_means(
    comparisons = list(c("submucosa", "muscle")),
    method = "wilcox.test",
    paired = TRUE,
    label = "p.format",
    tip.length = 0.01
  ) +
  labs(x = "TMR region", y = "Tumor-state proportion") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(output_dir, "F4c_expanded_tissue_depth_tumor_states.pdf"),
  p_all,
  width = 10,
  height = 4
)
write.csv(
  all_subtype_prop,
  file.path(output_dir, "F4c_expanded_tissue_depth_tumor_states.csv"),
  row.names = FALSE
)
message("Wrote expanded tissue-depth tumor-state figure to ", output_dir)
