# Standalone: CMETS proportion primary vs metastasis (paired / unpaired tests + boxplot).
# Run from any working directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(ggpubr)
  library(rstatix)
})

output_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Neuroendocrine_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumor_rds_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE_included.rds"
tumor_obj <- readRDS(tumor_rds_path)

sn_sample_info.df <- tumor_obj@meta.data %>%
  dplyr::select(orig.ident, Patient_ID, Tissue_Type, cell_type_all3)

sample_rename_map <- c(
  "HT413C1-Th1K2A2Nd1_2Bma1_1" = "HT413C1-Th1K2",
  "HT413C1-Th1K4A2Nd1_2Bma1_1" = "HT413C1-Th1K4",
  "CM1563C1-S1Y1" = "CM1563C1-S1",
  "CM1563C1-T1Y1" = "CM1563C1-T1",
  "CM478C2-T1Y2" = "CM478C2-T1",
  "CM478C1-T1Y2" = "CM478C1-T1"
)

merged_sample_map <- c(
  "HT413C1-Th1K2" = "HT413C1-Th1",
  "HT413C1-Th1K4" = "HT413C1-Th1"
)

metadata <- tumor_obj@meta.data %>%
  dplyr::select(orig.ident, cell_type_all3, Site_of_Origin, Tissue_Type) %>%
  mutate(
    orig.ident = dplyr::recode(orig.ident, !!!sample_rename_map),
    orig.ident = dplyr::recode(orig.ident, !!!merged_sample_map)
  )

tumor_tbl <- metadata %>%
  group_by(orig.ident) %>%
  summarise(
    Organ = dplyr::first(na.omit(Site_of_Origin)),
    tissue_type = dplyr::first(na.omit(Tissue_Type)),
    cell_count = n(),
    CMETS_count = sum(cell_type_all3 == "CMETS", na.rm = TRUE),
    CMETS_proportion = round(100 * CMETS_count / cell_count, 2),
    APCDD1_CRC_count = sum(cell_type_all3 == "APCDD1_CRC", na.rm = TRUE),
    APCDD1_CRC_proportion = round(100 * APCDD1_CRC_count / cell_count, 2),
    NE_count = sum(cell_type_all3 == "Neuroendocrine-like tumor", na.rm = TRUE),
    NE_proportion = round(100 * NE_count / cell_count, 2),
    CIP_count = sum(cell_type_all3 == "Canonical_CRC_Intestine_Proliferation", na.rm = TRUE),
    CIP_proportion = round(100 * CIP_count / cell_count, 2),
    CI_count = sum(cell_type_all3 == "Canonical_CRC_Intestine", na.rm = TRUE),
    CI_proportion = round(100 * CI_count / cell_count, 2),
    CS_count = sum(cell_type_all3 == "Canonical_CRC_Stem", na.rm = TRUE),
    CS_proportion = round(100 * CS_count / cell_count, 2),
    CSP_count = sum(cell_type_all3 == "Canonical_CRC_Stem_Proliferation", na.rm = TRUE),
    CSP_proportion = round(100 * CSP_count / cell_count, 2),
    .groups = "drop"
  )

sn_sample_info.df2 <- sn_sample_info.df %>%
  mutate(
    orig.ident = dplyr::recode(orig.ident, !!!sample_rename_map),
    orig.ident = dplyr::recode(orig.ident, !!!merged_sample_map)
  ) %>%
  distinct(orig.ident, .keep_all = TRUE)

tumor_tbl_merged <- tumor_tbl %>%
  left_join(sn_sample_info.df2, by = "orig.ident")

tissue_col <- c(primary = "#40409C", metastasis = "#BD398D")

tumor_tbl_merged$Patient_ID2 <- sub("C[0-9]+$", "", tumor_tbl_merged$Patient_ID)

plot_raw <- tumor_tbl_merged %>%
  filter(Tissue_Type %in% c("primary", "metastasis")) %>%
  mutate(Tissue_Type = factor(Tissue_Type, levels = c("primary", "metastasis"))) %>%
  filter(orig.ident != "HT413C1-Th1")

plot_pair <- plot_raw %>%
  group_by(Patient_ID2, Tissue_Type) %>%
  summarise(
    CMETS_proportion = mean(CMETS_proportion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Patient_ID2) %>%
  filter(n_distinct(Tissue_Type) == 2) %>%
  ungroup() %>%
  mutate(Tissue_Type = factor(Tissue_Type, levels = c("primary", "metastasis")))

stat_unpaired <- plot_raw %>%
  t_test(CMETS_proportion ~ Tissue_Type, paired = FALSE) %>%
  add_xy_position(x = "Tissue_Type")

stat_paired <- plot_pair %>%
  t_test(CMETS_proportion ~ Tissue_Type, paired = TRUE) %>%
  add_xy_position(x = "Tissue_Type")

ymax <- max(plot_raw$CMETS_proportion, na.rm = TRUE)
stat_unpaired$y.position <- ymax * 1.08
stat_paired$y.position <- ymax * 1.18

stat_unpaired$y.position <- 52
stat_paired$y.position <- 57

stat_unpaired$label <- paste0("unpaired t, p = ", signif(stat_unpaired$p, 3))
stat_paired$label <- paste0("paired t, p = ", signif(stat_paired$p, 3))

pbox <- ggplot(plot_raw, aes(x = Tissue_Type, y = CMETS_proportion)) +
  geom_boxplot(aes(fill = Tissue_Type), width = 0.5, outlier.shape = NA, color = "black") +
  geom_line(
    data = plot_pair,
    aes(x = Tissue_Type, y = CMETS_proportion, group = Patient_ID2),
    inherit.aes = FALSE,
    color = "grey70",
    alpha = 0.8
  ) +
  geom_point(
    data = plot_pair,
    aes(x = Tissue_Type, y = CMETS_proportion),
    inherit.aes = FALSE,
    color = "black",
    size = 2
  ) +
  scale_fill_manual(values = tissue_col) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(y = "CMETS proportion") +
  stat_pvalue_manual(stat_unpaired, label = "label", tip.length = 0.01) +
  stat_pvalue_manual(stat_paired, label = "label", tip.length = 0.01) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

options(repr.plot.width = 5, repr.plot.height = 4)
pbox

pdf(file.path(output_dir, "Barplot_mCRC_CMETS_primary_vs_metastasis.pdf"), width = 5, height = 4)
print(pbox)
dev.off()
