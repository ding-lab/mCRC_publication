# Standalone: CMETS vs APC / TP53 / KRAS mutation status (Wilcoxon + combined PDF).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(stringr)
})

support_r <- "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/support_scripts/jupyter_support_functions.R"
if (file.exists(support_r)) {
  source(support_r)
}

output_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Neuroendocrine_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumor_rds_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE_included.rds"
tumor_obj <- readRDS(tumor_rds_path)

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

tumor_tbl <- tumor_tbl %>%
  mutate(Case_ID = str_extract(orig.ident, ".*C(?=\\d)")) %>%
  filter(!is.na(Case_ID))

clinical_info_selected <- read.csv(
  "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/H_Organ_specific_adaption/Case_mutation_status_update.csv"
)

tumor_mut_tbl <- tumor_tbl %>%
  left_join(clinical_info_selected, by = "Case_ID") %>%
  filter(!(Case_ID %in% c("HT413C", "CM329C", "CM655C"))) %>%
  filter(Organ %in% c("colon", "liver", "lung", "rectum")) %>%
  mutate(
    APC_mut = case_when(APC == "wt" ~ "No", APC == "mut" ~ "Yes"),
    KRAS_mut = case_when(KRAS == "wt" ~ "No", KRAS == "mut" ~ "Yes"),
    TP53_mut = case_when(TP53 == "wt" ~ "No", TP53 == "mut" ~ "Yes")
  )

primary_mut_tbl <- tumor_mut_tbl %>% filter(tissue_type == "primary")
met_mut_tbl <- tumor_mut_tbl %>% filter(tissue_type == "metastasis")
liver_mut_tbl <- tumor_mut_tbl %>% filter(Organ == "liver")
head(tumor_mut_tbl, 3)

write.csv(tumor_mut_tbl, file.path(output_dir, "tumor_mut_tbl_NE_included.csv"), row.names = FALSE)

tumor_KRAS_tbl_filter <- tumor_mut_tbl %>% filter(!is.na(KRAS_mut) & cell_count > 50)

tumor_KRAS_tbl_filter %>%
  group_by(KRAS_mut) %>%
  summarise(
    mean_prop = mean(CMETS_proportion, na.rm = TRUE),
    median_prop = median(CMETS_proportion, na.rm = TRUE),
    n = n()
  )

wilcox.test(CMETS_proportion ~ KRAS_mut, data = tumor_KRAS_tbl_filter)

p_kras <- ggplot(tumor_KRAS_tbl_filter, aes(x = KRAS_mut, y = CMETS_proportion, fill = KRAS_mut)) +
  geom_boxplot(alpha = 0.6) +
  geom_point(alpha = 1) +
  labs(x = "KRAS mutation", y = "CMETS") +
  theme_mydefault() +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 50
  )

tumor_APC_tbl_filter <- tumor_mut_tbl %>% filter(!is.na(APC_mut) & cell_count > 50)

tumor_APC_tbl_filter %>%
  group_by(APC_mut) %>%
  summarise(
    mean_prop = mean(CMETS_proportion, na.rm = TRUE),
    median_prop = median(CMETS_proportion, na.rm = TRUE),
    n = n()
  )

wilcox.test(CMETS_proportion ~ APC_mut, data = tumor_APC_tbl_filter)

p_apc <- ggplot(tumor_APC_tbl_filter, aes(x = APC_mut, y = CMETS_proportion, fill = APC_mut)) +
  geom_boxplot(alpha = 0.6) +
  geom_point(alpha = 1) +
  labs(x = "APC mutation", y = "CMETS proportion") +
  theme_mydefault() +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 50
  )

tumor_TP53_tbl_filter <- tumor_mut_tbl %>% filter(!is.na(TP53_mut) & cell_count > 50)

tumor_TP53_tbl_filter %>%
  group_by(TP53_mut) %>%
  summarise(
    mean_prop = mean(CMETS_proportion, na.rm = TRUE),
    median_prop = median(CMETS_proportion, na.rm = TRUE),
    n = n()
  )

wilcox.test(CMETS_proportion ~ TP53_mut, data = tumor_TP53_tbl_filter)

p_tp53 <- ggplot(tumor_TP53_tbl_filter, aes(x = TP53_mut, y = CMETS_proportion, fill = TP53_mut)) +
  geom_boxplot(alpha = 0.6) +
  geom_point(alpha = 1) +
  labs(x = "TP53 mutation", y = "CMETS proportion") +
  theme_mydefault() +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = 50
  )

if (exists("set_size", mode = "function")) {
  set_size(12, 4)
}
(p_apc + p_tp53 + p_kras) + plot_layout(ncol = 3)

pdf(file = file.path(output_dir, "mCRC_snRNA_TP53_APC_KRAS_MUT_ALL_non_canonical_barplot.pdf"), width = 10, height = 5)
print((p_apc + p_tp53 + p_kras) + plot_layout(ncol = 3))
dev.off()
