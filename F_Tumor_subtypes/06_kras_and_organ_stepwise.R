# Standalone: CMETS proportion by organ and KRAS mutation (Kruskal-Wallis + boxplot).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
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
  ) %>%
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

tumor_KRAS_tbl_filtered <- tumor_mut_tbl %>%
  filter(Organ %in% c("colon", "liver", "lung", "rectum")) %>%
  filter(!is.na(KRAS_mut) & cell_count > 50) %>%
  filter(!(Case_ID %in% c("HT413C", "CM655C", "CM329C"))) %>%
  mutate(
    Organ = case_when(
      Organ %in% c("colon", "rectum") ~ "colorectum",
      Organ == "liver" ~ "liver",
      Organ == "lung" ~ "lung"
    ),
    Organ_KRAS_mut = paste0(Organ, "_", KRAS_mut)
  ) %>%
  dplyr::select(Case_ID, CMETS_proportion, Organ, KRAS_mut, Organ_KRAS_mut)

tumor_KRAS_tbl_filtered$Organ_KRAS_mut <- factor(
  tumor_KRAS_tbl_filtered$Organ_KRAS_mut,
  levels = c(
    "colorectum_No", "colorectum_Yes",
    "liver_No", "liver_Yes",
    "lung_No", "lung_Yes"
  )
)

kruskal.test(CMETS_proportion ~ Organ_KRAS_mut, data = tumor_KRAS_tbl_filtered)

p_merged <- ggplot(
  tumor_KRAS_tbl_filtered,
  aes(x = Organ, y = CMETS_proportion, fill = KRAS_mut)
) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    alpha = 1,
    size = 3
  ) +
  labs(
    x = "Organ",
    y = "CMETS proportion",
    fill = "KRAS mutation"
  ) +
  theme_mydefault()

options(repr.plot.width = 8, repr.plot.height = 6)
p_merged

pdf(file = file.path(output_dir, "mCRC_snRNA_KRAS_Organ_boxplot.pdf"), width = 8, height = 5)
print(p_merged)
dev.off()
