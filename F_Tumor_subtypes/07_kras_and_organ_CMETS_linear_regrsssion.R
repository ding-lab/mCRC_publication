# Standalone: multivariable linear model of CMETS proportion (forest plot).

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(stringr)
  library(broom)
  library(car)
  library(emmeans)
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

tumor_KRAS_tbl <- tumor_mut_tbl %>%
  mutate(
    KRAS_ROS1 = case_when(
      KRAS_mut == "Yes" & ROS1 == "mut" ~ "Both_mut",
      KRAS_mut == "Yes" ~ "KRAS_only",
      ROS1 == "mut" ~ "ROS1_only",
      KRAS_mut == "No" & ROS1 == "wt" ~ "Both_wt"
    )
  )

df <- tumor_KRAS_tbl %>%
  filter(cell_count > 50) %>%
  filter(!(Case_ID %in% c("HT413C", "CM655C", "CM329C"))) %>%
  dplyr::select(Case_ID, cell_count, CMETS_proportion, Sex, Race, Organ, KRAS_ROS1, Tx_in_6mo) %>%
  mutate(
    Organ = case_when(
      Organ %in% c("colon", "rectum") ~ "colorectum",
      TRUE ~ Organ
    ),
    KRAS_ROS1 = case_when(
      KRAS_ROS1 == "Both_wt" ~ "Both_wt",
      TRUE ~ "Either_mut"
    )
  ) %>%
  filter(complete.cases(.)) %>%
  filter(Organ %in% c("colorectum", "liver", "lung")) %>%
  mutate(
    Sex = as.factor(Sex),
    Race = as.factor(Race),
    Organ = as.factor(Organ),
    KRAS_ROS1 = as.factor(KRAS_ROS1),
    Tx_in_6mo = as.factor(Tx_in_6mo)
  )

if ("colorectum" %in% levels(df$Organ)) df$Organ <- stats::relevel(df$Organ, "colorectum")
if ("No" %in% levels(df$Tx_in_6mo)) df$Tx_in_6mo <- stats::relevel(df$Tx_in_6mo, "No")
if ("M" %in% levels(df$Sex)) df$Sex <- stats::relevel(df$Sex, "M")
if ("White" %in% levels(df$Race)) df$Race <- stats::relevel(df$Race, "White")
df$KRAS_ROS1 <- factor(df$KRAS_ROS1, levels = c("Both_wt", "Either_mut"))

fit <- lm(CMETS_proportion ~ Organ + KRAS_ROS1 + Tx_in_6mo + Sex + Race, data = df)

summary(fit)
car::Anova(fit, type = 2)
broom::tidy(fit, conf.int = TRUE)
car::vif(fit)

emmeans(fit, ~ Organ)
pairs(emmeans(fit, ~ Organ), adjust = "tukey")
emmeans(fit, ~ KRAS_ROS1)
emmeans(fit, ~ Tx_in_6mo)

par(mfrow = c(2, 2))
plot(fit)
par(mfrow = c(1, 1))

coef_df <- broom::tidy(fit, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    pretty = term |>
      gsub("^Organ", "Organ: ", x = _) |>
      gsub("^KRAS_ROS1", "KRAS/ROS1 mut: ", x = _) |>
      gsub("^Tx_in_6mo", "Tx in 6 mo: ", x = _) |>
      gsub("^Sex", "Sex: ", x = _) |>
      gsub("^Race", "Race: ", x = _)
  )

coef_df$pretty <- factor(
  coef_df$pretty,
  levels = c(
    "Sex: F",
    "Race: African American",
    "Tx in 6 mo: Yes",
    "KRAS/ROS1 mut: Either_mut",
    "Organ: liver",
    "Organ: lung"
  )
)

p_forest <- ggplot(coef_df, aes(x = estimate, y = pretty)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  labs(x = "Adjusted change in CMETS proportion", y = NULL) +
  theme_classic(base_size = 14)

if (exists("set_size", mode = "function")) {
  set_size(4, 4)
}
p_forest

pdf(file = file.path(output_dir, "mCRC_snRNA_KRAS_ROS1_Organ_forest.pdf"), width = 6, height = 5)
print(p_forest)
dev.off()

write.csv(coef_df, file.path(output_dir, "mCRC_snRNA_linear_regression_forest_plot_table.csv"), row.names = FALSE)
