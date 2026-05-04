# Standalone: liver metastasis cohort — Cox models and KM plots for CMETS / stem proportions.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(stringr)
  library(survival)
  library(survminer)
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

tumor_tbl_survival <- tumor_mut_tbl %>%
  mutate(
    Vital_status = case_when(
      Survival == "dead" ~ 1L,
      Survival == "alive" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

prepare_survival_data <- function(
    data,
    tissue_type_col = "tissue_type",
    tissue_type_value = "metastasis",
    time_col_days = "FU_days",
    proportions_cols = c(
      "CMETS_proportion",
      "NE_proportion",
      "APCDD1_CRC_proportion",
      "CS_proportion",
      "CIP_proportion",
      "CSP_proportion",
      "CI_proportion"
    )
) {
  filtered_data <- data %>% filter(.data[[tissue_type_col]] == tissue_type_value)

  medians <- sapply(proportions_cols, function(col) {
    median(filtered_data[[col]], na.rm = TRUE)
  })

  for (col in proportions_cols) {
    group_col_name <- paste0(col, "_tumor_group")
    med <- medians[[col]]
    message(col, " median: ", med)
    filtered_data[[group_col_name]] <- ifelse(filtered_data[[col]] <= med, "Low", "High")
  }

  time_col_months <- gsub("_days", "_months", time_col_days)
  time_col_years <- gsub("_days", "_years", time_col_days)
  filtered_data[[time_col_months]] <- filtered_data[[time_col_days]] / 30
  filtered_data[[time_col_years]] <- filtered_data[[time_col_days]] / 365.25

  filtered_data %>% filter(!is.na(.data[[time_col_days]]))
}

liver_sample_tbl <- prepare_survival_data(
  tumor_tbl_survival,
  tissue_type_col = "Organ",
  tissue_type_value = "liver",
  time_col_days = "FU_days"
)

liver_sample_tbl2 <- liver_sample_tbl %>%
  filter(Case_ID != "CM268C") %>%
  mutate(
    CMETS_low_stem = if_else(
      CMETS_proportion_tumor_group == "High" & CS_proportion_tumor_group == "Low",
      "High_CMETS_Low_stem",
      "Others"
    )
  )

table(liver_sample_tbl2$Vital_status)

fit_cmets <- coxph(Surv(FU_months, Vital_status) ~ CMETS_proportion, data = liver_sample_tbl2)
summary(fit_cmets)
exp(cbind(HR = coef(fit_cmets), confint(fit_cmets)))

fit_stem <- coxph(Surv(FU_months, Vital_status) ~ CS_proportion, data = liver_sample_tbl2)
summary(fit_stem)
exp(cbind(HR = coef(fit_stem), confint(fit_stem)))

fit_cmets_stem <- coxph(
  Surv(FU_months, Vital_status) ~ CMETS_proportion + CS_proportion,
  data = liver_sample_tbl2
)
summary(fit_cmets_stem)
exp(cbind(HR = coef(fit_cmets_stem), confint(fit_cmets_stem)))

liver_sample_tbl2$CMETS_low_stem <- factor(
  liver_sample_tbl2$CMETS_low_stem,
  levels = c("Others", "High_CMETS_Low_stem")
)

fit_cmets2 <- coxph(Surv(FU_months, Vital_status) ~ CMETS_low_stem, data = liver_sample_tbl2)
summary(fit_cmets2)
exp(cbind(HR = coef(fit_cmets2), confint(fit_cmets2)))

surv_obj <- Surv(time = liver_sample_tbl2$FU_years, event = liver_sample_tbl2$Vital_status)

fit <- survfit(surv_obj ~ CMETS_proportion_tumor_group, data = liver_sample_tbl2)
surv_diff <- survdiff(surv_obj ~ CMETS_proportion_tumor_group, data = liver_sample_tbl2)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(
  fit,
  data = liver_sample_tbl2,
  palette = c("#7570B3", "#E7298A"),
  xlab = "Years",
  ylab = "Survival Probability",
  xlim = c(0, 5),
  break.x.by = 1,
  censor = FALSE,
  risk.table = TRUE,
  pval = paste0("p value=", p_value),
  pval.coord = c(3, 0.9)
)

if (exists("set_size", mode = "function")) {
  set_size(4, 4)
}
print(surv_plot)

fit2 <- survfit(surv_obj ~ CMETS_low_stem, data = liver_sample_tbl2)
surv_diff2 <- survdiff(surv_obj ~ CMETS_low_stem, data = liver_sample_tbl2)
p_value2 <- round(1 - pchisq(surv_diff2$chisq, df = 1), 2)

surv_plot_b <- ggsurvplot(
  fit2,
  data = liver_sample_tbl2,
  palette = c("#7570B3", "#E7298A"),
  xlab = "Years",
  ylab = "Survival Probability",
  xlim = c(0, 5),
  break.x.by = 1,
  censor = FALSE,
  risk.table = TRUE,
  pval = paste0("p value=", p_value2),
  pval.coord = c(3, 0.9)
)

if (exists("set_size", mode = "function")) {
  set_size(4, 4)
}
print(surv_plot_b)

pval_data <- liver_sample_tbl2 %>%
  mutate(
    FU_months_trunc = pmin(FU_months, 36),
    Vital_status_trunc = if_else(FU_months <= 36 & Vital_status == 1, 1, 0)
  )

surv_obj_trunc <- Surv(time = pval_data$FU_months_trunc, event = pval_data$Vital_status_trunc)
surv_diff_trunc <- survdiff(surv_obj_trunc ~ CMETS_low_stem, data = pval_data)
p_value3 <- round(1 - pchisq(surv_diff_trunc$chisq, df = 1), 2)

surv_obj_full <- Surv(time = liver_sample_tbl2$FU_months, event = liver_sample_tbl2$Vital_status)
fit3 <- survfit(surv_obj_full ~ CMETS_low_stem, data = liver_sample_tbl2)

surv_plot2 <- ggsurvplot(
  fit3,
  data = liver_sample_tbl2,
  palette = c("#7570B3", "#E7298A"),
  xlab = "Months",
  ylab = "Survival Probability",
  xlim = c(0, 72),
  ylim = c(0.5, 1),
  break.x.by = 12,
  censor = FALSE,
  risk.table = TRUE,
  pval = paste0("p value (3-year) = ", p_value3),
  pval.coord = c(12, 0.9)
)

if (exists("set_size", mode = "function")) {
  set_size(4, 4)
}
print(surv_plot2)
