library(Seurat)
library(tidyverse)
library(googlesheets4)
library(survival)
library(survminer)
library(patchwork)
library(cowplot)
library(clinfun)
library(broom)
library(car)     
library(emmeans)
library(scales)
source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')
gs4_deauth()

output_dir = getwd()

rds_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/'
tumor_obj = readRDS(file.path(rds_path,'57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean5.rds'))
tumor_obj

metadata_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata = read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL

tumor_obj <- AddMetaData(tumor_obj, metadata)

tumor_metadata = tumor_obj@meta.data %>% select('orig.ident', 'Patient_ID', 'Site_of_Origin', 
                                                'Tissue_Type', 'Primary_Side', 'Tx_in_6mo', 'cell_type_xenium')

tumor_metadata$Case_ID <-  gsub("\\d+$", "", tumor_metadata$Patient_ID)

tumor_tbl <- tumor_metadata %>%
             group_by(Case_ID, Site_of_Origin) %>% 
             summarise(cell_count = n(),
                       APCDD1_tumor_count = sum(cell_type_xenium == 'APCDD1_CRC', na.rm = TRUE),
                       APCDD1_tumor_proportion = round(100*APCDD1_tumor_count/n(),2),
                       Non_cannonical_tumor_count = sum(cell_type_xenium == 'Non_Canonical_CRC_1', na.rm = TRUE),
                       Non_cannonical_tumor_proportion = round(100*Non_cannonical_tumor_count/n(),2),
                       Proliferative_tumor_count = sum(cell_type_xenium %in% c('Canonical_CRC_Intestine_Proliferation', 'Canonical_CRC_Stem_Proliferation'), na.rm = TRUE),
                       Proliferative_tumor_proportion = round(100*Proliferative_tumor_count/n(),2),
                       Stem_like_tumor_count = sum(cell_type_xenium == 'Canonical_CRC_Stem', na.rm = TRUE),
                       Stem_like_tumor_proportion = round(100*Stem_like_tumor_count/n(),2),
                       Intestine_like_tumor_count = sum(cell_type_xenium == 'Canonical_CRC_Intestine', na.rm = TRUE),
                       Intestine_tumor_proportion = round(100*Intestine_like_tumor_count/n(),2),                       
                         tissue_type = dplyr::first(Tissue_Type),
                         Organ = dplyr::first(Site_of_Origin),
                         Case_ID = dplyr::first(Case_ID),
                         Sample_ID = dplyr::first(Patient_ID),
                         Primary_Side = dplyr::first(Primary_Side),
                         Tx_in_6mo = dplyr::first(Tx_in_6mo),
                         .groups = 'drop')

write.csv(tumor_tbl, file.path(output_dir, 'mCRC_tumor_proportion_tbl.csv'))

clinical_info_selected <- read.csv(file.path(output_dir, 'snRNAseq_clinical_info.csv'),
                                   row.names = 1)
head(clinical_info_selected, 3)


tumor_tbl_survival <- tumor_tbl %>% 
                      left_join(clinical_info_selected, by = c("Case_ID" = "Case_ID")) %>% 
                      mutate(Vital_status = case_when(Survival == 'dead' ~ 1,
                                                      Survival == 'alive' ~ 0))


prepare_survival_data <- function(data, tissue_type_col = "tissue_type", tissue_type_value = "metastasis", 
                                  time_col_days = "FU_days", 
                                  proportions_cols = list("APCDD1_tumor_proportion", 
                                                          "Non_cannonical_tumor_proportion", 
                                                          "Proliferative_tumor_proportion", 
                                                          "Stem_like_tumor_proportion", 
                                                          "Intestine_tumor_proportion")) {
  
  # Filter the data for the specified tissue type
  filtered_data <- data %>% filter(!!sym(tissue_type_col) == tissue_type_value)
  
  # Calculate the median values for each specified tumor proportion column
  medians <- lapply(proportions_cols, function(col) {
    median(filtered_data[[col]], na.rm = TRUE)
  })
  names(medians) <- proportions_cols
  
  # Create grouping variables based on the median values
  for (col in proportions_cols) {
    group_col_name <- paste0(gsub("_tumor_proportion", "", col), "_tumor_group")
    print(paste(col, 'median:', medians[[col]]))  
    filtered_data[[group_col_name]] <- ifelse(filtered_data[[col]] <= medians[[col]], "Low", "High")
  }
  
  # Convert time from days to months
  time_col_months <- gsub("_days", "_months", time_col_days)
  filtered_data[[time_col_months]] <- filtered_data[[time_col_days]] / 30
  
  # Filter out rows with NA in the time column
  result_data <- filtered_data %>% filter(!is.na(filtered_data[[time_col_days]]))
  
  return(result_data)
}


lung_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                         tissue_type_col = "Site_of_Origin", 
                                         tissue_type_value = "lung",
                                         time_col_days = "FU_days"
                                         )


liver_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                          tissue_type_col = "Organ", 
                                          tissue_type_value = "liver",
                                          time_col_days = "FU_days"
                                         )

primary_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                            tissue_type_col = "tissue_type", 
                                            tissue_type_value = "primary",
                                            time_col_days = "FU_days"
                                         )

liver_sample_tbl2 <- liver_sample_tbl %>% 
                     filter(cell_count > 100) %>%  
                     filter(Case_ID!= 'CM268C') %>% # CM268C was only followed for 1 month
                     mutate(Non_cannonical_low_stem = if_else(Non_cannonical_tumor_group == 'High' & 
                                                              Stem_like_tumor_group == 'Low', 
                                                              'High_proinvasive_Low_stem',
                                                              'Others'))


# Filter only patients with follow-up ≤ 36 months for the p-value calculation
pval_data <- liver_sample_tbl2 %>%
  mutate(FU_months_trunc = pmin(FU_months, 60))

# Create survival object with truncated time
surv_obj_trunc <- Surv(time = pval_data$FU_months_trunc, event = pval_data$Vital_status)

# Fit and get p-value from truncated data
surv_diff_trunc <- survdiff(surv_obj_trunc ~ Non_cannonical_low_stem, data = pval_data)
p_value <- round(1 - pchisq(surv_diff_trunc$chisq, df = 1), 2)

# Fit survival model with original time for plotting
surv_obj_full <- Surv(time = liver_sample_tbl2$FU_months, event = liver_sample_tbl2$Vital_status)
fit <- survfit(surv_obj_full ~ Non_cannonical_low_stem, data = liver_sample_tbl2)

# Plot full curve but show truncated p-value
surv_plot1 <- ggsurvplot(
  fit, 
  data = liver_sample_tbl2, 
  palette = c("#7570B3", "#E7298A"),
  xlab = "Months",
  ylab = "Survival Probability",
  xlim = c(0, 72),
  ylim = c(0.5, 1),  
  break.x.by = 12,    
  censor = FALSE,
  risk.table = TRUE,
  pval = paste0("p value (5-year) = ", p_value),
  pval.coord = c(12, 0.9)
)

pdf(file.path(output_dir, "snRNA_mCRC_liver_KM_plot.pdf"), width=5, height=6)
surv_plot1$plot
dev.off()

primary_sample_tbl2 <- primary_sample_tbl %>% 
                       filter(cell_count > 100) %>%   
                       mutate(Non_cannonical_low_stem = if_else(Non_cannonical_tumor_group == 'High' & 
                                                                Stem_like_tumor_group == 'Low', 
                                                                'High_proinvasive_Low_stem',
                                                                'Others'))
pval_data <- primary_sample_tbl2 %>%
  mutate(FU_months_trunc = pmin(FU_months, 60))

surv_obj_trunc <- Surv(time = pval_data$FU_months_trunc, event = pval_data$Vital_status)

surv_diff_trunc <- survdiff(surv_obj_trunc ~ Non_cannonical_low_stem, data = pval_data)
p_value <- round(1 - pchisq(surv_diff_trunc$chisq, df = 1), 2)


surv_obj <- Surv(time = primary_sample_tbl2$FU_months, event = primary_sample_tbl2$'Vital_status')

fit <- survfit(surv_obj ~ Non_cannonical_low_stem, data = primary_sample_tbl2)

# surv_diff <- survdiff(surv_obj ~ Non_cannonical_low_stem, data = colon_sample_tbl2)
# p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot2 <- ggsurvplot(fit, 
                        data = primary_sample_tbl2, #colon_sample_tbl2, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(0, 84),
                        ylim = c(0.5, 1), 
                        break.x.by = 12,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(12, 0.9)  
                       )

pdf(file.path(output_dir, "snRNA_mCRC_primary_KM_plot.pdf"), width=5, height=6)
surv_plot2$plot
dev.off()


# KRAS tumor association
tumor_KRAS_tbl <- tumor_tbl %>% 
                  left_join(clinical_info_selected, by = c("Case_ID" = "Case_ID")) %>% 
                  mutate(KRAS_alt = case_when(
                         KRAS == 'wt' ~ 'No',
                         KRAS == 'mut' ~ 'Yes',
                         KRAS == 'gain' ~ 'Yes',
                         KRAS == 'mut/gain' ~ 'Yes'
                  )) %>% 
                  mutate(KRAS_mut = case_when(
                         KRAS == 'wt' ~ 'No',
                         KRAS == 'mut' ~ 'Yes',
                         KRAS == 'gain' ~ 'No',
                         KRAS == 'mut/gain' ~ 'Yes'
                  )) %>% 
                  mutate(KRAS_ROS1 = case_when(
                         KRAS_mut == 'Yes' & ROS1 == 'mut' ~ 'Both_mut',
                         KRAS_mut == 'Yes' & ROS1 == 'wt' ~ 'KRAS_only',
                         KRAS_mut == 'No' & ROS1 == 'mut' ~ 'ROS1_only',
                         KRAS_mut == 'No' & ROS1 == 'wt' ~ 'Both_wt'
                  ))
primary_KRAS_tbl <- tumor_KRAS_tbl %>% filter(tissue_type == 'primary')
met_KRAS_tbl <- tumor_KRAS_tbl %>% filter(tissue_type == 'metastasis')
liver_KRAS_tbl <- tumor_KRAS_tbl %>% filter(Site_of_Origin == 'liver')
head(tumor_KRAS_tbl, 3)

tumor_KRAS_tbl_filtered <- tumor_KRAS_tbl %>% 
                           filter(Organ %in% c('colon', 'liver', 'lung', 'rectum')) %>% 
                           filter(!is.na(KRAS_mut) & cell_count > 100) %>% 
                           mutate(Organ = case_when(
                                          Organ %in% c('colon', 'rectum') ~ 'colorectum',
                                          Organ == 'liver' ~ 'liver',
                                          Organ == 'lung' ~ 'lung'),
                                  Organ_KRAS_mut = paste0(Organ, '_', KRAS_mut)) %>% 
                           select(Case_ID, Non_cannonical_tumor_proportion, Organ, KRAS_mut, Organ_KRAS_mut)
tumor_KRAS_tbl_filtered$Organ_KRAS_mut <- factor(tumor_KRAS_tbl_filtered$Organ_KRAS_mut, 
                                                 level = c('colorectum_No', 'colorectum_Yes',
                                                           'liver_No', 'liver_Yes',
                                                           'lung_No', 'lung_Yes'
                                                          ))

kruskal.test(Non_cannonical_tumor_proportion ~ Organ_KRAS_mut, 
             data = tumor_KRAS_tbl_filtered)  

tumor_KRAS_tbl_filtered$Organ_KRAS_mut_num <- as.numeric(tumor_KRAS_tbl_filtered$Organ_KRAS_mut)
fit <- lm(Non_cannonical_tumor_proportion ~ Organ_KRAS_mut_num, 
          data = tumor_KRAS_tbl_filtered)
summary(lm(Non_cannonical_tumor_proportion ~ Organ_KRAS_mut_num, 
           data = tumor_KRAS_tbl_filtered))


organ_col <- c(colorectum = '#C2B280', liver = 'brown', lung = 'steelblue1')

ggplot(tumor_KRAS_tbl_filtered, 
       aes(x = Organ, y = Non_cannonical_tumor_proportion, fill = KRAS_mut)) +
  geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             alpha = 1, size = 1.5) +
  labs(
    x = "Organ",
    y = "Non-canonical tumor proportion",
    fill = "KRAS mutation"
  ) +
  theme_mydefault() -> p_merged

pdf(file.path(output_dir, "snRNA_KRAS_invasion_tumor_association.pdf"), width=8, height=4)
p_merged
dev.off()

# Multivariate linear regression
# Prep: keep needed cols, drop NAs, make factors + set sensible baselines
df <- tumor_KRAS_tbl %>%
  select(cell_count, Non_cannonical_tumor_proportion, Sex, Race, Organ, KRAS_mut, Tx_in_6mo) %>%
  mutate(Organ = case_when(Organ %in% c('colon', 'rectum') ~ 'colorectum',
                           TRUE ~ Organ)) %>% 
  filter(complete.cases(.)) %>%
  filter(cell_count > 100) %>%  
  filter(Organ %in% c('colorectum', 'liver', 'lung')) %>% 
  mutate(
    Sex = as.factor(Sex),
    Race = as.factor(Race),  
    Organ    = as.factor(Organ),
    KRAS_mut = as.factor(KRAS_mut),
    Tx_in_6mo = as.factor(Tx_in_6mo)
  )

# Optional: set reference levels if present
if ("colorectum" %in% levels(df$Organ)) df$Organ <- relevel(df$Organ, "colorectum")
if ("No" %in% levels(df$KRAS_mut))      df$KRAS_mut <- relevel(df$KRAS_mut, "No")
if ("No" %in% levels(df$Tx_in_6mo))     df$Tx_in_6mo <- relevel(df$Tx_in_6mo, "No")
if ("M" %in% levels(df$Sex))            df$Sex <- relevel(df$Sex, "M")
if ("White" %in% levels(df$Race))       df$Race <- relevel(df$Race, "White")


# Fit the multivariable linear regression
fit <- lm(Non_cannonical_tumor_proportion ~ Organ + KRAS_mut + Tx_in_6mo + Sex + Race, data = df)

# Main results
summary(fit)                 # coefficients, SEs, p-values, R^2
car::Anova(fit, type = 2)    # Type II tests (robust to unbalanced designs)
broom::tidy(fit, conf.int = TRUE)

# Multicollinearity check
car::vif(fit)

# Estimated marginal means (main-effect comparisons)
emmeans(fit, ~ Organ)
pairs(emmeans(fit, ~ Organ), adjust = "tukey")

emmeans(fit, ~ KRAS_mut)
emmeans(fit, ~ Tx_in_6mo)

coef_df <- tidy(fit, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(pretty = term |>
           gsub("^Organ",    "Organ: ",    x = _) |>
           gsub("^KRAS_mut", "KRAS mut: ", x = _) |>
           gsub("^Tx_in_6mo","Tx in 6 mo: ", x = _) |>
           gsub("^Sex",      "Sex: ",      x = _) |>
           gsub("^Race",     "Race: ",     x = _))

coef_df$pretty <- factor(coef_df$pretty, level = c('Sex: F', 'Race: African American', 'Tx in 6 mo: 1',
                                                   'KRAS mut: Yes', 'Organ: liver', 'Organ: lung'))

write.csv(coef_df, file.path(output_dir, 'snRNA_KRAS_tumor_linear_multivariate_tbl.csv'))


ggplot(coef_df, aes(x = estimate, y = pretty)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
  #scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Adjusted change in non-canonical proportion", y = NULL) +
  theme_classic(base_size = 14) -> p_forest

pdf(file.path(output_dir, "snRNA_KRAS_invasion_tumor_linear_forrest.pdf"), width=8, height=4)
p_forest
dev.off()


