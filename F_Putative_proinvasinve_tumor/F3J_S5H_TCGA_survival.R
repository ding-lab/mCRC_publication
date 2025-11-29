library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("dplyr")
library("cutpointr")
library("GSVA")
library("GSEABase")
library("biomaRt")
library("Matrix")
library("AUCell")
library("matrixStats")
library("DelayedMatrixStats")
source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')

output_dir <- getwd()

#------------------------------------------------------------
# 1. Load TCGA data
#------------------------------------------------------------
rds_path  <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/epeng/mCRC/TCGA'
tcga_crc  <- readRDS(file = file.path(rds_path, "tcga_coad.rds"))

clinical <- tcga_crc@colData

clin_df <- clinical[clinical$definition == "Primary solid Tumor",
                    c("patient",
                      "age_at_diagnosis",
                      "age_at_index",
                      "vital_status",
                      "days_to_death",
                      "days_to_last_follow_up",
                      "gender",
                      "ajcc_pathologic_stage",
                      "ajcc_pathologic_t"
                    )] |>
  as.data.frame()

#------------------------------------------------------------
# 2. Basic survival variables (OS and 3-year)
#------------------------------------------------------------
clin_df$deceased <- clin_df$vital_status == "Dead"

clin_df$overall_survival <- ifelse(
  clin_df$deceased,
  clin_df$days_to_death,
  clin_df$days_to_last_follow_up
)

# 3-year specific variables
clin_df$deceased_3yr <- dplyr::if_else(
  clin_df$deceased & clin_df$days_to_death < 1095,
  TRUE, FALSE
)

clin_df <- clin_df |>
  mutate(
    survival_3yr = case_when(
      deceased_3yr & days_to_death < 1095 ~ days_to_death,
      !deceased_3yr & days_to_last_follow_up >= 1095 ~ 1095,
      !deceased_3yr & days_to_last_follow_up < 1095 ~ days_to_last_follow_up
    ),
    stage = case_when(
      ajcc_pathologic_stage %in% c("Stage I", "Stage IA") ~ "Stage I",
      ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
      ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
      ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB") ~ "Stage IV",
      TRUE ~ ajcc_pathologic_stage
    ),
    T_stage = case_when(
      ajcc_pathologic_t %in% c("T4b", "T4", "T4a") ~ "T4",
      TRUE ~ ajcc_pathologic_t
    )
  )

#------------------------------------------------------------
# Load curated epithelial CRC gene sets
#------------------------------------------------------------
geneset_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_epithelial_genesets.rds"

geneset <- readRDS(geneset_path)

# Check structure
message("Loaded gene sets: ", length(geneset), " modules")
print(names(geneset))

# Example: preview Non-Canonical gene set
head(geneset$Non_Canonical_CRC_1)


#------------------------------------------------------------
# 3. AUCell: build rankings and AUC scores from expression
#------------------------------------------------------------
# expr: genes x samples
expr <- assay(tcga_crc, "fpkm_uq_unstrand")

sym <- rowData(tcga_crc)$gene_name
keep <- !is.na(sym) & sym != ""
expr <- expr[keep, , drop = FALSE]
sym  <- sym[keep]

# collapse duplicated symbols by sample-wise max
colMaxs2 <- function(m) {
  if (inherits(m, "DelayedMatrix")) {
    DelayedMatrixStats::colMaxs(m, na.rm = TRUE)
  } else {
    matrixStats::colMaxs(m, na.rm = TRUE)
  }
}

idx_list <- split(seq_len(nrow(expr)), sym)
expr_by_sym <- do.call(
  rbind,
  lapply(idx_list, function(ix) colMaxs2(expr[ix, , drop = FALSE]))
)
rownames(expr_by_sym) <- names(idx_list)

if (!is.matrix(expr_by_sym)) {
  expr_by_sym <- as.matrix(expr_by_sym)
}

# geneset: assumed to exist in your environment
# filter to genes present and drop tiny sets
gset_f <- lapply(geneset, function(g) intersect(g, rownames(expr_by_sym)))
gset_f <- gset_f[lengths(gset_f) >= 5]

set.seed(1)
rankings <- AUCell_buildRankings(
  expr_by_sym,
  nCores   = 1,
  plotStats = FALSE,
  verbose  = FALSE
)

auc <- AUCell_calcAUC(
  gset_f,
  rankings,
  aucMaxRank = ceiling(0.05 * nrow(expr_by_sym))
)

scores_mat <- as.matrix(getAUC(auc))  # gene sets x samples

#------------------------------------------------------------
# 4. Extract Non_Canonical_CRC_1 scores
#------------------------------------------------------------
prog_name <- "Non_Canonical_CRC_1"

if (!prog_name %in% rownames(scores_mat)) {
  stop(sprintf("Program '%s' not found in scores_mat rownames.", prog_name))
}

NC1_score <- as.numeric(scores_mat[prog_name, ])
names(NC1_score) <- colnames(scores_mat)

summary(NC1_score)

score_df <- data.frame(
  sample    = names(NC1_score),
  NC1_score = NC1_score,
  row.names = NULL
)

#------------------------------------------------------------
# 5. Merge survival + program scores
#------------------------------------------------------------
df <- clin_df |>
  mutate(sample = rownames(clin_df)) |>
  inner_join(score_df, by = "sample") |>
  mutate(
    NC1_group = ifelse(NC1_score >= median(NC1_score, na.rm = TRUE),
                       "High", "Low")
  )

keep_stages <- c("Stage I", "Stage II", "Stage III", "Stage IV")

df_sub <- df |>
  filter(stage %in% keep_stages) |>
  mutate(stage = factor(stage, levels = keep_stages))

# Prepare OS and 3-yr datasets
df_os <- df_sub |>
  filter(!is.na(overall_survival), !is.na(deceased)) |>
  mutate(overall_survival_yr = overall_survival / 365.25)

df_3yr <- df_sub |>
  filter(!is.na(survival_3yr), !is.na(deceased_3yr)) |>
  mutate(survival_3yr = survival_3yr / 365.25)

#------------------------------------------------------------
# 6. Kaplan–Meier fits
#------------------------------------------------------------
# Overall survival (time in days; x-axis labeled in days)
fit_os <- survfit(Surv(overall_survival, deceased) ~ NC1_group, data = df_os)
p_os <- ggsurvplot_facet(
  fit_os,
  data        = df_os,
  facet.by    = "stage",
  nrow        = 1, ncol = 4,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = TRUE,
  palette     = c("darkgreen", "firebrick"),
  legend.labs = c("Low", "High"),
  title       = "Overall Survival by Non-Canonical CRC-1 Program (Faceted by Stage)",
  xlab        = "Days",
  ylab        = "Survival probability"
)

# 3-year survival (time in years; censored at ~3 years)
fit_3yr <- survfit(Surv(survival_3yr, deceased_3yr) ~ NC1_group, data = df_3yr)
p_3yr <- ggsurvplot_facet(
  fit_3yr,
  data        = df_3yr,
  facet.by    = "stage",
  nrow        = 1, ncol = 4,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = TRUE,
  palette     = c("darkgreen", "firebrick"),
  legend.labs = c("Low", "High"),
  title       = "3-Year Survival by Non-Canonical CRC-1 Program (Faceted by Stage)",
  xlab        = "Years",
  ylab        = "Survival probability"
)

#------------------------------------------------------------
# 7. Save plots
#------------------------------------------------------------
# (A) 3-year survival
pdf(file.path(output_dir, "TCGA_CRC_NC1_3yr_survival_by_stage.pdf"),
    width = 14, height = 4)
print(p_3yr)
dev.off()

# (B) Overall survival (optional, comment out if not needed)
pdf(file.path(output_dir, "TCGA_CRC_NC1_overall_survival_by_stage.pdf"),
    width = 14, height = 4)
print(p_os)
dev.off()


