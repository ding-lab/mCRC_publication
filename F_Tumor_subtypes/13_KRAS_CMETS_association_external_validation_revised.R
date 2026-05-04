#!/usr/bin/env Rscript

## ============================================================
## KRAS/BRAF/ROS1 mutation association with CMETS score
## External validation in GSE39582 and TCGA COAD/READ
## ============================================================
## Generated from KRAS_CMETS_association_external_validation_revised.ipynb
## Notes:
## - TP53 testing has been removed.
## - Update `outdir` and `geneset_path` before running if needed.
## ============================================================


library(GEOquery)
library(tidyverse)
library(AUCell)
library(broom)
library(matrixStats)
library(ggplot2)

## GSE39582

## ============================================================
## 0) user-defined paths
## ============================================================

outdir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#geneset_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_epithelial_genesets.rds"
geneset_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_genesets_updated.rds"

## ============================================================
## 1) helper functions
## ============================================================

extract_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_

  ## handle GPL annotations like:
  ## NM_001142807 // ACOXL // acyl-CoA oxidase-like // ...
  has_double_slash <- grepl(" // ", x, fixed = TRUE)
  x[has_double_slash] <- vapply(
    strsplit(x[has_double_slash], " // ", fixed = TRUE),
    function(y) if (length(y) >= 2) y[2] else y[1],
    character(1)
  )

  ## if multiple mappings exist, keep the first symbol
  x <- vapply(
    strsplit(x, " /// ", fixed = TRUE),
    function(y) y[1],
    character(1)
  )

  x <- trimws(x)
  x[x %in% c("", "---", "NA", "N/A")] <- NA_character_
  x
}

clean_mutation <- function(x) {
  x <- toupper(trimws(as.character(x)))
  y <- case_when(
    x %in% c("M", "MUT", "MUTANT", "MUTATED") ~ "Mutant",
    x %in% c("WT", "WILD TYPE", "WILD-TYPE", "WILDTYPE") ~ "WT",
    TRUE ~ NA_character_
  )
  factor(y, levels = c("WT", "Mutant"))
}

clean_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "UNKNOWN", "NOT AVAILABLE")] <- NA_character_

  y <- case_when(
    str_detect(x, "IV|4")  ~ "IV",
    str_detect(x, "III|3") ~ "III",
    str_detect(x, "II|2")  ~ "II",
    str_detect(x, "(^|[^V])I($|[^I])|1") ~ "I",
    TRUE ~ NA_character_
  )

  factor(y, levels = c("I", "II", "III", "IV"))
}

empty_binary_result <- function(label, n = 0) {
  empty_lm_main <- tibble(
    term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, gene = label, n = n, wilcox_p = NA_real_
  )

  empty_glm_main <- tibble(
    term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, gene = label, n = n
  )

  empty_group_summary <- tibble(
    status = character(), n = integer(), mean_CMETS = numeric(),
    median_CMETS = numeric(), sd_CMETS = numeric(), gene = character()
  )

  list(
    dat = tibble(),
    group_tab = table(factor(character())),
    stage_tab = table(factor(character()), factor(character())),
    wilcox_p = NA_real_,
    fit_lm = NULL,
    lm_tidy = tibble(),
    lm_main = empty_lm_main,
    fit_glm = NULL,
    glm_tidy = tibble(),
    glm_main = empty_glm_main,
    group_summary = empty_group_summary
  )
}

run_binary_cmets_test <- function(df, varname, label, adjust_var = "TNM_stage") {
  dat <- df %>%
    select(sample, CMETS, all_of(varname), all_of(adjust_var)) %>%
    filter(
      !is.na(CMETS),
      !is.na(.data[[varname]]),
      !is.na(.data[[adjust_var]])
    )

  dat[[varname]] <- droplevels(dat[[varname]])
  dat[[adjust_var]] <- droplevels(dat[[adjust_var]])

  if (nrow(dat) == 0) {
    warning("No non-missing samples for: ", label)
    return(empty_binary_result(label, n = 0))
  }

  if (nlevels(dat[[varname]]) != 2) {
    warning(
      "Variable ", varname,
      " does not have exactly two observed levels after filtering. Observed levels: ",
      paste(levels(dat[[varname]]), collapse = ", ")
    )
    out <- empty_binary_result(label, n = nrow(dat))
    out$dat <- dat
    out$group_tab <- table(dat[[varname]], useNA = "ifany")
    return(out)
  }

  group_tab <- table(dat[[varname]], useNA = "ifany")
  stage_tab <- table(dat[[varname]], dat[[adjust_var]])

  if (any(group_tab < 3)) {
    warning(
      "Small group size for ", label, ": ",
      paste(names(group_tab), group_tab, sep = "=", collapse = ", ")
    )
  }

  dat_wilcox <- dat %>% mutate(status = .data[[varname]])
  wilcox_p <- tryCatch(
    wilcox.test(CMETS ~ status, data = dat_wilcox)$p.value,
    error = function(e) NA_real_
  )

  lm_formula <- as.formula(paste0("CMETS ~ ", varname, " + ", adjust_var))
  fit_lm <- lm(lm_formula, data = dat)
  lm_tidy <- broom::tidy(fit_lm, conf.int = TRUE)

  target_term <- paste0(varname, levels(dat[[varname]])[2])
  lm_main <- lm_tidy %>%
    filter(term == target_term) %>%
    mutate(
      gene = label,
      n = nrow(dat),
      wilcox_p = wilcox_p
    )

  dat2 <- dat %>% mutate(CMETS_z = as.numeric(scale(CMETS)))

  glm_formula <- as.formula(paste0(varname, " ~ CMETS_z + ", adjust_var))
  fit_glm <- tryCatch(
    glm(glm_formula, data = dat2, family = binomial()),
    error = function(e) NULL
  )

  if (is.null(fit_glm)) {
    glm_tidy <- tibble()
    glm_main <- empty_binary_result(label, nrow(dat2))$glm_main
  } else {
    glm_tidy <- broom::tidy(fit_glm, conf.int = TRUE, exponentiate = TRUE)
    glm_main <- glm_tidy %>%
      filter(term == "CMETS_z") %>%
      mutate(
        gene = label,
        n = nrow(dat2)
      )
  }

  group_summary <- dat %>%
    group_by(status = .data[[varname]]) %>%
    summarise(
      n = n(),
      mean_CMETS = mean(CMETS, na.rm = TRUE),
      median_CMETS = median(CMETS, na.rm = TRUE),
      sd_CMETS = sd(CMETS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(gene = label)

  list(
    dat = dat,
    group_tab = group_tab,
    stage_tab = stage_tab,
    wilcox_p = wilcox_p,
    fit_lm = fit_lm,
    lm_tidy = lm_tidy,
    lm_main = lm_main,
    fit_glm = fit_glm,
    glm_tidy = glm_tidy,
    glm_main = glm_main,
    group_summary = group_summary
  )
}

plot_binary_box <- function(df, xvar, title_txt, file, facet_stage = FALSE) {
  dat <- df %>%
    filter(!is.na(CMETS), !is.na(.data[[xvar]]), !is.na(TNM_stage))

  p <- ggplot(dat, aes(x = .data[[xvar]], y = CMETS)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.5, size = 1.1) +
    labs(title = title_txt, x = NULL, y = "CMETS AUCell score") +
    theme_bw(base_size = 12)

  if (facet_stage) {
    p <- p + facet_wrap(~TNM_stage, nrow = 1)
  }

  ggsave(
    file.path(outdir, file),
    p,
    width = ifelse(facet_stage, 10, 4.5),
    height = ifelse(facet_stage, 3.5, 4),
    useDingbats = FALSE
  )
  p
}

## ============================================================
## 2) load CMETS gene set
## ============================================================

geneset_raw <- readRDS(geneset_path)

CMETS_genes <- unique(geneset_raw$CMETS)

if (is.null(CMETS_genes) || length(CMETS_genes) == 0) {
  stop("CMETS gene set not found in the RDS object.")
}

## ============================================================
## 3) load GSE39582
## ============================================================

gse <- getGEO("GSE39582", GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno_df <- pData(eset)
feature_df <- fData(eset)

cat("Expression matrix:", dim(expr), "\n")
cat("Phenotype data:", dim(pheno_df), "\n")
cat("Feature data:", dim(feature_df), "\n")

# ============================================================
## 4) collapse probes to gene symbols
## ============================================================

if (!"Gene Symbol" %in% colnames(feature_df)) {
  stop("'Gene Symbol' column not found in feature data.")
}

sym <- extract_symbol(feature_df[["Gene Symbol"]])

keep <- !is.na(sym) & sym != ""
expr2 <- expr[keep, , drop = FALSE]
sym2 <- sym[keep]

probe_list <- split(seq_len(nrow(expr2)), sym2)

expr_by_sym <- do.call(
  rbind,
  lapply(probe_list, function(ix) {
    matrixStats::colMaxs(expr2[ix, , drop = FALSE], na.rm = TRUE)
  })
)

rownames(expr_by_sym) <- names(probe_list)
expr_by_sym <- expr_by_sym[rownames(expr_by_sym) != "", , drop = FALSE]

cat("Collapsed expression matrix:", dim(expr_by_sym), "\n")

## ============================================================
## 5) score CMETS by AUCell
## ============================================================

cmets_overlap <- intersect(CMETS_genes, rownames(expr_by_sym))
cat("Number of overlapping CMETS genes:", length(cmets_overlap), "\n")

if (length(cmets_overlap) < 5) {
  stop("Too few overlapping CMETS genes after symbol mapping.")
}

set.seed(1)

rankings <- AUCell_buildRankings(
  expr_by_sym,
  nCores = 1,
  plotStats = FALSE,
  verbose = FALSE
)

auc <- AUCell_calcAUC(
  geneSets = list(CMETS = cmets_overlap),
  rankings = rankings,
  aucMaxRank = ceiling(0.05 * nrow(expr_by_sym))
)

cmets_score <- as.numeric(getAUC(auc)["CMETS", ])

score_df <- data.frame(
  sample = colnames(expr_by_sym),
  CMETS = cmets_score,
  stringsAsFactors = FALSE
)

write.csv(score_df, file.path(outdir, "GSE39582_CMETS_scores.csv"), row.names = FALSE)

## ============================================================
## 6) build analysis dataframe
## ============================================================

required_cols <- c(
  "kras.mutation:ch1",
  "braf.mutation:ch1",
  "tnm.stage:ch1"
)

missing_cols <- setdiff(required_cols, colnames(pheno_df))
if (length(missing_cols) > 0) {
  stop("Missing columns in pheno_df: ", paste(missing_cols, collapse = ", "))
}

analysis_df <- pheno_df %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  transmute(
    sample,
    KRAS_raw = `kras.mutation:ch1`,
    BRAF_raw = `braf.mutation:ch1`,
    TNM_raw  = `tnm.stage:ch1`
  ) %>%
  inner_join(score_df, by = "sample") %>%
  mutate(
    KRAS_status = clean_mutation(KRAS_raw),
    BRAF_status = clean_mutation(BRAF_raw),
    TNM_stage   = clean_stage(TNM_raw)
  )

write.csv(analysis_df, file.path(outdir, "GSE39582_CMETS_analysis_df_raw.csv"), row.names = FALSE)

cat("\nKRAS raw:\n")
print(table(analysis_df$KRAS_raw, useNA = "ifany"))

cat("\nBRAF raw:\n")
print(table(analysis_df$BRAF_raw, useNA = "ifany"))

cat("\nTNM raw:\n")
print(table(analysis_df$TNM_raw, useNA = "ifany"))

cat("\nCleaned KRAS:\n")
print(table(analysis_df$KRAS_status, useNA = "ifany"))

cat("\nCleaned BRAF:\n")
print(table(analysis_df$BRAF_status, useNA = "ifany"))

cat("\nCleaned TNM stage:\n")
print(table(analysis_df$TNM_stage, useNA = "ifany"))

## ============================================================
## 7) main binary analyses: KRAS and BRAF
## ============================================================

res_kras <- run_binary_cmets_test(analysis_df, "KRAS_status", "KRAS")
res_braf <- run_binary_cmets_test(analysis_df, "BRAF_status", "BRAF")

individual_lm_results <- bind_rows(
  res_kras$lm_main,
  res_braf$lm_main
) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  select(gene, n, estimate, conf.low, conf.high, p.value, wilcox_p, FDR)

individual_glm_results <- bind_rows(
  res_kras$glm_main,
  res_braf$glm_main
) %>%
  mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  select(gene, n, estimate, conf.low, conf.high, p.value, FDR)

individual_group_summary <- bind_rows(
  res_kras$group_summary,
  res_braf$group_summary
)

write.csv(individual_lm_results, file.path(outdir, "GSE39582_CMETS_individual_linear_models.csv"), row.names = FALSE)
write.csv(individual_glm_results, file.path(outdir, "GSE39582_CMETS_individual_reverse_logistic_models.csv"), row.names = FALSE)
write.csv(individual_group_summary, file.path(outdir, "GSE39582_CMETS_individual_group_summary.csv"), row.names = FALSE)

write.csv(res_kras$lm_tidy, file.path(outdir, "GSE39582_CMETS_KRAS_lm_full.csv"), row.names = FALSE)
write.csv(res_braf$lm_tidy, file.path(outdir, "GSE39582_CMETS_BRAF_lm_full.csv"), row.names = FALSE)

write.csv(res_kras$glm_tidy, file.path(outdir, "GSE39582_CMETS_KRAS_glm_full.csv"), row.names = FALSE)
write.csv(res_braf$glm_tidy, file.path(outdir, "GSE39582_CMETS_BRAF_glm_full.csv"), row.names = FALSE)

## ============================================================
## 8) KRAS/BRAF binary combined variable
## ============================================================

analysis_df2 <- analysis_df %>%
  mutate(
    MAPK_KRAS_BRAF = case_when(
      KRAS_status == "Mutant" | BRAF_status == "Mutant" ~ "Altered",
      KRAS_status == "WT" & BRAF_status == "WT" ~ "WT",
      TRUE ~ NA_character_
    ),
    MAPK_KRAS_BRAF = factor(MAPK_KRAS_BRAF, levels = c("WT", "Altered"))
  )

cat("\nKRAS/BRAF binary combined:\n")
print(table(analysis_df2$MAPK_KRAS_BRAF, useNA = "ifany"))

res_mapk_binary <- run_binary_cmets_test(analysis_df2, "MAPK_KRAS_BRAF", "KRAS/BRAF_pathway")

mapk_binary_lm <- res_mapk_binary$lm_main %>%
  mutate(FDR = p.value) %>%
  select(gene, n, estimate, conf.low, conf.high, p.value, wilcox_p, FDR)

mapk_binary_glm <- res_mapk_binary$glm_main %>%
  mutate(FDR = p.value) %>%
  select(gene, n, estimate, conf.low, conf.high, p.value, FDR)

write.csv(mapk_binary_lm, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_binary_linear_model.csv"), row.names = FALSE)
write.csv(mapk_binary_glm, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_binary_reverse_logistic.csv"), row.names = FALSE)
write.csv(res_mapk_binary$group_summary, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_binary_group_summary.csv"), row.names = FALSE)
write.csv(res_mapk_binary$lm_tidy, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_binary_lm_full.csv"), row.names = FALSE)
write.csv(res_mapk_binary$glm_tidy, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_binary_glm_full.csv"), row.names = FALSE)

mapk_binary_lm
mapk_binary_glm

## ============================================================
## 9) KRAS/BRAF 3-group variable
## ============================================================

analysis_df3 <- analysis_df %>%
  mutate(
    KRAS_BRAF_3group = case_when(
      KRAS_status == "WT" & BRAF_status == "WT" ~ "WT",
      KRAS_status == "Mutant" & BRAF_status == "WT" ~ "KRAS_mut",
      KRAS_status == "WT" & BRAF_status == "Mutant" ~ "BRAF_mut",
      KRAS_status == "Mutant" & BRAF_status == "Mutant" ~ "Double_mut",
      TRUE ~ NA_character_
    ),
    KRAS_BRAF_3group = factor(
      KRAS_BRAF_3group,
      levels = c("WT", "KRAS_mut", "BRAF_mut", "Double_mut")
    )
  )

cat("\nKRAS/BRAF 3-group:\n")
print(table(analysis_df3$KRAS_BRAF_3group, useNA = "ifany"))

dat3 <- analysis_df3 %>%
  filter(!is.na(CMETS), !is.na(KRAS_BRAF_3group), !is.na(TNM_stage))

if (sum(dat3$KRAS_BRAF_3group == "Double_mut", na.rm = TRUE) < 5) {
  dat3 <- dat3 %>%
    filter(KRAS_BRAF_3group != "Double_mut") %>%
    mutate(KRAS_BRAF_3group = droplevels(KRAS_BRAF_3group))
}

fit_3group <- lm(CMETS ~ KRAS_BRAF_3group + TNM_stage, data = dat3)
fit_3group_tidy <- broom::tidy(fit_3group, conf.int = TRUE)
## Stage-adjusted omnibus test for KRAS/BRAF 3-group effect.
omnibus_3group <- drop1(fit_3group, test = "F")

group3_summary <- dat3 %>%
  group_by(KRAS_BRAF_3group) %>%
  summarise(
    n = n(),
    mean_CMETS = mean(CMETS, na.rm = TRUE),
    median_CMETS = median(CMETS, na.rm = TRUE),
    sd_CMETS = sd(CMETS, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(fit_3group_tidy, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_3group_lm.csv"), row.names = FALSE)
write.csv(as.data.frame(omnibus_3group), file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_3group_drop1_stage_adjusted.csv"), row.names = TRUE)
write.csv(group3_summary, file.path(outdir, "GSE39582_CMETS_KRAS_BRAF_3group_summary.csv"), row.names = FALSE)

## ============================================================
## 10) combined summary tables
## ============================================================

main_linear_results <- bind_rows(
  individual_lm_results %>% mutate(model = "CMETS ~ mutation + TNM_stage"),
  mapk_binary_lm %>% mutate(model = "CMETS ~ KRAS/BRAF_pathway + TNM_stage")
)

main_reverse_logistic_results <- bind_rows(
  individual_glm_results %>% mutate(model = "mutation ~ scaled_CMETS + TNM_stage"),
  mapk_binary_glm %>% mutate(model = "KRAS/BRAF_pathway ~ scaled_CMETS + TNM_stage")
)

write.csv(main_linear_results, file.path(outdir, "GSE39582_CMETS_main_linear_results_summary.csv"), row.names = FALSE)
write.csv(main_reverse_logistic_results, file.path(outdir, "GSE39582_CMETS_main_reverse_logistic_results_summary.csv"), row.names = FALSE)

## ============================================================
## 11) plots
## ============================================================

p_kras <- plot_binary_box(analysis_df, "KRAS_status", "GSE39582: CMETS by KRAS status", "GSE39582_CMETS_by_KRAS.pdf")
p_braf <- plot_binary_box(analysis_df, "BRAF_status", "GSE39582: CMETS by BRAF status", "GSE39582_CMETS_by_BRAF.pdf")
p_mapk <- plot_binary_box(analysis_df2, "MAPK_KRAS_BRAF", "GSE39582: CMETS by KRAS/BRAF pathway status", "GSE39582_CMETS_by_KRAS_BRAF_pathway.pdf")

p_kras_stage <- plot_binary_box(analysis_df, "KRAS_status", "GSE39582: CMETS by KRAS status within TNM stage", "GSE39582_CMETS_by_KRAS_facet_stage.pdf", facet_stage = TRUE)

p_3group <- ggplot(dat3, aes(x = KRAS_BRAF_3group, y = CMETS)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.1) +
  labs(title = "GSE39582: CMETS by KRAS/BRAF 3-group", x = NULL, y = "CMETS AUCell score") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "GSE39582_CMETS_by_KRAS_BRAF_3group.pdf"), p_3group, width = 5.5, height = 4, useDingbats = FALSE)

## TCGA COAD+READ

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(AUCell)
  library(broom)
  library(ggplot2)
})

geneset_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_genesets_updated.rds"
geneset_raw <- readRDS(geneset_path)
CMETS_genes <- unique(geneset_raw$CMETS)

## ============================================================
## 1) helper functions
## ============================================================

clean_stage <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x)] <- NA_character_
  x[x %in% c("", "NA", "N/A", "NOT REPORTED", "UNKNOWN", "STAGE X", "X")] <- NA_character_

  out <- case_when(
    str_detect(x, "STAGE IV|^IV$|IVA|IVB|IVC") ~ "IV",
    str_detect(x, "STAGE III|^III$|IIIA|IIIB|IIIC") ~ "III",
    str_detect(x, "STAGE II|^II$|IIA|IIB|IIC") ~ "II",
    str_detect(x, "STAGE I|^I$|IA|IB|IC") ~ "I",
    TRUE ~ NA_character_
  )

  factor(out, levels = c("I", "II", "III", "IV"))
}

first_nonmissing <- function(...) {
  x <- list(...)
  out <- x[[1]]
  for (i in seq_along(x)) {
    idx <- is.na(out) | out == ""
    out[idx] <- x[[i]][idx]
  }
  out
}

barcode_patient <- function(x) substr(as.character(x), 1, 12)
barcode_sample  <- function(x) substr(as.character(x), 1, 16)
sample_type_code <- function(x) substr(as.character(x), 14, 15)

clean_mutation_binary <- function(x) {
  x <- as.character(x)
  factor(ifelse(!is.na(x) & x == "Mutant", "Mutant", "WT"),
         levels = c("WT", "Mutant"))
}

empty_binary_result_tcga <- function(label, n = 0) {
  empty_lm_main <- tibble(
    term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, gene = label, n = n, wilcox_p = NA_real_
  )

  empty_glm_main <- tibble(
    term = NA_character_, estimate = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_,
    conf.high = NA_real_, gene = label, n = n
  )

  list(
    dat = tibble(),
    lm_main = empty_lm_main,
    glm_main = empty_glm_main,
    lm_tidy = tibble(),
    glm_tidy = tibble(),
    group_summary = tibble(),
    wilcox = NA
  )
}

make_binary_summary <- function(df, varname, label) {
  dat <- df %>%
    select(patient, CMETS, stage, project, all_of(varname)) %>%
    filter(
      !is.na(CMETS),
      !is.na(.data[[varname]]),
      !is.na(stage),
      !is.na(project)
    )

  dat[[varname]] <- droplevels(dat[[varname]])
  dat$stage <- droplevels(dat$stage)
  dat$project <- droplevels(dat$project)

  if (nrow(dat) == 0) {
    warning("No non-missing samples for: ", label)
    return(empty_binary_result_tcga(label, n = 0))
  }

  if (nlevels(dat[[varname]]) != 2) {
    warning(
      "Variable ", varname,
      " does not have exactly two observed levels after filtering. Observed levels: ",
      paste(levels(dat[[varname]]), collapse = ", ")
    )
    out <- empty_binary_result_tcga(label, n = nrow(dat))
    out$dat <- dat
    out$group_summary <- dat %>%
      group_by(status = .data[[varname]]) %>%
      summarise(
        n = n(),
        mean_CMETS = mean(CMETS, na.rm = TRUE),
        median_CMETS = median(CMETS, na.rm = TRUE),
        sd_CMETS = sd(CMETS, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(gene = label)
    return(out)
  }

  dat_wilcox <- dat %>% mutate(status = .data[[varname]])
  wlx_p <- tryCatch(
    wilcox.test(CMETS ~ status, data = dat_wilcox)$p.value,
    error = function(e) NA_real_
  )

  fit_lm <- lm(as.formula(paste0("CMETS ~ ", varname, " + stage + project")), data = dat)
  lm_tidy <- broom::tidy(fit_lm, conf.int = TRUE)

  mut_term <- paste0(varname, levels(dat[[varname]])[2])
  lm_main <- lm_tidy %>%
    filter(term == mut_term) %>%
    mutate(
      gene = label,
      n = nrow(dat),
      wilcox_p = wlx_p
    )

  dat2 <- dat %>% mutate(CMETS_z = as.numeric(scale(CMETS)))
  fit_glm <- tryCatch(
    glm(
      as.formula(paste0(varname, " ~ CMETS_z + stage + project")),
      data = dat2,
      family = binomial()
    ),
    error = function(e) NULL
  )

  if (is.null(fit_glm)) {
    glm_tidy <- tibble()
    glm_main <- empty_binary_result_tcga(label, nrow(dat2))$glm_main
  } else {
    glm_tidy <- broom::tidy(fit_glm, conf.int = TRUE, exponentiate = TRUE)
    glm_main <- glm_tidy %>%
      filter(term == "CMETS_z") %>%
      mutate(gene = label, n = nrow(dat2))
  }

  group_summary <- dat %>%
    group_by(status = .data[[varname]]) %>%
    summarise(
      n = n(),
      mean_CMETS = mean(CMETS, na.rm = TRUE),
      median_CMETS = median(CMETS, na.rm = TRUE),
      sd_CMETS = sd(CMETS, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(gene = label)

  list(
    dat = dat,
    lm_main = lm_main,
    glm_main = glm_main,
    lm_tidy = lm_tidy,
    glm_tidy = glm_tidy,
    group_summary = group_summary,
    wilcox = wlx_p
  )
}

## ============================================================
## 2) download TCGA expression data
## ============================================================
## Notes:
## - Use one workflow consistently.
## - If your GDC release does not return "STAR - Counts",
##   inspect available workflows and switch to the current one.

query_exp <- GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp)
se_exp <- GDCprepare(query_exp)

## ============================================================
## 3) extract expression matrix
## ============================================================

assay_names <- assayNames(se_exp)
print(assay_names)

## For bulk RNA-seq gene-set scoring, prefer normalized expression
## TPM is preferred if available
expr_assay_name <- if ("tpm_unstrand" %in% assay_names) {
  "tpm_unstrand"
} else if ("fpkm_uq_unstrand" %in% assay_names) {
  "fpkm_uq_unstrand"
} else if ("fpkm_unstrand" %in% assay_names) {
  "fpkm_unstrand"
} else if ("unstranded" %in% assay_names) {
  "unstranded"
} else if ("stranded_first" %in% assay_names) {
  "stranded_first"
} else if ("stranded_second" %in% assay_names) {
  "stranded_second"
} else {
  stop("Could not find a suitable assay in the prepared SummarizedExperiment.")
}

cat("Using assay:", expr_assay_name, "\n")
expr_mat <- assay(se_exp, expr_assay_name)

## log-transform normalized expression matrices
if (grepl("tpm|fpkm", expr_assay_name, ignore.case = TRUE)) {
  expr_mat <- log2(expr_mat + 1)
}

gene_annot <- rowData(se_exp) %>% as.data.frame()

## Try to find gene symbol column
symbol_col <- c("gene_name", "external_gene_name", "Symbol", "symbol")
symbol_col <- symbol_col[symbol_col %in% colnames(gene_annot)][1]

if (is.na(symbol_col)) {
  stop("Could not find gene symbol column in rowData(se_exp).")
}

gene_symbols <- gene_annot[[symbol_col]]
keep <- !is.na(gene_symbols) & gene_symbols != ""

expr_mat <- expr_mat[keep, , drop = FALSE]
gene_symbols <- gene_symbols[keep]

## collapse duplicated symbols by max
split_idx <- split(seq_len(nrow(expr_mat)), gene_symbols)

expr_by_sym <- do.call(
  rbind,
  lapply(split_idx, function(ix) {
    if (length(ix) == 1) {
      as.numeric(expr_mat[ix, ])
    } else {
      apply(expr_mat[ix, , drop = FALSE], 2, max, na.rm = TRUE)
    }
  })
)

rownames(expr_by_sym) <- names(split_idx)
colnames(expr_by_sym) <- colnames(expr_mat)

expr_by_sym <- as.matrix(expr_by_sym)
mode(expr_by_sym) <- "numeric"

## ============================================================
## 4) keep primary tumor samples only and deduplicate by patient
## ============================================================

sample_df <- data.frame(
  barcode = colnames(expr_by_sym),
  patient = barcode_patient(colnames(expr_by_sym)),
  sample_barcode = barcode_sample(colnames(expr_by_sym)),
  sample_type = sample_type_code(colnames(expr_by_sym)),
  stringsAsFactors = FALSE
) %>%
  mutate(project = NA_character_)

meta_exp <- as.data.frame(colData(se_exp))
meta_exp$barcode <- colnames(expr_by_sym)

## recover TCGA project if available in colData
project_col <- c("project_id", "project")
project_col <- project_col[project_col %in% colnames(meta_exp)][1]

if (!is.na(project_col)) {
  sample_df <- sample_df %>%
    left_join(
      meta_exp %>% select(barcode, project_from_metadata = all_of(project_col)),
      by = "barcode"
    ) %>%
    mutate(project = project_from_metadata) %>%
    select(-project_from_metadata)
}

if (all(is.na(sample_df$project))) {
  warning("Project metadata was not found in colData(se_exp); setting project to 'TCGA-CRC' for all samples.")
  sample_df$project <- "TCGA-CRC"
}

## primary tumor only: sample type 01
sample_df <- sample_df %>%
  filter(sample_type == "01")

expr_primary <- expr_by_sym[, sample_df$barcode, drop = FALSE]

## if multiple primary aliquots per patient, keep the first aliquot
keep_barcodes <- sample_df %>%
  arrange(patient, barcode) %>%
  group_by(patient) %>%
  dplyr::slice(1) %>%
  ungroup()

expr_primary <- expr_primary[, keep_barcodes$barcode, drop = FALSE]
sample_df <- keep_barcodes

cat("Primary tumor samples kept:", ncol(expr_primary), "\n")
print(table(sample_df$project, useNA = "ifany"))

## ============================================================
## 5) CMETS score with AUCell
## ============================================================

cmets_overlap <- intersect(CMETS_genes, rownames(expr_primary))
length(cmets_overlap)

if (length(cmets_overlap) < 5) {
  stop("Too few overlapping CMETS genes in TCGA expression matrix.")
}

set.seed(1)
rankings <- AUCell_buildRankings(
  expr_primary,
  nCores = 1,
  plotStats = FALSE,
  verbose = FALSE
)

auc <- AUCell_calcAUC(
  geneSets = list(CMETS = cmets_overlap),
  rankings = rankings,
  aucMaxRank = ceiling(0.05 * nrow(expr_primary))
)

cmets_scores <- as.numeric(getAUC(auc)["CMETS", ])

score_df <- data.frame(
  barcode = colnames(expr_primary),
  patient = barcode_patient(colnames(expr_primary)),
  CMETS = cmets_scores,
  stringsAsFactors = FALSE
)

write.csv(score_df, file.path(outdir, "TCGA_CRC_CMETS_scores.csv"), row.names = FALSE)

## ============================================================
## 6) clinical data
## ============================================================

clin_coad <- GDCquery_clinic("TCGA-COAD", type = "clinical")
clin_read <- GDCquery_clinic("TCGA-READ", type = "clinical")

clin <- bind_rows(clin_coad, clin_read)

## Harmonize patient id and stage
## Common columns can vary across releases, so use fallbacks
candidate_patient_cols <- c("submitter_id", "bcr_patient_barcode", "case_submitter_id")
patient_col <- candidate_patient_cols[candidate_patient_cols %in% colnames(clin)][1]
if (is.na(patient_col)) stop("No recognizable patient id column in clinical table.")

candidate_stage_cols <- c("ajcc_pathologic_stage", "tumor_stage", "pathologic_stage", "ajcc_stage")
stage_col_present <- candidate_stage_cols[candidate_stage_cols %in% colnames(clin)]

clin2 <- clin %>%
  mutate(
    patient = .data[[patient_col]],
    stage_raw = if (length(stage_col_present) == 0) NA_character_ else
      first_nonmissing(!!!rlang::syms(stage_col_present)),
    stage = clean_stage(stage_raw)
  ) %>%
  select(patient, stage_raw, stage) %>%
  distinct()

## ============================================================
## 7) mutation data
## ============================================================

## patients already present in the expression cohort
all_patients <- score_df %>% dplyr::distinct(patient)

## Query TCGA COAD + READ masked somatic mutation data
query_maf <- TCGAbiolinks::GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

TCGAbiolinks::GDCdownload(query_maf)
maf_all <- TCGAbiolinks::GDCprepare(query_maf)

## GDCprepare for mutation data should yield a data.frame-like object
maf_all <- as.data.frame(maf_all)

## Optional mutation-barcode sanity checks
barcode_candidates <- c("Tumor_Sample_Barcode", "Tumor_Sample_UUID")
barcode_candidates <- barcode_candidates[barcode_candidates %in% colnames(maf_all)]

cat("Mutation barcode columns found:\n")
print(barcode_candidates)

if ("Tumor_Sample_Barcode" %in% barcode_candidates) {
  cat("\nExample tumor sample barcodes:\n")
  print(head(maf_all[["Tumor_Sample_Barcode"]], 10))
}

cat("\nExpression cohort patients:", nrow(all_patients), "\n")

## Keep only non-silent coding events
nonsilent_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
  "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
  "Missense_Mutation"
)

if (!"Tumor_Sample_Barcode" %in% colnames(maf_all)) {
  stop("Could not find Tumor_Sample_Barcode in the MAF table.")
}

genes_to_test <- c("KRAS", "BRAF", "ROS1")

maf_use <- maf_all %>%
  dplyr::filter(Variant_Classification %in% nonsilent_classes) %>%
  dplyr::mutate(
    patient = barcode_patient(.data[["Tumor_Sample_Barcode"]]),
    Hugo_Symbol = as.character(Hugo_Symbol)
  ) %>%
  ## keep only patients that are present in the expression cohort
  dplyr::semi_join(all_patients, by = "patient")

mut_wide <- maf_use %>%
  dplyr::filter(Hugo_Symbol %in% genes_to_test) %>%
  dplyr::distinct(patient, Hugo_Symbol) %>%
  dplyr::mutate(mut = "Mutant") %>%
  tidyr::pivot_wider(
    names_from = Hugo_Symbol,
    values_from = mut,
    values_fill = "WT"
  )

## make sure all genes_to_test exist as columns, even if no mutations are detected
for (g in genes_to_test) {
  if (!g %in% colnames(mut_wide)) {
    mut_wide[[g]] <- "WT"
  }
}

mut_df <- mut_wide %>%
  dplyr::mutate(
    KRAS = clean_mutation_binary(KRAS),
    BRAF = clean_mutation_binary(BRAF),
    ROS1 = clean_mutation_binary(ROS1),
    KRAS_BRAF_pathway = factor(
      ifelse(KRAS == "Mutant" | BRAF == "Mutant", "Mutant", "WT"),
      levels = c("WT", "Mutant")
    ),
    KRAS_BRAF_ROS1_combo = factor(
      ifelse(KRAS == "Mutant" | BRAF == "Mutant" | ROS1 == "Mutant", "Mutant", "WT"),
      levels = c("WT", "Mutant")
    ),
    KRAS_BRAF_3group = dplyr::case_when(
      KRAS == "WT" & BRAF == "WT" ~ "WT",
      KRAS == "Mutant" & BRAF == "WT" ~ "KRAS_mut",
      KRAS == "WT" & BRAF == "Mutant" ~ "BRAF_mut",
      KRAS == "Mutant" & BRAF == "Mutant" ~ "Double_mut",
      TRUE ~ NA_character_
    ),
    KRAS_BRAF_3group = factor(
      KRAS_BRAF_3group,
      levels = c("WT", "KRAS_mut", "BRAF_mut", "Double_mut")
    )
  )

## retain WT-only patients
mut_df <- all_patients %>%
  dplyr::left_join(mut_df, by = "patient") %>%
  dplyr::mutate(
    KRAS = factor(ifelse(is.na(KRAS), "WT", as.character(KRAS)), levels = c("WT", "Mutant")),
    BRAF = factor(ifelse(is.na(BRAF), "WT", as.character(BRAF)), levels = c("WT", "Mutant")),
    ROS1 = factor(ifelse(is.na(ROS1), "WT", as.character(ROS1)), levels = c("WT", "Mutant")),
    KRAS_BRAF_pathway = factor(
      ifelse(KRAS == "Mutant" | BRAF == "Mutant", "Mutant", "WT"),
      levels = c("WT", "Mutant")
    ),
    KRAS_BRAF_ROS1_combo = factor(
      ifelse(KRAS == "Mutant" | BRAF == "Mutant" | ROS1 == "Mutant", "Mutant", "WT"),
      levels = c("WT", "Mutant")
    ),
    KRAS_BRAF_3group = dplyr::case_when(
      KRAS == "WT" & BRAF == "WT" ~ "WT",
      KRAS == "Mutant" & BRAF == "WT" ~ "KRAS_mut",
      KRAS == "WT" & BRAF == "Mutant" ~ "BRAF_mut",
      KRAS == "Mutant" & BRAF == "Mutant" ~ "Double_mut",
      TRUE ~ NA_character_
    ),
    KRAS_BRAF_3group = factor(
      KRAS_BRAF_3group,
      levels = c("WT", "KRAS_mut", "BRAF_mut", "Double_mut")
    )
  )

## sanity check
print(table(mut_df$KRAS, useNA = "ifany"))
print(table(mut_df$BRAF, useNA = "ifany"))
print(table(mut_df$ROS1, useNA = "ifany"))
print(table(mut_df$KRAS_BRAF_pathway, useNA = "ifany"))
print(table(mut_df$KRAS_BRAF_ROS1_combo, useNA = "ifany"))

## ============================================================
## 8) build final analysis dataframe
## ============================================================

analysis_df <- score_df %>%
  dplyr::left_join(
    sample_df %>% dplyr::select(barcode, patient, project),
    by = c("barcode", "patient")
  ) %>%
  dplyr::left_join(clin2, by = "patient") %>%
  dplyr::left_join(mut_df, by = "patient") %>%
  dplyr::distinct(patient, .keep_all = TRUE) %>%
  dplyr::mutate(
    project = factor(project),
    stage = factor(stage, levels = c("I", "II", "III", "IV"))
  )

write.csv(analysis_df, file.path(outdir, "TCGA_CRC_CMETS_analysis_df.csv"), row.names = FALSE)

## optional sanity checks
print(table(analysis_df$KRAS, useNA = "ifany"))
print(table(analysis_df$BRAF, useNA = "ifany"))
print(table(analysis_df$ROS1, useNA = "ifany"))
print(table(analysis_df$KRAS_BRAF_pathway, useNA = "ifany"))
print(table(analysis_df$KRAS_BRAF_ROS1_combo, useNA = "ifany"))
print(table(analysis_df$stage, useNA = "ifany"))
print(table(analysis_df$project, useNA = "ifany"))

## ============================================================
## 9) main binary analyses
## ============================================================

res_kras <- make_binary_summary(analysis_df, "KRAS", "KRAS")
res_braf <- make_binary_summary(analysis_df, "BRAF", "BRAF")
res_ros1 <- make_binary_summary(analysis_df, "ROS1", "ROS1")
res_mapk <- make_binary_summary(analysis_df, "KRAS_BRAF_pathway", "KRAS/BRAF_pathway")
res_kras_braf_ros1 <- make_binary_summary(
  analysis_df,
  "KRAS_BRAF_ROS1_combo",
  "KRAS/BRAF/ROS1_combo"
)

main_linear_results <- dplyr::bind_rows(
  res_kras$lm_main,
  res_braf$lm_main,
  res_ros1$lm_main,
  res_mapk$lm_main,
  res_kras_braf_ros1$lm_main
) %>%
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  dplyr::select(gene, n, estimate, conf.low, conf.high, p.value, wilcox_p, FDR)

main_reverse_logistic_results <- dplyr::bind_rows(
  res_kras$glm_main,
  res_braf$glm_main,
  res_ros1$glm_main,
  res_mapk$glm_main,
  res_kras_braf_ros1$glm_main
) %>%
  dplyr::mutate(FDR = p.adjust(p.value, method = "BH")) %>%
  dplyr::select(gene, n, estimate, conf.low, conf.high, p.value, FDR)

write.csv(main_linear_results, file.path(outdir, "TCGA_CRC_CMETS_main_linear_results.csv"), row.names = FALSE)
write.csv(main_reverse_logistic_results, file.path(outdir, "TCGA_CRC_CMETS_main_reverse_logistic_results.csv"), row.names = FALSE)

## ============================================================
## 10) KRAS/BRAF 3-group model
## ============================================================

dat3 <- analysis_df %>%
  dplyr::filter(!is.na(CMETS), !is.na(KRAS_BRAF_3group), !is.na(stage), !is.na(project))

if (sum(dat3$KRAS_BRAF_3group == "Double_mut", na.rm = TRUE) < 5) {
  dat3 <- dat3 %>%
    dplyr::filter(KRAS_BRAF_3group != "Double_mut") %>%
    dplyr::mutate(KRAS_BRAF_3group = droplevels(KRAS_BRAF_3group))
}

fit_3group <- lm(CMETS ~ KRAS_BRAF_3group + stage + project, data = dat3)
fit_3group_tidy <- broom::tidy(fit_3group, conf.int = TRUE)
## Stage- and project-adjusted omnibus test for KRAS/BRAF 3-group effect.
omnibus_3group <- drop1(fit_3group, test = "F")

group3_summary <- dat3 %>%
  dplyr::group_by(KRAS_BRAF_3group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_CMETS = mean(CMETS, na.rm = TRUE),
    median_CMETS = median(CMETS, na.rm = TRUE),
    sd_CMETS = sd(CMETS, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(fit_3group_tidy, file.path(outdir, "TCGA_CRC_CMETS_KRAS_BRAF_3group_lm.csv"), row.names = FALSE)
write.csv(as.data.frame(omnibus_3group), file.path(outdir, "TCGA_CRC_CMETS_KRAS_BRAF_3group_drop1_stage_project_adjusted.csv"), row.names = TRUE)
write.csv(group3_summary, file.path(outdir, "TCGA_CRC_CMETS_KRAS_BRAF_3group_summary.csv"), row.names = FALSE)

## ============================================================
## 11) plots
## ============================================================

plot_box <- function(df, xvar, title_txt, file) {
  p <- ggplot(
    df %>% dplyr::filter(!is.na(CMETS), !is.na(.data[[xvar]]), !is.na(stage), !is.na(project)),
    aes(x = .data[[xvar]], y = CMETS)
  ) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.5, size = 1.1) +
    labs(title = title_txt, x = NULL, y = "CMETS AUCell score") +
    theme_bw(base_size = 12)

  ggsave(file.path(outdir, file), p, width = 4.8, height = 4, useDingbats = FALSE)
  p
}

plot_box(analysis_df, "KRAS", "TCGA CRC: CMETS by KRAS", "TCGA_CRC_CMETS_by_KRAS.pdf")
plot_box(analysis_df, "BRAF", "TCGA CRC: CMETS by BRAF", "TCGA_CRC_CMETS_by_BRAF.pdf")
plot_box(analysis_df, "ROS1", "TCGA CRC: CMETS by ROS1", "TCGA_CRC_CMETS_by_ROS1.pdf")
plot_box(analysis_df, "KRAS_BRAF_pathway", "TCGA CRC: CMETS by KRAS/BRAF pathway", "TCGA_CRC_CMETS_by_KRAS_BRAF_pathway.pdf")
plot_box(
  analysis_df,
  "KRAS_BRAF_ROS1_combo",
  "TCGA CRC: CMETS by KRAS/BRAF/ROS1 combined status",
  "TCGA_CRC_CMETS_by_KRAS_BRAF_ROS1_combo.pdf"
)

p3 <- ggplot(dat3, aes(x = KRAS_BRAF_3group, y = CMETS)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.1) +
  labs(title = "TCGA CRC: CMETS by KRAS/BRAF 3-group", x = NULL, y = "CMETS AUCell score") +
  theme_bw(base_size = 12)

ggsave(file.path(outdir, "TCGA_CRC_CMETS_by_KRAS_BRAF_3group.pdf"), p3, width = 5.5, height = 4, useDingbats = FALSE)

## Summary

## ------------------------------------
## user options
## ------------------------------------
include_tcga <- TRUE
include_ros1 <- TRUE
include_kras_braf_ros1_combo <- TRUE

## ------------------------------------
## input files from your current script
## ------------------------------------
gse_file  <- file.path(outdir, "GSE39582_CMETS_main_linear_results_summary.csv")
tcga_file <- file.path(outdir, "TCGA_CRC_CMETS_main_linear_results.csv")

## ------------------------------------
## helper functions
## ------------------------------------
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  ifelse(
    is.na(x), NA_character_,
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 4))
  )
}

fmt_ci <- function(low, high, digits = 4) {
  paste0(fmt_num(low, digits), " to ", fmt_num(high, digits))
}

clean_variable <- function(x) {
  case_when(
    x %in% c("KRAS") ~ "KRAS mutation",
    x %in% c("BRAF") ~ "BRAF mutation",
    x %in% c("ROS1") ~ "ROS1 mutation",
    x %in% c("KRAS/BRAF_pathway", "KRAS/BRAF pathway", "MAPK_KRAS_BRAF") ~ "KRAS/BRAF altered",
    x %in% c("KRAS/BRAF/ROS1_combo", "KRAS_BRAF_ROS1_combo") ~ "KRAS/BRAF/ROS1 altered",
    TRUE ~ x
  )
}

add_note <- function(beta, p) {
  case_when(
    is.na(beta) | is.na(p) ~ NA_character_,
    p < 0.05 & beta > 0 ~ "Positive association",
    p < 0.05 & beta < 0 ~ "Negative association",
    p >= 0.05 & beta > 0 ~ "Positive trend",
    p >= 0.05 & beta < 0 ~ "Negative trend",
    TRUE ~ "No clear association"
  )
}

## ------------------------------------
## read GSE39582 results
## ------------------------------------
gse_tbl <- read_csv(gse_file, show_col_types = FALSE) %>%
  mutate(
    Cohort = "GSE39582",
    Variable = clean_variable(gene),
    Model = "Linear model: CMETS ~ variable + TNM stage",
    `Adjusted beta` = fmt_num(estimate, 4),
    `95% CI` = fmt_ci(conf.low, conf.high, 4),
    `Adjusted P` = fmt_p(p.value),
    `Unadjusted P` = fmt_p(wilcox_p),
    FDR = fmt_p(FDR),
    Note = add_note(estimate, p.value)
  ) %>%
  select(
    Cohort,
    Variable,
    N = n,
    Model,
    `Adjusted beta`,
    `95% CI`,
    `Adjusted P`,
    FDR,
    `Unadjusted P`,
    Note
  )

## ------------------------------------
## read TCGA results
## ------------------------------------
if (include_tcga) {
  tcga_tbl <- read_csv(tcga_file, show_col_types = FALSE) %>%
    mutate(
      Cohort = "TCGA-COAD/READ",
      Variable = clean_variable(gene),
      Model = "Linear model: CMETS ~ variable + stage + project",
      `Adjusted beta` = fmt_num(estimate, 4),
      `95% CI` = fmt_ci(conf.low, conf.high, 4),
      `Adjusted P` = fmt_p(p.value),
      `Unadjusted P` = fmt_p(wilcox_p),
      FDR = fmt_p(FDR),
      Note = add_note(estimate, p.value)
    ) %>%
    select(
      Cohort,
      Variable,
      N = n,
      Model,
      `Adjusted beta`,
      `95% CI`,
      `Adjusted P`,
      FDR,
      `Unadjusted P`,
      Note
    )

  supp_tbl <- bind_rows(gse_tbl, tcga_tbl)
} else {
  supp_tbl <- gse_tbl
}

## ------------------------------------
## choose rows to keep
## ------------------------------------
keep_vars <- c(
  "KRAS mutation",
  "BRAF mutation",
  "KRAS/BRAF altered"
)

if (include_ros1) {
  keep_vars <- c(keep_vars, "ROS1 mutation")
}

if (include_kras_braf_ros1_combo) {
  keep_vars <- c(keep_vars, "KRAS/BRAF/ROS1 altered")
}

supp_tbl <- supp_tbl %>%
  filter(Variable %in% keep_vars) %>%
  mutate(
    Variable = factor(Variable, levels = keep_vars),
    Cohort = factor(Cohort, levels = c("GSE39582", "TCGA-COAD/READ"))
  ) %>%
  arrange(Cohort, Variable) %>%
  mutate(
    Variable = as.character(Variable),
    Cohort = as.character(Cohort)
  )

## ------------------------------------
## write output
## ------------------------------------
outfile <- if (include_tcga) {
  file.path(outdir, "Supplementary_Table_CMETS_external_bulk_validation.csv")
} else {
  file.path(outdir, "Supplementary_Table_CMETS_GSE39582_external_validation.csv")
}

write_csv(supp_tbl, outfile)

cat("\nSaved supplementary table to:\n", outfile, "\n")
print(supp_tbl)

## Forrest plot

## Build table

library(dplyr)
library(ggplot2)
library(stringr)

## order you want
plot_vars <- c(
  "KRAS mutation",
  "BRAF mutation",
  "ROS1 mutation",
  "KRAS/BRAF altered",
  "KRAS/BRAF/ROS1 altered"
)

forest_df <- supp_tbl %>%
  filter(Variable %in% plot_vars) %>%
  mutate(
    estimate = suppressWarnings(as.numeric(`Adjusted beta`)),
    ci_low_chr = str_extract(`95% CI`, "^[^ ]+"),
    ci_high_chr = str_extract(`95% CI`, "(?<=to ).*$"),
    conf.low = suppressWarnings(as.numeric(ci_low_chr)),
    conf.high = suppressWarnings(as.numeric(ci_high_chr)),
    Cohort = factor(Cohort, levels = c("GSE39582", "TCGA-COAD/READ"))
  ) %>%
  filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
  mutate(
    Variable = factor(Variable, levels = rev(plot_vars))
  )

p_forest <- ggplot(
  forest_df,
  aes(x = estimate, y = Variable)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.20, linewidth = 0.7) +
  geom_point(size = 2.8) +
  facet_wrap(~ Cohort, ncol = 1, scales = "free_y") +
  scale_x_continuous(limits = c(-0.03, NA)) +
  labs(
    x = "Adjusted difference in CMETS AUCell score (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11)
  )

print(p_forest)

## ------------------------------------
## save plot
## ------------------------------------
ggsave(
  filename = file.path(outdir, "CMETS_external_bulk_forest_plot.pdf"),
  plot = p_forest,
  width = 4,
  height = 6,
  useDingbats = FALSE
)
