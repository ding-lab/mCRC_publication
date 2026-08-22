#!/usr/bin/env Rscript

# CMETS external survival validation using AUCell scores.
# Revision #2 adds MSI/MMR-adjusted Cox models (Supp Table 3f).
# Required local inputs:
#   1) WASHU snRNA CRC gene-set RDS
#   2) selected mCRC gene-set RDS containing Hallmark_EMT
#   3) Hallmark hypoxia gene-set RDS containing HALLMARK_HYPOXIA
#   4) TCGA COAD/READ SummarizedExperiment RDS files: tcga_coad.rds, tcga_read.rds
# GEO cohorts are downloaded with GEOquery at runtime.

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
  library(AUCell)
  library(survival)
  library(matrixStats)
  library(SummarizedExperiment)
  library(DelayedMatrixStats)
})

options(stringsAsFactors = FALSE)

# -----------------------------
# Paths and configuration
# -----------------------------
# For public/reproducible use, set these environment variables or edit the defaults below.
output_dir <- Sys.getenv(
  "CMETS_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes"
)

geneset_path <- Sys.getenv(
  "CMETS_GENESET_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_genesets_updated.rds"
)
selected_geneset_path <- Sys.getenv(
  "CMETS_SELECTED_GENESET_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/mCRC_selected_genesets.rds"
)
hypoxia_geneset_path <- Sys.getenv(
  "CMETS_HYPOXIA_GENESET_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/Hallmark_hpoxia_genesets.rds"
)
tcga_dir <- Sys.getenv(
  "CMETS_TCGA_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/epeng/mCRC/TCGA"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  geneset_path,
  selected_geneset_path,
  hypoxia_geneset_path,
  file.path(tcga_dir, "tcga_coad.rds"),
  file.path(tcga_dir, "tcga_read.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}
geneset <- readRDS(geneset_path)
geneset2 <- readRDS(selected_geneset_path)
geneset3 <- readRDS(hypoxia_geneset_path)


# -----------------------------
# Generic helpers
# -----------------------------
extract_symbol <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_

  # Handle GPL annotations such as:
  # NM_001142807 // ACOXL // acyl-CoA oxidase-like // ...
  has_double_slash <- grepl(" // ", x, fixed = TRUE)
  x[has_double_slash] <- vapply(
    strsplit(x[has_double_slash], " // ", fixed = TRUE),
    function(y) if (length(y) >= 2) y[2] else y[1],
    character(1)
  )

  # If multiple mappings exist, keep the first symbol.
  x <- vapply(
    strsplit(x, " /// ", fixed = TRUE),
    function(y) y[1],
    character(1)
  )

  x <- trimws(x)
  x[x %in% c("", "---", "NA", "N/A", "NULL")] <- NA_character_
  x
}

require_named_element <- function(x, nm, object_name) {
  if (!nm %in% names(x) || is.null(x[[nm]]) || length(x[[nm]]) == 0) {
    stop("Missing gene set '", nm, "' in ", object_name, ".")
  }
  unique(as.character(x[[nm]]))
}

geneset_use <- list(
  CMETS = require_named_element(geneset, "CMETS", "geneset"),
  Canonical_CRC_Stem = require_named_element(geneset, "Canonical_CRC_Stem", "geneset"),
  EMT = require_named_element(geneset2, "Hallmark_EMT", "geneset2"),
  Hypoxia = require_named_element(geneset3, "HALLMARK_HYPOXIA", "geneset3")
)

collapse_expr_by_symbol <- function(expr, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]

  idx <- split(seq_len(nrow(expr)), symbols)
  expr_by_sym <- do.call(
    rbind,
    lapply(idx, function(ix) matrixStats::colMaxs(expr[ix, , drop = FALSE], na.rm = TRUE))
  )
  rownames(expr_by_sym) <- names(idx)
  expr_by_sym[rownames(expr_by_sym) != "", , drop = FALSE]
}

colMaxs2 <- function(m) {
  if (inherits(m, "DelayedMatrix")) {
    DelayedMatrixStats::colMaxs(m, na.rm = TRUE)
  } else {
    matrixStats::colMaxs(m, na.rm = TRUE)
  }
}

score_signatures_aucell <- function(expr_by_sym, genesets, min_genes = 5) {
  # Remove genes with no variance across samples to reduce noisy rankings
  if (ncol(expr_by_sym) > 1) {
    keep_var <- matrixStats::rowSds(as.matrix(expr_by_sym), na.rm = TRUE) > 0
    expr_by_sym <- expr_by_sym[keep_var, , drop = FALSE]
  }

  gset_f <- lapply(genesets, function(g) intersect(g, rownames(expr_by_sym)))
  gset_f <- gset_f[lengths(gset_f) >= min_genes]

  if (length(gset_f) == 0) stop("No gene sets retained after overlap filtering.")

  message("Retained sets: ", paste(names(gset_f), collapse = ", "))
  print(sapply(gset_f, length))

  rankings <- AUCell::AUCell_buildRankings(
    expr_by_sym,
    plotStats = FALSE,
    verbose = FALSE
  )

  auc <- AUCell::AUCell_calcAUC(
    gset_f,
    rankings,
    aucMaxRank = ceiling(0.05 * nrow(expr_by_sym))
  )

  scores <- as.matrix(AUCell::getAUC(auc))
  as.data.frame(t(scores)) %>% rownames_to_column("sample")
}

add_signature_zscores <- function(df) {
  req <- c("CMETS", "Canonical_CRC_Stem", "EMT", "Hypoxia")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) stop("Missing signature score columns: ", paste(miss, collapse = ", "))

  df %>%
    mutate(
      CMETS_z = as.numeric(scale(CMETS)),
      Stem_z = as.numeric(scale(Canonical_CRC_Stem)),
      EMT_z = as.numeric(scale(EMT)),
      Hypoxia_z = as.numeric(scale(Hypoxia))
    )
}

extract_cox <- function(fit, model_name) {
  s <- summary(fit)
  data.frame(
    model = model_name,
    term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    CI_low = s$conf.int[, "lower .95"],
    CI_high = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"],
    n = fit$n,
    events = fit$nevent,
    row.names = NULL
  )
}

fit_standard_models <- function(df, time_col, event_col, stage_col = NULL) {
  core <- c("CMETS_z", "Stem_z", "EMT_z", "Hypoxia_z", time_col, event_col)
  dat <- df %>%
    filter(if_all(all_of(core), ~ !is.na(.x))) %>%
    filter(.data[[time_col]] > 0, .data[[event_col]] %in% c(0, 1))

  if (nrow(dat) == 0) stop("No complete survival records for ", time_col, "/", event_col, ".")
  if (sum(dat[[event_col]] == 1) == 0) stop("No events available for ", time_col, "/", event_col, ".")

  surv_formula <- function(rhs) as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs))

  fits <- list(
    "CMETS" = coxph(surv_formula("CMETS_z"), data = dat),
    "Stem" = coxph(surv_formula("Stem_z"), data = dat),
    "EMT" = coxph(surv_formula("EMT_z"), data = dat),
    "Hypoxia" = coxph(surv_formula("Hypoxia_z"), data = dat),
    "CMETS + Stem" = coxph(surv_formula("CMETS_z + Stem_z"), data = dat),
    "CMETS + EMT" = coxph(surv_formula("CMETS_z + EMT_z"), data = dat),
    "CMETS + Hypoxia" = coxph(surv_formula("CMETS_z + Hypoxia_z"), data = dat),
    "CMETS + EMT + Hypoxia" = coxph(surv_formula("CMETS_z + EMT_z + Hypoxia_z"), data = dat),
    "CMETS + Stem + EMT + Hypoxia" = coxph(surv_formula("CMETS_z + Stem_z + EMT_z + Hypoxia_z"), data = dat)
  )

  if (!is.null(stage_col)) {
    dat_stage <- dat %>% filter(!is.na(.data[[stage_col]]))
    dat_stage[[stage_col]] <- droplevels(factor(dat_stage[[stage_col]]))

    if (nrow(dat_stage) > 0 && nlevels(dat_stage[[stage_col]]) >= 2 && sum(dat_stage[[event_col]] == 1) > 0) {
      fits <- c(
        fits,
        list(
          "CMETS + stage" = coxph(surv_formula(paste0("CMETS_z + ", stage_col)), data = dat_stage),
          "CMETS + Stem + stage" = coxph(surv_formula(paste0("CMETS_z + Stem_z + ", stage_col)), data = dat_stage),
          "CMETS + EMT + stage" = coxph(surv_formula(paste0("CMETS_z + EMT_z + ", stage_col)), data = dat_stage),
          "CMETS + Hypoxia + stage" = coxph(surv_formula(paste0("CMETS_z + Hypoxia_z + ", stage_col)), data = dat_stage),
          "CMETS + EMT + Hypoxia + stage" = coxph(
            surv_formula(paste0("CMETS_z + EMT_z + Hypoxia_z + ", stage_col)),
            data = dat_stage
          ),
          "CMETS + Stem + EMT + Hypoxia + stage" = coxph(
            surv_formula(paste0("CMETS_z + Stem_z + EMT_z + Hypoxia_z + ", stage_col)),
            data = dat_stage
          )
        )
      )
    } else {
      warning("Skipping stage-adjusted models for ", time_col, "/", event_col, ": insufficient non-missing stage levels or events.")
    }
  }

  fits
}

fit_msi_models <- function(df, time_col, event_col, msi_col, stage_col = NULL) {
  core <- c("CMETS_z", "Stem_z", "EMT_z", "Hypoxia_z", time_col, event_col, msi_col)
  dat <- df %>%
    filter(if_all(all_of(core), ~ !is.na(.x))) %>%
    filter(.data[[time_col]] > 0, .data[[event_col]] %in% c(0, 1))
  dat[[msi_col]] <- droplevels(factor(dat[[msi_col]]))

  if (nrow(dat) == 0 || nlevels(dat[[msi_col]]) < 2 || sum(dat[[event_col]] == 1) == 0) {
    warning("Skipping MSI/MMR-adjusted models for ", time_col, "/", event_col, ".")
    return(list())
  }

  surv_formula <- function(rhs) as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs))
  msi_rhs <- msi_col
  fits <- list()
  fits[[paste0("CMETS + ", msi_col)]] <- coxph(surv_formula(paste0("CMETS_z + ", msi_rhs)), data = dat)
  fits[[paste0("CMETS + Stem + ", msi_col)]] <- coxph(surv_formula(paste0("CMETS_z + Stem_z + ", msi_rhs)), data = dat)
  fits[[paste0("CMETS + EMT + ", msi_col)]] <- coxph(surv_formula(paste0("CMETS_z + EMT_z + ", msi_rhs)), data = dat)
  fits[[paste0("CMETS + Hypoxia + ", msi_col)]] <- coxph(surv_formula(paste0("CMETS_z + Hypoxia_z + ", msi_rhs)), data = dat)
  fits[[paste0("CMETS + Stem + EMT + Hypoxia + ", msi_col)]] <- coxph(
    surv_formula(paste0("CMETS_z + Stem_z + EMT_z + Hypoxia_z + ", msi_rhs)),
    data = dat
  )

  if (!is.null(stage_col)) {
    dat_stage <- dat %>% filter(!is.na(.data[[stage_col]]))
    dat_stage[[stage_col]] <- droplevels(factor(dat_stage[[stage_col]]))
    if (nrow(dat_stage) > 0 && nlevels(dat_stage[[stage_col]]) >= 2 && sum(dat_stage[[event_col]] == 1) > 0) {
      fits[[paste0("CMETS + stage + ", msi_col)]] <- coxph(
        surv_formula(paste0("CMETS_z + ", stage_col, " + ", msi_rhs)),
        data = dat_stage
      )
      fits[[paste0("CMETS + Stem + EMT + Hypoxia + stage + ", msi_col)]] <- coxph(
        surv_formula(paste0("CMETS_z + Stem_z + EMT_z + Hypoxia_z + ", stage_col, " + ", msi_rhs)),
        data = dat_stage
      )
    }
  }
  fits
}

build_result_tbl <- function(fit_list, cohort, endpoint) {
  bind_rows(lapply(names(fit_list), function(nm) extract_cox(fit_list[[nm]], nm))) %>%
    as_tibble() %>%
    mutate(cohort = cohort, endpoint = endpoint) %>%
    select(cohort, endpoint, everything())
}

get_tcga_assay_name <- function(se) {
  assay_priority <- c("tpm_unstrand", "fpkm_unstrand", "fpkm_uq_unstrand", "unstranded")
  assay_hit <- assay_priority[assay_priority %in% assayNames(se)][1]
  if (is.na(assay_hit) || length(assay_hit) == 0) {
    stop("None of expected TCGA assays found. Available assays: ", paste(assayNames(se), collapse = ", "))
  }
  assay_hit
}

# -----------------------------
# Cohort-specific preparation
# -----------------------------
prep_gse39582 <- function(genesets) {
  eset <- getGEO("GSE39582", GSEMatrix = TRUE)[[1]]
  expr <- exprs(eset)
  pheno <- pData(eset)
  feature <- fData(eset)

  expr_by_sym <- collapse_expr_by_symbol(expr, extract_symbol(feature$`Gene Symbol`))
  score_df <- score_signatures_aucell(expr_by_sym, genesets)

  pheno_df <- pheno %>%
    rownames_to_column("sample") %>%
    filter(`dataset:ch1` %in% c("discovery", "validation")) %>%
    mutate(
      OS_time = suppressWarnings(as.numeric(`os.delay (months):ch1`)),
      OS_event = case_when(`os.event:ch1` == "1" ~ 1, `os.event:ch1` == "0" ~ 0, TRUE ~ NA_real_),
      stage = case_when(
        as.character(`tnm.stage:ch1`) == "1" ~ "Stage I",
        as.character(`tnm.stage:ch1`) == "2" ~ "Stage II",
        as.character(`tnm.stage:ch1`) == "3" ~ "Stage III",
        as.character(`tnm.stage:ch1`) == "4" ~ "Stage IV",
        TRUE ~ NA_character_
      ),
      stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
      mmr_raw = as.character(`mmr.status:ch1`),
      mmr_status = case_when(
        mmr_raw == "pMMR" ~ "pMMR",
        mmr_raw == "dMMR" ~ "dMMR",
        TRUE ~ NA_character_
      ),
      mmr_status = factor(mmr_status, levels = c("pMMR", "dMMR"))
    )

  pheno_df %>% inner_join(score_df, by = "sample") %>% add_signature_zscores()
}

prep_gse17536 <- function(genesets) {
  eset <- getGEO("GSE17536", GSEMatrix = TRUE)[[1]]
  expr <- exprs(eset)
  pheno <- pData(eset)
  feature <- fData(eset)

  expr_by_sym <- collapse_expr_by_symbol(expr, extract_symbol(feature$`Gene Symbol`))
  score_df <- score_signatures_aucell(expr_by_sym, genesets)

  pheno_df <- pheno %>%
    rownames_to_column("sample") %>%
    mutate(
      OS_time = suppressWarnings(as.numeric(`overall survival follow-up time:ch1`)),
      OS_event = case_when(
        `overall_event (death from any cause):ch1` == "death" ~ 1,
        `overall_event (death from any cause):ch1` == "no death" ~ 0,
        TRUE ~ NA_real_
      ),
      stage = case_when(
        as.character(`ajcc_stage:ch1`) == "1" ~ "Stage I",
        as.character(`ajcc_stage:ch1`) == "2" ~ "Stage II",
        as.character(`ajcc_stage:ch1`) == "3" ~ "Stage III",
        as.character(`ajcc_stage:ch1`) == "4" ~ "Stage IV",
        TRUE ~ NA_character_
      ),
      stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
    )

  pheno_df %>% inner_join(score_df, by = "sample") %>% add_signature_zscores()
}

prep_gse41258 <- function(genesets) {
  eset <- getGEO("GSE41258", GSEMatrix = TRUE)[[1]]
  expr <- exprs(eset)
  pheno <- pData(eset)
  feature <- fData(eset)

  expr_by_sym <- collapse_expr_by_symbol(log2(expr + 1), extract_symbol(feature$`Gene Symbol`))
  score_df <- score_signatures_aucell(expr_by_sym, genesets)

  pheno_df <- pheno %>%
    rownames_to_column("sample") %>%
    mutate(
      tissue = as.character(`tissue:ch1`),
      OS_time = suppressWarnings(as.numeric(`fup interval:ch1`)),
      fup_status = toupper(as.character(`fup status (ned,  no evidence of disease, awd, alive with disease, aun, alive unknown, dod, dead of disease, doc, dead of other cause, dun, dead, cause unknown):ch1`)),
      OS_event = case_when(
        fup_status %in% c("DOD", "DOC", "DUN") ~ 1,
        fup_status %in% c("NED", "AWD", "AUN") ~ 0,
        TRUE ~ NA_real_
      ),
      stage = case_when(
        as.character(`group stage:ch1`) %in% c("I", "1", "Stage I", "stage i") ~ "Stage I",
        as.character(`group stage:ch1`) %in% c("II", "2", "Stage II", "stage ii") ~ "Stage II",
        as.character(`group stage:ch1`) %in% c("III", "3", "Stage III", "stage iii") ~ "Stage III",
        as.character(`group stage:ch1`) %in% c("IV", "4", "Stage IV", "stage iv") ~ "Stage IV",
        TRUE ~ NA_character_
      ),
      stage = factor(stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
      msi_raw = as.character(`microsattelite instability:ch1`),
      msi_status = case_when(
        msi_raw %in% c("MSS", "MSI-low") ~ "MSS/MSI-low",
        msi_raw %in% c("MSI-high") ~ "MSI-high",
        TRUE ~ NA_character_
      ),
      msi_status = factor(msi_status, levels = c("MSS/MSI-low", "MSI-high"))
    ) %>%
    filter(tissue == "Primary Tumor")

  pheno_df %>% inner_join(score_df, by = "sample") %>% add_signature_zscores()
}

prep_gse159216 <- function(genesets) {
  eset <- getGEO("GSE159216", GSEMatrix = TRUE)[[1]]
  expr <- exprs(eset)
  pheno <- pData(eset)
  feature <- fData(eset)

  is_control <- grepl("Housekeeping|Control", feature$gene_assignment, ignore.case = TRUE)
  is_control[is.na(is_control)] <- FALSE
  symbols <- sapply(strsplit(as.character(feature$gene_assignment), " /// ", fixed = TRUE), function(blocks) {
    first_block <- blocks[1]
    parts <- strsplit(first_block, " // ", fixed = TRUE)[[1]]
    if (length(parts) >= 2) parts[2] else NA_character_
  })
  symbols <- trimws(symbols)
  symbols[symbols %in% c("", "---", "NA", "NULL")] <- NA
  symbols[is_control] <- NA

  expr_by_sym <- collapse_expr_by_symbol(expr, symbols)
  score_df <- score_signatures_aucell(expr_by_sym, genesets)

  pheno %>%
    rownames_to_column("sample") %>%
    mutate(
      css60_status = suppressWarnings(as.numeric(`60 months cancer-specific survival, status:ch1`)),
      css60_time = suppressWarnings(as.numeric(`60 months cancer-specific survival, time:ch1`)),
      css24_time = pmin(css60_time, 24),
      css24_status = ifelse(css60_status == 1 & css60_time <= 24, 1, 0),
      css36_time = pmin(css60_time, 36),
      css36_status = ifelse(css60_status == 1 & css60_time <= 36, 1, 0),
      css48_time = pmin(css60_time, 48),
      css48_status = ifelse(css60_status == 1 & css60_time <= 48, 1, 0),
      tx_raw = as.character(`chemotherapy prior to tumor sampling:ch1`),
      tx_binary = case_when(
        tx_raw == "Naive" ~ "Naive",
        tx_raw %in% c("NeoAdjuvant", "Previous_chemotherapy", "Previous_chemotherapy/NeoAdjuvant") ~ "Treated",
        TRUE ~ NA_character_
      ),
      tx_binary = factor(tx_binary, levels = c("Naive", "Treated"))
    ) %>%
    inner_join(score_df, by = "sample") %>%
    filter(!is.na(css60_time), !is.na(css60_status)) %>%
    add_signature_zscores()
}

prep_tcga_surv <- function(tcga_obj, cancer_type, genesets) {
  clinical <- as.data.frame(colData(tcga_obj))
  clin_df <- clinical[
    clinical$definition == "Primary solid Tumor",
    c("patient", "vital_status", "days_to_death", "days_to_last_follow_up", "ajcc_pathologic_stage"),
    drop = FALSE
  ] %>%
    mutate(
      days_to_death = suppressWarnings(as.numeric(days_to_death)),
      days_to_last_follow_up = suppressWarnings(as.numeric(days_to_last_follow_up)),
      deceased = as.integer(vital_status == "Dead"),
      overall_survival = ifelse(deceased == 1, days_to_death, days_to_last_follow_up),
      stage = case_when(
        ajcc_pathologic_stage %in% c("Stage I", "Stage IA") ~ "1",
        ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "2",
        ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "3",
        ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB") ~ "4",
        TRUE ~ NA_character_
      ),
      stage = factor(stage, levels = c("1", "2", "3", "4"), ordered = FALSE)
    )

  assay_use <- get_tcga_assay_name(tcga_obj)
  expr <- assay(tcga_obj, assay_use)
  sym <- rowData(tcga_obj)$gene_name
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]
  sym <- sym[keep]

  idx_list <- split(seq_len(nrow(expr)), sym)
  expr_by_sym <- do.call(rbind, lapply(idx_list, function(ix) colMaxs2(expr[ix, , drop = FALSE])))
  rownames(expr_by_sym) <- names(idx_list)
  if (!is.matrix(expr_by_sym)) expr_by_sym <- as.matrix(expr_by_sym)

  score_df <- score_signatures_aucell(expr_by_sym, genesets)

  cd <- as.data.frame(colData(tcga_obj))
  cd$sample_id_internal <- rownames(cd)

  if ("sample" %in% colnames(score_df)) {
    score_df <- dplyr::rename(score_df, sample_id_internal = sample)
  } else if ("rowname" %in% colnames(score_df)) {
    score_df <- dplyr::rename(score_df, sample_id_internal = rowname)
  } else {
    stop("score_signatures_aucell() output is missing the sample identifier column.")
  }

  score_df <- score_df[match(cd$sample_id_internal, score_df$sample_id_internal), , drop = FALSE]
  stopifnot(identical(cd$sample_id_internal, score_df$sample_id_internal))

  tcga_df <- cd %>% inner_join(score_df, by = "sample_id_internal")
  if (!"paper_MSI_status" %in% colnames(tcga_df)) {
    tcga_df$paper_MSI_status <- NA_character_
  }
  clin_df2 <- clin_df %>%
    rownames_to_column("sample_id_internal") %>%
    select(sample_id_internal, deceased, overall_survival, stage)

  tcga_df %>%
    left_join(clin_df2, by = "sample_id_internal") %>%
    mutate(
      cancer_type = cancer_type,
      paper_MSI_status = as.character(paper_MSI_status),
      msi_status = case_when(
        paper_MSI_status %in% c("MSS", "MSI-L") ~ "MSS/MSI-L",
        paper_MSI_status == "MSI-H" ~ "MSI-H",
        TRUE ~ NA_character_
      ),
      msi_status = factor(msi_status, levels = c("MSS/MSI-L", "MSI-H"))
    ) %>%
    filter(
      !is.na(overall_survival),
      !is.na(deceased),
      !is.na(CMETS),
      !is.na(Canonical_CRC_Stem),
      !is.na(EMT),
      !is.na(Hypoxia)
    ) %>%
    mutate(sample = sample_id_internal) %>%
    add_signature_zscores()
}

# -----------------------------
# Run analyses
# -----------------------------
gse39582_df <- prep_gse39582(geneset_use)
gse17536_df <- prep_gse17536(geneset_use)
gse41258_df <- prep_gse41258(geneset_use)
gse159216_df <- prep_gse159216(geneset_use)
tcga_coad <- readRDS(file.path(tcga_dir, "tcga_coad.rds"))
tcga_read <- readRDS(file.path(tcga_dir, "tcga_read.rds"))
tcga_coad_df <- prep_tcga_surv(tcga_coad, "COAD", geneset_use)
tcga_read_df <- prep_tcga_surv(tcga_read, "READ", geneset_use)
tcga_df <- bind_rows(tcga_coad_df, tcga_read_df) %>%
  mutate(cancer_type = factor(cancer_type, levels = c("COAD", "READ")))

fits_39582 <- fit_standard_models(gse39582_df, time_col = "OS_time", event_col = "OS_event", stage_col = "stage")
fits_17536 <- fit_standard_models(gse17536_df, time_col = "OS_time", event_col = "OS_event", stage_col = "stage")
fits_41258 <- fit_standard_models(gse41258_df, time_col = "OS_time", event_col = "OS_event", stage_col = "stage")
fits_tcga <- fit_standard_models(tcga_df, time_col = "overall_survival", event_col = "deceased", stage_col = "stage")

# Revision #2 Supp Table 3f: CMETS Cox models adjusted for MSI / MMR.
fits_39582_mmr <- fit_msi_models(gse39582_df, "OS_time", "OS_event", "mmr_status", "stage")
fits_41258_msi <- fit_msi_models(gse41258_df, "OS_time", "OS_event", "msi_status", "stage")
fits_tcga_msi <- fit_msi_models(tcga_df, "overall_survival", "deceased", "msi_status", "stage")

css_dat <- gse159216_df %>%
  filter(!is.na(CMETS_z), !is.na(Stem_z), !is.na(EMT_z), !is.na(Hypoxia_z))
css_dat_tx <- css_dat %>% filter(!is.na(tx_binary))

fits_159216 <- list(
  "CSS60: CMETS" = coxph(Surv(css60_time, css60_status) ~ CMETS_z, data = css_dat),
  "CSS24: CMETS" = coxph(Surv(css24_time, css24_status) ~ CMETS_z, data = css_dat),
  "CSS36: CMETS" = coxph(Surv(css36_time, css36_status) ~ CMETS_z, data = css_dat),
  "CSS48: CMETS" = coxph(Surv(css48_time, css48_status) ~ CMETS_z, data = css_dat),
  "CSS24: CMETS + EMT + Hypoxia" = coxph(
    Surv(css24_time, css24_status) ~ CMETS_z + EMT_z + Hypoxia_z,
    data = css_dat
  ),
  "CSS24: CMETS + tx" = coxph(
    Surv(css24_time, css24_status) ~ CMETS_z + tx_binary,
    data = css_dat_tx
  ),
  "CSS36: CMETS + EMT + Hypoxia" = coxph(
    Surv(css36_time, css36_status) ~ CMETS_z + EMT_z + Hypoxia_z,
    data = css_dat
  ),
  "CSS36: CMETS + tx" = coxph(
    Surv(css36_time, css36_status) ~ CMETS_z + tx_binary,
    data = css_dat_tx
  ),
  "CSS48: CMETS + EMT + Hypoxia" = coxph(
    Surv(css48_time, css48_status) ~ CMETS_z + EMT_z + Hypoxia_z,
    data = css_dat
  ),
  "CSS48: CMETS + tx" = coxph(
    Surv(css48_time, css48_status) ~ CMETS_z + tx_binary,
    data = css_dat_tx
  ),
  "CSS60: Stem" = coxph(Surv(css60_time, css60_status) ~ Stem_z, data = css_dat),
  "CSS60: EMT" = coxph(Surv(css60_time, css60_status) ~ EMT_z, data = css_dat),
  "CSS60: Hypoxia" = coxph(Surv(css60_time, css60_status) ~ Hypoxia_z, data = css_dat),
  "CSS60: CMETS + Stem" = coxph(Surv(css60_time, css60_status) ~ CMETS_z + Stem_z, data = css_dat),
  "CSS60: CMETS + EMT" = coxph(Surv(css60_time, css60_status) ~ CMETS_z + EMT_z, data = css_dat),
  "CSS60: CMETS + Hypoxia" = coxph(Surv(css60_time, css60_status) ~ CMETS_z + Hypoxia_z, data = css_dat),
  "CSS60: CMETS + EMT + Hypoxia" = coxph(
    Surv(css60_time, css60_status) ~ CMETS_z + EMT_z + Hypoxia_z,
    data = css_dat
  ),
  "CSS60: CMETS + Stem + EMT + Hypoxia" = coxph(
    Surv(css60_time, css60_status) ~ CMETS_z + Stem_z + EMT_z + Hypoxia_z,
    data = css_dat
  ),
  "CSS60: CMETS + tx" = coxph(Surv(css60_time, css60_status) ~ CMETS_z + tx_binary, data = css_dat_tx),
  "CSS60: CMETS + Stem + tx" = coxph(
    Surv(css60_time, css60_status) ~ CMETS_z + Stem_z + tx_binary,
    data = css_dat_tx
  ),
  "CSS60: CMETS + Stem + EMT + Hypoxia + tx" = coxph(
    Surv(css60_time, css60_status) ~ CMETS_z + Stem_z + EMT_z + Hypoxia_z + tx_binary,
    data = css_dat_tx
  )
)

summary_all_terms <- bind_rows(
  build_result_tbl(fits_tcga, "TCGA", "OS"),
  build_result_tbl(fits_39582, "GSE39582", "OS"),
  build_result_tbl(fits_17536, "GSE17536", "OS"),
  build_result_tbl(fits_41258, "GSE41258", "OS"),
  build_result_tbl(fits_159216, "GSE159216", "CSS")
)

summary_msi_mmr <- bind_rows(
  if (length(fits_tcga_msi)) build_result_tbl(fits_tcga_msi, "TCGA", "OS") else NULL,
  if (length(fits_39582_mmr)) build_result_tbl(fits_39582_mmr, "GSE39582", "OS") else NULL,
  if (length(fits_41258_msi)) build_result_tbl(fits_41258_msi, "GSE41258", "OS") else NULL
)
if (is.null(summary_msi_mmr) || !is.data.frame(summary_msi_mmr)) {
  summary_msi_mmr <- tibble()
}

summary_cmets <- summary_all_terms %>%
  filter(term == "CMETS_z") %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
    p_fmt = ifelse(p_value < 0.001, formatC(p_value, format = "e", digits = 2), sprintf("%.3f", p_value))
  ) %>%
  select(cohort, endpoint, model, HR_CI, p_fmt, n, events) %>%
  arrange(cohort, model)

write.table(
  summary_all_terms,
  file = file.path(output_dir, "summary_all_terms_aucell_clean.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  summary_cmets,
  file = file.path(output_dir, "summary_cmets_aucell_clean.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (nrow(summary_msi_mmr) > 0) {
  write.table(
    summary_msi_mmr,
    file = file.path(output_dir, "summary_cmets_aucell_msi_mmr.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

message("Done.")
message("Wrote: ", file.path(output_dir, "summary_all_terms_aucell_clean.tsv"))
message("Wrote: ", file.path(output_dir, "summary_cmets_aucell_clean.tsv"))
if (nrow(summary_msi_mmr) > 0) {
  message("Wrote: ", file.path(output_dir, "summary_cmets_aucell_msi_mmr.tsv"))
}
