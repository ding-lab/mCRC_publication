#!/usr/bin/env Rscript
# InferCNV HMM tumor gain calls per gene (pred_cnv_genes.dat) × sample.
# Subclone = cell group "Tumor.Tumor_s*" (Leiden subcluster within Tumor).
# Default genes: KRAS + MAPK/RTK-related panel (PIK3CA included per user list).
# All genes: gain if HMM state >= INFERCNV_GAIN_MIN_STATE (default 4).
#
# Output:
#   Primary CSV (long): one row per sample × gene — gain_in_tumor, subclones, cell proportions, etc.
#   Legacy CSV: KRAS-only rows with original column names (kras_*) for existing notebooks.
#
# Run: Rscript infercnv_kras_gain_tumor_table.R
# Or:  INFERCNV_GENES="KRAS,EGFR" INFERCNV_GAIN_MIN_STATE=4 OUTPUT_CSV=out.csv Rscript ...

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Install dplyr: install.packages(\"dplyr\")", call. = FALSE)
  }
})

## ---- user parameters ---------------------------------------------------------
INFERCNV_RESULTS_DIR <- Sys.getenv(
  "INFERCNV_RESULTS_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20241001_inferCNV_results"
)

PRED_FILE <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_genes.dat"
CELL_GROUPINGS_FILE <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"

# Gain (all genes): pred HMM state >= GAIN_MIN_STATE (override with INFERCNV_GAIN_MIN_STATE).
GAIN_MIN_STATE <- suppressWarnings(as.integer(Sys.getenv("INFERCNV_GAIN_MIN_STATE", unset = "4")))
if (!is.finite(GAIN_MIN_STATE)) {
  stop("INFERCNV_GAIN_MIN_STATE must be a single integer (default 4)", call. = FALSE)
}

# Genes to summarize (HGNC symbols as in pred_cnv_genes.dat). Override: INFERCNV_GENES=KRAS,EGFR,...
GENES_DEFAULT <- c(
  "KRAS",
  "AURKA",
  "EGFR",
  "MYC",
  "CDK8",
  "CDX2",
  "ERBB2",
  "PIK3CA",
  "ROS1"
)
genes_env <- Sys.getenv("INFERCNV_GENES", unset = "")
GENES_OF_INTEREST <- if (nzchar(genes_env)) {
  trimws(strsplit(genes_env, ",", fixed = TRUE)[[1L]])
} else {
  GENES_DEFAULT
}
GENES_OF_INTEREST <- unique(GENES_OF_INTEREST[nzchar(GENES_OF_INTEREST)])
if (!length(GENES_OF_INTEREST)) {
  stop("No genes in INFERCNV_GENES / default list", call. = FALSE)
}

OUTPUT_CSV <- Sys.getenv(
  "OUTPUT_CSV",
  unset = file.path(getwd(), "infercnv_gene_gain_tumor_by_sample.csv")
)

# Same schema as older script — KRAS rows only, for notebook merge scripts.
OUTPUT_CSV_KRAS_LEGACY <- Sys.getenv(
  "OUTPUT_CSV_KRAS_LEGACY",
  unset = file.path(dirname(OUTPUT_CSV), "kras_gain_tumor_by_sample.csv")
)

## ---- helpers -----------------------------------------------------------------
sample_id_from_path <- function(path) {
  d <- dirname(path)
  basename(dirname(d))
}

read_pred_cnv_genes <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
}

read_cell_groupings <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = ""
  )
}

tumor_gene_gain_cell_props <- function(cg, tumor_groups_with_gain) {
  empty_counts <- function() {
    list(
      n_tumor_cells_total = NA_integer_,
      n_tumor_cells_gain = NA_integer_,
      prop_tumor_cells_gain = NA_real_
    )
  }
  if (is.null(cg) || !nrow(cg)) {
    return(empty_counts())
  }
  req <- c("cell_group_name", "cell")
  if (!all(req %in% names(cg))) {
    return(empty_counts())
  }
  tum <- cg[grepl("^Tumor\\.", cg$cell_group_name), , drop = FALSE]
  n_tot <- nrow(tum)
  if (n_tot == 0L) {
    return(list(
      n_tumor_cells_total = 0L,
      n_tumor_cells_gain = 0L,
      prop_tumor_cells_gain = NA_real_
    ))
  }
  grp_gain <- unique(tumor_groups_with_gain)
  n_gain <- sum(tum$cell_group_name %in% grp_gain)
  list(
    n_tumor_cells_total = as.integer(n_tot),
    n_tumor_cells_gain = as.integer(n_gain),
    prop_tumor_cells_gain = n_gain / n_tot
  )
}

tumor_subclone_label <- function(cell_group_name) {
  sub <- sub("^Tumor\\.", "", cell_group_name)
  ifelse(grepl("^Tumor_s", sub), sub, cell_group_name)
}

append_cell_prop_cols <- function(df, props, notes_add = NULL) {
  df$n_tumor_cells_total <- props$n_tumor_cells_total
  df$n_tumor_cells_gain <- props$n_tumor_cells_gain
  df$prop_tumor_cells_gain <- props$prop_tumor_cells_gain
  if (length(notes_add) && nzchar(notes_add)) {
    df$notes <- if (nzchar(df$notes)) paste(df$notes, notes_add, sep = ";") else notes_add
  }
  df
}

empty_gene_row <- function(sid, gene, notes_txt, pred_ok, n_tumor_total, cg) {
  props <- tumor_gene_gain_cell_props(cg, character(0))
  note2 <- if (is.null(cg) || !nrow(cg)) "missing_cell_groupings" else ""
  df <- data.frame(
    sample = sid,
    gene = gene,
    pred_file_found = pred_ok,
    n_tumor_subclones_total = n_tumor_total,
    gain_in_tumor = NA,
    tumor_subclones_with_gain = NA_character_,
    n_subclones_with_gain = NA_integer_,
    gain_is_subclonal = NA,
    gain_max_state = NA_integer_,
    notes = notes_txt,
    stringsAsFactors = FALSE
  )
  append_cell_prop_cols(df, props, note2)
}

summarize_one_gene <- function(sid, gene, cg, tab, n_tumor_total) {
  note_cg <- if (is.null(cg) || !nrow(cg)) "missing_cell_groupings" else ""

  k <- tab[tab$gene == gene & grepl("^Tumor\\.", tab$cell_group_name), , drop = FALSE]

  if (!nrow(k)) {
    props <- tumor_gene_gain_cell_props(cg, character(0))
    df <- data.frame(
      sample = sid,
      gene = gene,
      pred_file_found = TRUE,
      n_tumor_subclones_total = n_tumor_total,
      gain_in_tumor = FALSE,
      tumor_subclones_with_gain = NA_character_,
      n_subclones_with_gain = 0L,
      gain_is_subclonal = FALSE,
      gain_max_state = NA_integer_,
      notes = "no_gene_row_in_tumor_groups",
      stringsAsFactors = FALSE
    )
    return(append_cell_prop_cols(df, props, note_cg))
  }

  k_gain <- k[!is.na(k$state) & k$state >= GAIN_MIN_STATE, , drop = FALSE]
  groups_gain_full <- unique(k_gain$cell_group_name)
  labels_gain <- sort(unique(tumor_subclone_label(k_gain$cell_group_name)))

  max_state <- max(k$state, na.rm = TRUE)
  if (!is.finite(max_state)) {
    max_state <- NA_integer_
  }

  n_gain <- length(labels_gain)
  has_gain <- n_gain > 0L
  subclonal <- isTRUE(n_tumor_total > 1L && has_gain && n_gain < n_tumor_total)

  props <- tumor_gene_gain_cell_props(cg, groups_gain_full)

  df <- data.frame(
    sample = sid,
    gene = gene,
    pred_file_found = TRUE,
    n_tumor_subclones_total = n_tumor_total,
    gain_in_tumor = has_gain,
    tumor_subclones_with_gain = if (n_gain) paste(labels_gain, collapse = ";") else NA_character_,
    n_subclones_with_gain = as.integer(n_gain),
    gain_is_subclonal = subclonal,
    gain_max_state = as.integer(max_state),
    notes = "",
    stringsAsFactors = FALSE
  )
  append_cell_prop_cols(df, props, note_cg)
}

summarize_one_sample_path <- function(path) {
  sid <- sample_id_from_path(path)
  cg_path <- file.path(dirname(path), CELL_GROUPINGS_FILE)
  cg <- read_cell_groupings(cg_path)
  tab <- read_pred_cnv_genes(path)

  if (is.null(tab) || !nrow(tab)) {
    return(dplyr::bind_rows(lapply(GENES_OF_INTEREST, function(g) {
      empty_gene_row(sid, g, "missing_or_empty_pred_file", FALSE, NA_integer_, cg)
    })))
  }

  req <- c("cell_group_name", "gene_region_name", "state", "gene")
  if (!all(req %in% names(tab))) {
    return(dplyr::bind_rows(lapply(GENES_OF_INTEREST, function(g) {
      empty_gene_row(sid, g, "unexpected_column_names", TRUE, NA_integer_, cg)
    })))
  }

  tab$state <- suppressWarnings(as.integer(tab$state))
  tumor_groups <- unique(tab$cell_group_name[grepl("^Tumor\\.", tab$cell_group_name)])
  n_tumor_total <- length(tumor_groups)

  dplyr::bind_rows(lapply(GENES_OF_INTEREST, function(g) {
    summarize_one_gene(sid, g, cg, tab, n_tumor_total)
  }))
}

missing_pred_gene_rows <- function(sid, cg) {
  dplyr::bind_rows(lapply(GENES_OF_INTEREST, function(g) {
    props <- tumor_gene_gain_cell_props(cg, character(0))
    note <- "missing_pred_file"
    if (is.null(cg) || !nrow(cg)) {
      note <- paste(note, "missing_cell_groupings", sep = ";")
    }
    df <- data.frame(
      sample = sid,
      gene = g,
      pred_file_found = FALSE,
      n_tumor_subclones_total = NA_integer_,
      gain_in_tumor = NA,
      tumor_subclones_with_gain = NA_character_,
      n_subclones_with_gain = NA_integer_,
      gain_is_subclonal = NA,
      gain_max_state = NA_integer_,
      notes = note,
      stringsAsFactors = FALSE
    )
    append_cell_prop_cols(df, props, NULL)
  }))
}

to_kras_legacy <- function(long_df) {
  k <- long_df[long_df$gene == "KRAS", , drop = FALSE]
  if (!nrow(k)) {
    return(k)
  }
  data.frame(
    sample = k$sample,
    pred_file_found = k$pred_file_found,
    n_tumor_subclones_total = k$n_tumor_subclones_total,
    kras_gain_in_tumor = k$gain_in_tumor,
    tumor_subclones_with_kras_gain = k$tumor_subclones_with_gain,
    n_tumor_subclones_with_kras_gain = k$n_subclones_with_gain,
    kras_gain_is_subclonal = k$gain_is_subclonal,
    kras_gain_max_state = k$gain_max_state,
    notes = k$notes,
    n_tumor_cells_total = k$n_tumor_cells_total,
    n_tumor_cells_kras_gain = k$n_tumor_cells_gain,
    prop_tumor_cells_kras_gain = k$prop_tumor_cells_gain,
    stringsAsFactors = FALSE
  )
}

## ---- main --------------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))

paths <- sort(Sys.glob(file.path(INFERCNV_RESULTS_DIR, "*", "out", PRED_FILE)))

if (!length(paths)) {
  stop("No pred_cnv_genes files found under: ", INFERCNV_RESULTS_DIR)
}

out <- dplyr::bind_rows(lapply(paths, summarize_one_sample_path))

sample_dirs <- Sys.glob(file.path(INFERCNV_RESULTS_DIR, "*"))
sample_dirs <- sample_dirs[dir.exists(sample_dirs)]
all_ids <- basename(sample_dirs)
missing_pred <- setdiff(all_ids, unique(out$sample))

if (length(missing_pred)) {
  extra_rows <- lapply(missing_pred, function(sid) {
    cg_path <- file.path(INFERCNV_RESULTS_DIR, sid, "out", CELL_GROUPINGS_FILE)
    cg <- read_cell_groupings(cg_path)
    missing_pred_gene_rows(sid, cg)
  })
  out <- dplyr::bind_rows(out, dplyr::bind_rows(extra_rows))
}

out <- out %>%
  dplyr::arrange(sample, gene) %>%
  dplyr::mutate(
    tumor_subclones_with_gain = dplyr::if_else(
      is.na(tumor_subclones_with_gain) | tumor_subclones_with_gain == "",
      NA_character_,
      tumor_subclones_with_gain
    )
  )

utils::write.csv(out, OUTPUT_CSV, row.names = FALSE, na = "")
message("Wrote (long, all genes): ", OUTPUT_CSV)

kras_legacy <- to_kras_legacy(out)
utils::write.csv(kras_legacy, OUTPUT_CSV_KRAS_LEGACY, row.names = FALSE, na = "")
message("Wrote (KRAS legacy): ", OUTPUT_CSV_KRAS_LEGACY)

message(
  "Rows: ", nrow(out), " (samples × genes) | genes: ",
  paste(GENES_OF_INTEREST, collapse = ","),
  " | gain: state>=", GAIN_MIN_STATE
)
