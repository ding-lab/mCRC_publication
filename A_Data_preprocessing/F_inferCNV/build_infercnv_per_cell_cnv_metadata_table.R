#!/usr/bin/env Rscript
# Build a two-column table (cell ID + InferCNV CNV call) from
#   20241001_inferCNV_state_summary_all_groups/<sample>/
#     InferCNV_labels_all_groups.tsv
#     InferCNV_states_all_groups.tsv
#
# Rule (per cell, using the gene × state column for that cell's InferCNV group):
#   If any gene has HMM state in {1, 2, 4, 5, 6}  ->  "CNV_change"
#   If all genes are 3 or NA                       ->  "Neutral"
#   Cells with no InferCNV_state in labels         ->  NA (column value literal "NA" or empty)
#   Other / unknown integer states                 ->  NA
#
# Output (default): STATE_SUMMARY_ROOT/InferCNV_per_cell_cnv_call.tsv
#   cell_id, infercnv_cnv_call
#
# Add to Seurat:
#   m <- read.delim("InferCNV_per_cell_cnv_call.tsv", row.names = 1)
#   obj <- AddMetaData(obj, m)   # or match by cell name / colnames
#
# Run:
#   Rscript build_infercnv_per_cell_cnv_metadata_table.R
# Env: STATE_SUMMARY_ROOT, OUT_TSV, SAMPLES (optional comma-list; default = all sample subdirs)
# 20241002 cohort summaries:
#   bash build_infercnv_per_cell_cnv_call_20241002.sh
#   # or: STATE_SUMMARY_ROOT=.../20241002_inferCNV_state_summary_all_groups Rscript ...
# 20260403 cohort summaries:
#   bash build_infercnv_per_cell_cnv_call_20260403.sh
#   # or: STATE_SUMMARY_ROOT=.../20260403_inferCNV_state_summary_all_groups Rscript ...

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

CNV_STATES <- c(1L, 2L, 4L, 5L, 6L)
NEUTRAL_STATES <- 3L

STATE_SUMMARY_ROOT <- Sys.getenv(
  "STATE_SUMMARY_ROOT",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20241001_inferCNV_state_summary_all_groups"
)

OUT_TSV <- Sys.getenv(
  "OUT_TSV",
  unset = file.path(STATE_SUMMARY_ROOT, "InferCNV_per_cell_cnv_call.tsv")
)

LABELS_FILE <- "InferCNV_labels_all_groups.tsv"
STATES_FILE <- "InferCNV_states_all_groups.tsv"

discover_samples <- function(root) {
  if (!dir.exists(root)) {
    stop("STATE_SUMMARY_ROOT does not exist: ", root, call. = FALSE)
  }
  dirs <- list.dirs(root, full.names = FALSE, recursive = FALSE)
  dirs <- dirs[nzchar(dirs)]
  ok <- vapply(dirs, function(sid) {
    file.exists(file.path(root, sid, LABELS_FILE)) &&
      file.exists(file.path(root, sid, STATES_FILE))
  }, logical(1L))
  sort(dirs[ok])
}

samples_env <- Sys.getenv("SAMPLES", unset = "")
SAMPLES <- if (nzchar(samples_env)) {
  trimws(strsplit(samples_env, ",", fixed = TRUE)[[1L]])
} else {
  discover_samples(STATE_SUMMARY_ROOT)
}

#' Classify one inferCNV group column in the wide states table.
classify_infercnv_column <- function(states_df, colname) {
  if (!nzchar(colname) || !colname %in% names(states_df)) {
    return(NA_character_)
  }
  v <- suppressWarnings(as.integer(states_df[[colname]]))
  if (any(v %in% CNV_STATES, na.rm = TRUE)) {
    return("CNV_change")
  }
  if (all(is.na(v) | v == NEUTRAL_STATES)) {
    return("Neutral")
  }
  NA_character_
}

process_one_sample <- function(sample_id) {
  lab_path <- file.path(STATE_SUMMARY_ROOT, sample_id, LABELS_FILE)
  st_path <- file.path(STATE_SUMMARY_ROOT, sample_id, STATES_FILE)
  labels <- readr::read_tsv(lab_path, show_col_types = FALSE)
  states <- readr::read_tsv(st_path, show_col_types = FALSE)

  if (!"Gene" %in% names(states)) {
    stop("Missing Gene column in ", st_path, call. = FALSE)
  }
  if (!all(c("barcode", "InferCNV_state") %in% names(labels))) {
    stop("Unexpected columns in ", lab_path, call. = FALSE)
  }

  grp <- unique(labels$InferCNV_state)
  grp <- grp[!is.na(grp) & nzchar(as.character(grp))]
  call_map <- vapply(grp, function(g) classify_infercnv_column(states, g), character(1L))

  out <- labels %>%
    dplyr::transmute(
      cell_id = .data$barcode,
      infercnv_cnv_call = dplyr::case_when(
        is.na(.data$InferCNV_state) | .data$InferCNV_state == "" ~ NA_character_,
        TRUE ~ unname(call_map[as.character(.data$InferCNV_state)])
      )
    )
  out
}

if (!length(SAMPLES)) {
  stop("No samples found under ", STATE_SUMMARY_ROOT, call. = FALSE)
}

message("STATE_SUMMARY_ROOT: ", STATE_SUMMARY_ROOT)
message("Samples: ", length(SAMPLES))

all_rows <- dplyr::bind_rows(lapply(SAMPLES, process_one_sample))

# Single combined table for merged Seurat (cell barcodes should be unique)
if (any(duplicated(all_rows$cell_id))) {
  warning(
    "Duplicate cell_id in output (first wins): ",
    sum(duplicated(all_rows$cell_id)),
    " duplicates"
  )
  all_rows <- all_rows %>% dplyr::distinct(.data$cell_id, .keep_all = TRUE)
}

readr::write_tsv(all_rows, OUT_TSV, na = "")

message("Wrote ", nrow(all_rows), " rows -> ", OUT_TSV)
message("Done.")
