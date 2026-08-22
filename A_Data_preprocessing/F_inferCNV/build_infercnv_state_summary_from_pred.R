#!/usr/bin/env Rscript
# Build per-sample InferCNV "state summary" tables analogous to
#   /diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20220421_InferCNV_state_summary/
# from 20241001 HMM Leiden subcluster outputs:
#   <sample>/out/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat
#   <sample>/out/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings
#
# Outputs per sample (under OUT_ROOT/<sample_id>/):
#   InferCNV_states.tsv     — genes × columns {sample}_Tumor_s* (HMM state integers)
#   InferCNV_labels.tsv   — barcode, seurat_clusters, InferCNV_state ({sample}_Tumor_s* or NA)
#   sn_CNV_cluster_count.tsv — contingency: seurat_clusters (rows) × InferCNV_state (cols); tumor-assigned cells only
#
# Run:
#   Rscript build_infercnv_state_summary_from_pred.R
# Or set env: SAMPLES, INFERCNV_ROOT, SEURAT_RDS, OUT_ROOT
# Per-sample inferCNV run directory (optional):
#   INFERCNV_ROOT_OVERRIDES="HT225C1-Th1=/path/to/20241002_inferCNV_results"
#   (comma-separated list; paths are the parent folder that contains <sample_id>/out/...)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

PRED_BASENAME <- paste0(
  "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.",
  "Pnorm_0.5.pred_cnv_genes.dat"
)
CG_BASENAME <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"
# HMM state when a tumor subcluster (or gene) has no row in pred (inferCNV neutral).
NEUTRAL_HMM_STATE <- as.integer(Sys.getenv("INFERCNV_NEUTRAL_HMM_STATE", unset = "3"))

samples_env <- Sys.getenv("SAMPLES", unset = "")
SAMPLES <- if (nzchar(samples_env)) {
  trimws(strsplit(samples_env, ",", fixed = TRUE)[[1L]])
} else {
  c(
    "CM329C1-S1",
    "CM349C1-S1",
    "HT225C1-Th1",
    "HT230C1-Th1",
    "HT260C1-Th1",
    "CM426C1-Th1"
  )
}

INFERCNV_ROOT <- Sys.getenv(
  "INFERCNV_ROOT",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20241001_inferCNV_results"
)

parse_infercnv_root_overrides <- function(raw) {
  if (!nzchar(raw)) {
    return(setNames(character(), character()))
  }
  out <- character()
  for (p in trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])) {
    idx <- regexpr("=", p, fixed = TRUE)
    if (idx < 1L) next
    sid <- trimws(substring(p, 1L, idx - 1L))
    root <- trimws(substring(p, idx + 1L, nchar(p)))
    if (nzchar(sid) && nzchar(root)) out[sid] <- root
  }
  out
}

INFERCNV_ROOT_OVERRIDE_MAP <- parse_infercnv_root_overrides(
  Sys.getenv("INFERCNV_ROOT_OVERRIDES", unset = "")
)

infercnv_results_root <- function(sample_id) {
  if (sample_id %in% names(INFERCNV_ROOT_OVERRIDE_MAP)) {
    INFERCNV_ROOT_OVERRIDE_MAP[[sample_id]]
  } else {
    INFERCNV_ROOT
  }
}

SEURAT_RDS <- Sys.getenv(
  "SEURAT_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset.rds"
)

OUT_ROOT <- Sys.getenv(
  "OUT_ROOT",
  unset = file.path(getwd(), "infercnv_state_summary_fisher5_20241001")
)

infercnv_state_colname <- function(sample_id, cell_group_name) {
  if (!nzchar(cell_group_name) || !grepl("^Tumor\\.", cell_group_name)) {
    return(NA_character_)
  }
  sub <- sub("^Tumor\\.", "", cell_group_name)
  paste0(sample_id, "_", sub)
}

process_one_sample <- function(sample_id, so) {
  root <- infercnv_results_root(sample_id)
  pred_path <- file.path(root, sample_id, "out", PRED_BASENAME)
  cg_path <- file.path(root, sample_id, "out", CG_BASENAME)

  if (!file.exists(pred_path)) {
    warning("Missing pred file: ", pred_path)
    return(invisible(NULL))
  }
  if (!file.exists(cg_path)) {
    warning("Missing cell_groupings: ", cg_path)
    return(invisible(NULL))
  }

  pred <- read.delim(pred_path, stringsAsFactors = FALSE, check.names = FALSE)
  cg <- read.delim(cg_path, stringsAsFactors = FALSE, check.names = FALSE)

  reqp <- c("cell_group_name", "gene", "state")
  if (!all(reqp %in% names(pred))) {
    stop("Unexpected pred columns in ", pred_path)
  }
  pred$state <- suppressWarnings(as.integer(pred$state))

  tum <- pred %>%
    dplyr::filter(grepl("^Tumor\\.", .data$cell_group_name)) %>%
    dplyr::mutate(
      colname = vapply(
        .data$cell_group_name,
        function(g) infercnv_state_colname(sample_id, g),
        character(1L)
      )
    ) %>%
    dplyr::filter(!is.na(.data$colname))

  ## Genes: all genes seen for any Tumor.* row in pred (wide matrix rows).
  gene_universe <- if (nrow(tum)) {
    unique(tum$gene[!is.na(tum$gene)])
  } else {
    warning("No Tumor.* rows in pred for ", sample_id, "; using genes from full pred table")
    unique(pred$gene[!is.na(pred$gene)])
  }

  ## Expected tumor columns from cell_groupings (subclusters assigned to cells).
  cg_tumor_groups <- unique(cg$cell_group_name[grepl("^Tumor\\.", cg$cell_group_name)])
  expected_colnames <- unique(stats::na.omit(vapply(
    cg_tumor_groups,
    function(g) infercnv_state_colname(sample_id, g),
    character(1L)
  )))

  if (!length(expected_colnames) && !nrow(tum)) {
    warning("No tumor cell groups in cg and no Tumor.* pred rows for ", sample_id)
    return(invisible(NULL))
  }

  if (nrow(tum)) {
    states_wide <- tum %>%
      dplyr::select("gene", "colname", "state") %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = "colname", values_from = "state") %>%
      dplyr::rename(Gene = "gene")
  } else {
    states_wide <- tibble::tibble(Gene = gene_universe)
  }

  ## Subclusters in cell_groupings but absent from pred → neutral at all genes.
  missing_cols <- setdiff(expected_colnames, setdiff(names(states_wide), "Gene"))
  for (cl in missing_cols) {
    states_wide[[cl]] <- NEUTRAL_HMM_STATE
  }

  ## Subcluster present but gene missing in pred for that column → neutral.
  for (cl in intersect(expected_colnames, setdiff(names(states_wide), "Gene"))) {
    idx_na <- is.na(states_wide[[cl]])
    if (any(idx_na)) {
      states_wide[[cl]][idx_na] <- NEUTRAL_HMM_STATE
    }
  }

  if (!nrow(states_wide)) {
    warning("No genes in state matrix for ", sample_id)
    return(invisible(NULL))
  }

  ## column order: Tumor_s1, Tumor_s2, …
  cn <- setdiff(names(states_wide), "Gene")
  if (!length(cn)) {
    warning("No tumor state columns for ", sample_id)
    return(invisible(NULL))
  }
  tum_sub <- sub(paste0("^", sample_id, "_"), "", cn)
  ord <- order(as.integer(sub("^Tumor_s", "", tum_sub)))
  states_wide <- states_wide[, c("Gene", cn[ord]), drop = FALSE]

  out_dir <- file.path(OUT_ROOT, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_tsv(states_wide, file.path(out_dir, "InferCNV_states.tsv"), na = "")

  ## Labels + Seurat clusters
  if (!"cell" %in% names(cg) || !"cell_group_name" %in% names(cg)) {
    stop("Unexpected cell_groupings columns in ", cg_path)
  }

  labels <- tibble::tibble(
    barcode = cg$cell,
    infercnv_cell_group = cg$cell_group_name,
    InferCNV_state = vapply(
      cg$cell_group_name,
      function(g) infercnv_state_colname(sample_id, g),
      character(1L)
    )
  )

  md <- tibble::as_tibble(so@meta.data, rownames = "cell")
  if (!all(labels$barcode %in% md$cell)) {
    n_miss <- sum(!labels$barcode %in% md$cell)
    warning(sample_id, ": ", n_miss, " cell_groupings barcodes not in Seurat object")
  }

  if (!"seurat_clusters" %in% names(md)) {
    stop("Seurat metadata missing seurat_clusters")
  }

  labels <- labels %>%
    dplyr::left_join(
      md %>% dplyr::select("cell", "seurat_clusters"),
      by = c("barcode" = "cell")
    )

  labels_out <- labels %>%
    dplyr::transmute(
      barcode = .data$barcode,
      seurat_clusters = as.character(.data$seurat_clusters),
      InferCNV_state = dplyr::if_else(
        is.na(.data$InferCNV_state) | .data$InferCNV_state == "",
        NA_character_,
        .data$InferCNV_state
      )
    )

  readr::write_tsv(labels_out, file.path(out_dir, "InferCNV_labels.tsv"), na = "")

  ## Contingency: only cells with a tumor InferCNV_state (matches state matrix columns)
  lab_tumor <- labels %>%
    dplyr::filter(!is.na(.data$InferCNV_state), nzchar(.data$InferCNV_state))

  if (nrow(lab_tumor)) {
    tab <- table(
      seurat_cluster = as.character(lab_tumor$seurat_clusters),
      InferCNV_state = lab_tumor$InferCNV_state
    )
    ct <- as.data.frame.matrix(tab)
    ct <- tibble::rownames_to_column(ct, "seurat_clusters")
    readr::write_tsv(ct, file.path(out_dir, "sn_CNV_cluster_count.tsv"), na = "0")
  } else {
    message("No tumor-labeled cells for contingency: ", sample_id)
  }

  message("Wrote state summary for ", sample_id, " -> ", out_dir)
  invisible(NULL)
}

if (!file.exists(SEURAT_RDS)) {
  stop("Seurat RDS not found: ", SEURAT_RDS, call. = FALSE)
}

so <- readRDS(SEURAT_RDS)

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
message("OUT_ROOT: ", OUT_ROOT)
message("INFERCNV_ROOT: ", INFERCNV_ROOT)
if (length(INFERCNV_ROOT_OVERRIDE_MAP)) {
  message(
    "INFERCNV_ROOT_OVERRIDES: ",
    paste(names(INFERCNV_ROOT_OVERRIDE_MAP), INFERCNV_ROOT_OVERRIDE_MAP, sep = " -> ", collapse = "; ")
  )
}
message("Samples: ", paste(SAMPLES, collapse = ", "))

for (sid in SAMPLES) {
  process_one_sample(sid, so)
}

message("Done.")
