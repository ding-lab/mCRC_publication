#!/usr/bin/env Rscript
# Like build_infercnv_state_summary_from_pred_all_samples.R but includes *all*
# InferCNV cell_group_name entries (tumor + reference / immune / stromal Leiden groups),
# not only Tumor.* rows.
#
# pred_cnv_genes.dat often has states only for some groups (e.g. tumor + one reference);
# any group listed in cell_groupings but missing from pred gets NEUTRAL_HMM_STATE at all genes.
#
# Outputs per sample (under OUT_ROOT/<sample_id>/):
#   InferCNV_states_all_groups.tsv
#   InferCNV_labels_all_groups.tsv
#   sn_CNV_cluster_count_all_groups.tsv
#
# Run:
#   Rscript build_infercnv_state_summary_from_pred_all_groups.R
# Env: SAMPLES, INFERCNV_ROOT, SEURAT_RDS, OUT_ROOT, INFERCNV_NEUTRAL_HMM_STATE (same as all_samples)
#   INFERCNV_ROOT_OVERRIDES="HT225C1-Th1=/diskmnt/.../20241002_inferCNV_results"  (comma-separated)

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
NEUTRAL_HMM_STATE <- as.integer(Sys.getenv("INFERCNV_NEUTRAL_HMM_STATE", unset = "3"))

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
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20241001_inferCNV_state_summary_all_groups"
)

discover_infercnv_samples <- function(root) {
  if (!dir.exists(root)) {
    stop("INFERCNV_ROOT does not exist: ", root, call. = FALSE)
  }
  cand <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  ok <- vapply(cand, function(d) {
    file.exists(file.path(d, "out", PRED_BASENAME)) &&
      file.exists(file.path(d, "out", CG_BASENAME))
  }, logical(1L))
  sort(basename(cand)[ok])
}

samples_env <- Sys.getenv("SAMPLES", unset = "")
SAMPLES <- if (nzchar(samples_env)) {
  trimws(strsplit(samples_env, ",", fixed = TRUE)[[1L]])
} else {
  sort(unique(c(
    discover_infercnv_samples(INFERCNV_ROOT),
    names(INFERCNV_ROOT_OVERRIDE_MAP)
  )))
}

#' Column name for wide state matrix / InferCNV_state label (every assigned group).
infercnv_state_colname_all_groups <- function(sample_id, cell_group_name) {
  if (!nzchar(cell_group_name)) {
    return(NA_character_)
  }
  if (grepl("^Tumor\\.", cell_group_name)) {
    paste0(sample_id, "_", sub("^Tumor\\.", "", cell_group_name))
  } else {
    sfx <- gsub("[^A-Za-z0-9]+", "_", cell_group_name)
    sfx <- gsub("_+", "_", sfx)
    sfx <- sub("^_|_$", "", sfx)
    paste0(sample_id, "_", sfx)
  }
}

#' Order columns: Tumor_s* numerically, then other groups alphabetically.
order_state_colnames <- function(sample_id, cn) {
  pfx <- paste0(sample_id, "_")
  npx <- nchar(pfx)
  is_tumor_col <- vapply(cn, function(c) {
    startsWith(c, pfx) &&
      grepl("^Tumor_s[0-9]+$", substring(c, npx + 1L), perl = TRUE)
  }, logical(1L))
  cn_tum <- cn[is_tumor_col]
  cn_oth <- sort(cn[!is_tumor_col])
  if (length(cn_tum)) {
    rest <- substring(cn_tum, npx + 1L)
    nums <- as.integer(sub("^Tumor_s", "", rest))
    cn_tum <- cn_tum[order(nums)]
  }
  c(cn_tum, cn_oth)
}

process_one_sample_all_groups <- function(sample_id, so) {
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

  pred_use <- pred %>%
    dplyr::filter(!is.na(.data$cell_group_name), nzchar(.data$cell_group_name)) %>%
    dplyr::mutate(
      colname = vapply(
        .data$cell_group_name,
        function(g) infercnv_state_colname_all_groups(sample_id, g),
        character(1L)
      )
    ) %>%
    dplyr::filter(!is.na(.data$colname))

  dup_map <- pred_use %>%
    dplyr::distinct(.data$cell_group_name, .data$colname)
  if (any(duplicated(dup_map$colname))) {
    stop(
      "Duplicate state column names after sanitization for ",
      sample_id, ": ",
      paste(dup_map$colname[duplicated(dup_map$colname)], collapse = ", "),
      call. = FALSE
    )
  }

  gene_universe <- if (nrow(pred_use)) {
    unique(pred_use$gene[!is.na(pred_use$gene)])
  } else {
    warning("No usable pred rows for ", sample_id, "; using genes from full pred table")
    unique(pred$gene[!is.na(pred$gene)])
  }

  cg_groups <- unique(cg$cell_group_name[!is.na(cg$cell_group_name) & nzchar(cg$cell_group_name)])
  expected_colnames <- unique(stats::na.omit(vapply(
    cg_groups,
    function(g) infercnv_state_colname_all_groups(sample_id, g),
    character(1L)
  )))

  if (!length(expected_colnames) && !nrow(pred_use)) {
    warning("No cell groups in cg and no pred rows for ", sample_id)
    return(invisible(NULL))
  }

  if (nrow(pred_use)) {
    states_wide <- pred_use %>%
      dplyr::select("gene", "colname", "state") %>%
      dplyr::distinct() %>%
      tidyr::pivot_wider(names_from = "colname", values_from = "state") %>%
      dplyr::rename(Gene = "gene")
  } else {
    states_wide <- tibble::tibble(Gene = gene_universe)
  }

  missing_cols <- setdiff(expected_colnames, setdiff(names(states_wide), "Gene"))
  for (cl in missing_cols) {
    states_wide[[cl]] <- NEUTRAL_HMM_STATE
  }

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

  cn <- setdiff(names(states_wide), "Gene")
  if (!length(cn)) {
    warning("No state columns for ", sample_id)
    return(invisible(NULL))
  }
  cn_ord <- order_state_colnames(sample_id, cn)
  states_wide <- states_wide[, c("Gene", cn_ord), drop = FALSE]

  out_dir <- file.path(OUT_ROOT, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_tsv(
    states_wide,
    file.path(out_dir, "InferCNV_states_all_groups.tsv"),
    na = ""
  )

  if (!"cell" %in% names(cg) || !"cell_group_name" %in% names(cg)) {
    stop("Unexpected cell_groupings columns in ", cg_path)
  }

  labels <- tibble::tibble(
    barcode = cg$cell,
    infercnv_cell_group = cg$cell_group_name,
    InferCNV_state = vapply(
      cg$cell_group_name,
      function(g) infercnv_state_colname_all_groups(sample_id, g),
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

  readr::write_tsv(
    labels_out,
    file.path(out_dir, "InferCNV_labels_all_groups.tsv"),
    na = ""
  )

  lab_with_state <- labels %>%
    dplyr::filter(!is.na(.data$InferCNV_state), nzchar(.data$InferCNV_state))

  if (nrow(lab_with_state)) {
    tab <- table(
      seurat_cluster = as.character(lab_with_state$seurat_clusters),
      InferCNV_state = lab_with_state$InferCNV_state
    )
    ct <- as.data.frame.matrix(tab)
    ct <- tibble::rownames_to_column(ct, "seurat_clusters")
    readr::write_tsv(
      ct,
      file.path(out_dir, "sn_CNV_cluster_count_all_groups.tsv"),
      na = "0"
    )
  } else {
    message("No labeled cells for contingency: ", sample_id)
  }

  message("Wrote all-groups state summary for ", sample_id, " -> ", out_dir)
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
message("N samples: ", length(SAMPLES))
if (length(SAMPLES) <= 20L) {
  message("Samples: ", paste(SAMPLES, collapse = ", "))
} else {
  message(
    "First/last: ", SAMPLES[1L], " ... ", SAMPLES[length(SAMPLES)]
  )
}

for (sid in SAMPLES) {
  process_one_sample_all_groups(sid, so)
}

message("Done.")
