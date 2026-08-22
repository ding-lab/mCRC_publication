#!/usr/bin/env Rscript
# Genome-scale InferCNV HMM heatmaps (ComplexHeatmap), adapted from:
#   MetNet .../InferCNV/src/archive/inferCNV_states_summary_ver1.1.R
# Inputs per sample (first pair found under INFERCNV_STATE_SUMMARY_DIR/<sample>/):
#   Tumor-only: InferCNV_states.tsv + sn_CNV_cluster_count.tsv
#   All groups: InferCNV_states_all_groups.tsv + sn_CNV_cluster_count_all_groups.tsv
#   Override basenames: INFERCNV_STATES_BASENAME, INFERCNV_CLUSTER_COUNT_BASENAME (both required).
#
# HEATMAP_GROUPS=tumor (default): only columns matching <sample_id>_Tumor_s* (CMETS bar + labels).
# HEATMAP_GROUPS=all: all InferCNV state columns; row labels parsed from <sample_id>_<suffix>.
#
# Bottom gene band: genes from TCGA_CNV_candidate_genes_v3.txt only (same file as cytoband track).
#
# Does not require metnet helper.R or googlesheets: embeds average_in_window
# from https://jokergoo.github.io/2020/10/29/make-genome-scale-heatmap
#
# Usage:
#   Rscript plot_infercnv_state_summary_heatmap.R
# Optional env:
#   INFERCNV_STATE_SUMMARY_DIR (default: 20260403 all-groups state summary on MetNet disk)
#   HEATMAP_SAMPLES — comma-separated library IDs (default: six multi-subclone CMETS samples).
#     Set to ALL to plot every sample under INFERCNV_STATE_SUMMARY_DIR.
#   OUT_DIR  MIN_CLONE_FRAC (default 0.01 = 1% of sample cells per InferCNV column)  WINDOW_BP
#   HEATMAP_GROUPS (tumor|all)
#   HEATMAP_CLUSTER_ROWS — false (default): rows blocked by sample, then Tumor_s0, s1, s2, ...
#     true: hierarchical clustering by CNV profile (same sample can split apart).
#   OUT_PDF_FILE — full path for PDF; if unset, name includes basename(INFERCNV_STATE_SUMMARY_DIR)
#   GENCODE_GENE_POS  TCGA_GENES_TSV (v3 table: Cytoband, chr, start, end, Amp_Del, Candidate_gene)
#   CMETS_FISHER_CSV — use cmets_infercnv_subclone_enrichment_fisher.csv (all subclones).
#     Not cmets_infercnv_fisher_sig_with_kras_gain.csv (Fisher-significant subset only — will miss rows).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(circlize)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(colorspace)
  library(grid)
})

average_in_window <- function(window, gr, v, method = "weighted", empty_v = NA) {
  opw <- options(warn = -1)
  on.exit(options(opw), add = TRUE)
  if (missing(v)) v <- rep(1, length(gr))
  if (is.null(v)) v <- rep(1, length(gr))
  if (is.atomic(v) && is.vector(v)) v <- cbind(v)
  v <- as.matrix(v)
  if (is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  if (length(empty_v) == 1) {
    empty_v <- rep(empty_v, ncol(v))
  }
  u <- matrix(rep(empty_v, each = length(window)),
    nrow = length(window), ncol = ncol(v)
  )
  mtch <- as.matrix(findOverlaps(window, gr))
  intersect <- pintersect(window[mtch[, 1]], gr[mtch[, 2]])
  w <- width(intersect)
  v <- v[mtch[, 2], , drop = FALSE]
  n <- nrow(v)
  ind_list <- split(seq_len(n), mtch[, 1])
  window_index <- as.numeric(names(ind_list))
  window_w <- width(window)
  if (is.character(v)) {
    for (i in seq_along(ind_list)) {
      ind <- ind_list[[i]]
      if (is.function(method)) {
        u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
      } else {
        tb <- tapply(w[ind], v[ind], sum)
        u[window_index[i], ] <- names(tb[which.max(tb)])
      }
    }
  } else {
    if (method == "w0") {
      gr2 <- reduce(gr, min.gapwidth = 0)
      mtch2 <- as.matrix(findOverlaps(window, gr2))
      intersect2 <- pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      width_intersect <- tapply(width(intersect2), mtch2[, 1], sum)
      ind <- unique(mtch2[, 1])
      width_setdiff <- width(window[ind]) - width_intersect
      w2 <- width(window[ind])
      for (i in seq_along(ind_list)) {
        ind <- ind_list[[i]]
        x <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
        u[window_index[i], ] <- (x * width_intersect[i] + empty_v * width_setdiff[i]) / w2[i]
      }
    } else if (method == "absolute") {
      for (i in seq_along(ind_list)) {
        u[window_index[i], ] <- colMeans(v[ind_list[[i]], , drop = FALSE])
      }
    } else if (method == "weighted") {
      for (i in seq_along(ind_list)) {
        ind <- ind_list[[i]]
        u[window_index[i], ] <- colSums(v[ind, , drop = FALSE] * w[ind]) / sum(w[ind])
      }
    } else {
      if (is.function(method)) {
        for (i in seq_along(ind_list)) {
          ind <- ind_list[[i]]
          u[window_index[i], ] <- method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  u
}

## Same idea as EnrichedHeatmap::makeWindows (no EnrichedHeatmap dependency)
make_genome_windows <- function(chr_gr, w) {
  lst <- lapply(seq_len(length(chr_gr)), function(i) {
    st <- start(chr_gr[i])
    en <- end(chr_gr[i])
    ch <- as.character(seqnames(chr_gr)[i])
    breaks <- seq.int(st, en, by = w)
    GRanges(
      seqnames = ch,
      ranges = IRanges(start = breaks, end = pmin(breaks + as.integer(w) - 1L, en))
    )
  })
  do.call(c, lst)
}

state_summary_dir <- Sys.getenv(
  "INFERCNV_STATE_SUMMARY_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20260403_inferCNV_state_summary_all_groups"
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/CMETS/InferCNV"
)
min_clone_frac <- as.numeric(Sys.getenv("MIN_CLONE_FRAC", unset = "0.01"))
window_bp <- as.numeric(Sys.getenv("WINDOW_BP", unset = "5e6"))
heatmap_groups <- tolower(trimws(Sys.getenv("HEATMAP_GROUPS", unset = "tumor")))
if (!heatmap_groups %in% c("tumor", "all")) {
  stop("HEATMAP_GROUPS must be 'tumor' or 'all', got: ", heatmap_groups, call. = FALSE)
}

heatmap_cluster_rows_env <- tolower(trimws(Sys.getenv("HEATMAP_CLUSTER_ROWS", unset = "false")))
heatmap_cluster_rows <- heatmap_cluster_rows_env %in% c("1", "true", "yes", "on")

states_basename_env <- trimws(Sys.getenv("INFERCNV_STATES_BASENAME", unset = ""))
ct_basename_env <- trimws(Sys.getenv("INFERCNV_CLUSTER_COUNT_BASENAME", unset = ""))
if (nzchar(states_basename_env) != nzchar(ct_basename_env)) {
  stop(
    "Set both INFERCNV_STATES_BASENAME and INFERCNV_CLUSTER_COUNT_BASENAME, or neither.",
    call. = FALSE
  )
}
use_explicit_basenames <- nzchar(states_basename_env) && nzchar(ct_basename_env)

#' Resolve states + contingency paths for one sample (tumor-only vs all-groups).
resolve_state_summary_inputs <- function(root, sample_id) {
  d <- file.path(root, sample_id)
  if (use_explicit_basenames) {
    st <- file.path(d, states_basename_env)
    ct <- file.path(d, ct_basename_env)
    if (file.exists(st) && file.exists(ct)) {
      return(list(states_fp = st, ct_fp = ct))
    }
    return(NULL)
  }
  st1 <- file.path(d, "InferCNV_states.tsv")
  ct1 <- file.path(d, "sn_CNV_cluster_count.tsv")
  if (file.exists(st1) && file.exists(ct1)) {
    return(list(states_fp = st1, ct_fp = ct1))
  }
  st2 <- file.path(d, "InferCNV_states_all_groups.tsv")
  ct2 <- file.path(d, "sn_CNV_cluster_count_all_groups.tsv")
  if (file.exists(st2) && file.exists(ct2)) {
    return(list(states_fp = st2, ct_fp = ct2))
  }
  NULL
}

gencode_path <- Sys.getenv(
  "GENCODE_GENE_POS",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/src/gencode_v32_gene_name_with_band.txt"
)
tcga_path <- Sys.getenv(
  "TCGA_GENES_TSV",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/src/TCGA_CNV_candidate_genes_v3.txt"
)
cmets_fisher_path <- Sys.getenv(
  "CMETS_FISHER_CSV",
  unset = file.path(
    "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/CMETS/InferCNV",
    "cmets_infercnv_subclone_enrichment_fisher.csv"
  )
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- list.dirs(state_summary_dir, full.names = FALSE, recursive = FALSE)
sample_dirs <- sample_dirs[sample_dirs != ""]
sample_dirs <- sample_dirs[vapply(
  sample_dirs,
  function(s) !is.null(resolve_state_summary_inputs(state_summary_dir, s)),
  logical(1L)
)]
if (!length(sample_dirs)) {
  stop(
    "No sample subfolders with a valid states + cluster-count pair under ",
    state_summary_dir,
    call. = FALSE
  )
}

## Default: multi-subclone CMETS panel (>1% rule via MIN_CLONE_FRAC); override or use ALL for full cohort.
default_heatmap_samples <- paste(
  c(
    "CM349C1-S1",
    "CM426C1-Th1",
    "HT225C1-Th1",
    "HT230C1-Th1",
    "HT260C1-Th1"
  ),
  collapse = ","
)
heatmap_samples_raw <- trimws(Sys.getenv("HEATMAP_SAMPLES", unset = default_heatmap_samples))
if (nzchar(heatmap_samples_raw) && toupper(heatmap_samples_raw) != "ALL") {
  allow <- trimws(strsplit(heatmap_samples_raw, ",", fixed = TRUE)[[1L]])
  allow <- allow[nzchar(allow)]
  miss <- setdiff(allow, sample_dirs)
  if (length(miss)) {
    warning(
      "HEATMAP_SAMPLES not found under state summary dir (skipped): ",
      paste(miss, collapse = ", ")
    )
  }
  sample_dirs <- allow[allow %in% sample_dirs]
}
if (!length(sample_dirs)) {
  stop(
    "No samples left after HEATMAP_SAMPLES filter. Set HEATMAP_SAMPLES=all to use all libraries.",
    call. = FALSE
  )
}

message("State summary dir: ", state_summary_dir)
message("HEATMAP_GROUPS: ", heatmap_groups)
message("HEATMAP_CLUSTER_ROWS: ", heatmap_cluster_rows)
message("MIN_CLONE_FRAC: ", min_clone_frac)
message("Samples: ", paste(sample_dirs, collapse = ", "))

summary_states_tbl <- data.frame()

for (sample_id in sample_dirs) {
  inp <- resolve_state_summary_inputs(state_summary_dir, sample_id)
  if (is.null(inp)) {
    warning("Skip ", sample_id, " (missing states or contingency)")
    next
  }
  states_fp <- inp$states_fp
  ct_fp <- inp$ct_fp
  states_tbl <- read_tsv(states_fp, show_col_types = FALSE) %>%
    column_to_rownames("Gene") %>%
    as.data.frame()
  contingency_tbl <- read_tsv(ct_fp, show_col_types = FALSE) %>%
    column_to_rownames("seurat_clusters")
  contingency_tbl[] <- lapply(contingency_tbl, function(x) as.numeric(as.character(x)))
  total_cells <- sum(as.matrix(contingency_tbl), na.rm = TRUE)
  cs <- colSums(as.matrix(contingency_tbl), na.rm = TRUE)
  filtered_subgroups <- names(cs)[cs >= (min_clone_frac * total_cells)]
  if (!length(filtered_subgroups)) {
    warning(sample_id, ": no subclone columns pass MIN_CLONE_FRAC; using all.")
    filtered_subgroups <- colnames(states_tbl)
  }
  keep <- intersect(filtered_subgroups, colnames(states_tbl))
  if (heatmap_groups == "tumor") {
    pfx <- paste0(sample_id, "_")
    npx <- nchar(pfx)
    is_tumor_col <- vapply(keep, function(cn) {
      startsWith(cn, pfx) &&
        grepl("^Tumor_s[0-9]+$", substring(cn, npx + 1L), perl = TRUE)
    }, logical(1L))
    keep <- keep[is_tumor_col]
  }
  if (!length(keep)) {
    warning(sample_id, ": no columns after HEATMAP_GROUPS/min_frac filter; skip.")
    next
  }
  filtered_states_tbl <- states_tbl[, keep, drop = FALSE]
  filtered_states_tbl <- t(filtered_states_tbl) %>% as.data.frame()
  summary_states_tbl <- bind_rows(summary_states_tbl, filtered_states_tbl)
}

if (!nrow(summary_states_tbl)) {
  stop("No data combined; check inputs.")
}

summary_states_tbl <- summary_states_tbl %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Gene")

infercnv_hmm <- summary_states_tbl %>%
  column_to_rownames("Gene") %>%
  mutate_all(function(x) {
    dplyr::case_when(
      is.na(x) ~ 3,
      TRUE ~ as.numeric(x)
    )
  }) %>%
  rownames_to_column("Gene")

gene_pos <- read_tsv(gencode_path, col_names = FALSE, show_col_types = FALSE)
if (ncol(gene_pos) >= 4L) {
  colnames(gene_pos)[1:4] <- c("Gene", "Chr", "Start", "End")
} else {
  stop("GENCODE file needs >=4 columns: Gene Chr Start End[, ...]")
}
gene_pos <- gene_pos[, 1:4]

std_chr <- paste0("chr", 1:22)
keep_std_chr <- function(gr) {
  gr[as.character(GenomicRanges::seqnames(gr)) %in% std_chr]
}

infercnv_hmm_bed <- infercnv_hmm %>%
  left_join(gene_pos, by = "Gene") %>%
  relocate(Gene, Chr, Start, End) %>%
  filter(!is.na(Chr), !is.na(Start), !is.na(End)) %>%
  filter(as.character(Chr) %in% std_chr) %>%
  column_to_rownames("Gene") %>%
  mutate_at(vars(-Chr, -Start, -End), as.numeric)

infercnv_hmm_gr <- GRanges(
  seqnames = infercnv_hmm_bed[, 1],
  ranges = IRanges(infercnv_hmm_bed[, 2], infercnv_hmm_bed[, 3])
)

chr_df <- read.chromInfo()$df
chr_df <- chr_df[chr_df$chr %in% paste0("chr", 1:22), ]
chr_gr <- GRanges(
  seqnames = chr_df[, 1],
  ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3])
)
chr_window <- make_genome_windows(chr_gr, w = as.integer(window_bp))

num_mat <- average_in_window(
  chr_window, infercnv_hmm_gr,
  infercnv_hmm_bed[, -(1:3), drop = FALSE],
  method = "absolute"
)
colnames(num_mat) <- colnames(infercnv_hmm_bed[, -(1:3), drop = FALSE])

## Row meta: matrix columns = inferCNV state columns, e.g. CM329C1-S1_Tumor_s2
infercnv_colnames <- colnames(num_mat)
sample_ids_by_len <- sample_dirs[order(-nchar(sample_dirs), sample_dirs)]
parse_infercnv_matrix_colname <- function(col, sids_ordered) {
  for (sid in sids_ordered) {
    pfx <- paste0(sid, "_")
    if (startsWith(col, pfx)) {
      suf <- substring(col, nchar(pfx) + 1L)
      return(list(sample = sid, suffix = suf))
    }
  }
  list(sample = NA_character_, suffix = col)
}
parsed_cols <- lapply(infercnv_colnames, parse_infercnv_matrix_colname, sample_ids_by_len)
sample_per_infercnv_row <- vapply(parsed_cols, function(z) z$sample, character(1L))
clone_per_infercnv_row <- vapply(parsed_cols, function(z) z$suffix, character(1L))
if (any(is.na(sample_per_infercnv_row))) {
  warning(
    "Some heatmap columns did not match any sample_id prefix (labels may be wrong): ",
    paste(infercnv_colnames[is.na(sample_per_infercnv_row)], collapse = ", ")
  )
}

## Normalize Tumor_s* for joining (Tumor_s01 -> Tumor_s1)
norm_tumor_subclone <- function(x) {
  x <- trimws(as.character(x))
  i <- grepl("^Tumor_s[0-9]+$", x)
  if (any(i)) {
    n <- as.integer(sub("^Tumor_s", "", x[i]))
    x[i] <- paste0("Tumor_s", n)
  }
  x
}

prop_cmet_per_row <- rep(NA_real_, length(infercnv_colnames))
if (file.exists(cmets_fisher_path)) {
  fisher_tbl <- read_csv(cmets_fisher_path, show_col_types = FALSE)
  reqf <- c("orig.ident", "infercnv_tumor_clone", "prop_noncanonical_in_clone")
  if (!all(reqf %in% names(fisher_tbl))) {
    warning(
      "CMETS_FISHER_CSV missing columns ",
      paste(setdiff(reqf, names(fisher_tbl)), collapse = ", "),
      "; skipping barplot values."
    )
  } else {
    fisher_tbl <- fisher_tbl %>%
      mutate(
        .sid = trimws(as.character(.data$orig.ident)),
        .clone = norm_tumor_subclone(.data$infercnv_tumor_clone),
        .key = paste0(.data$.sid, "\x1f", .data$.clone),
        .prop = as.numeric(.data$prop_noncanonical_in_clone)
      )
    rk <- paste0(
      trimws(sample_per_infercnv_row),
      "\x1f",
      norm_tumor_subclone(clone_per_infercnv_row)
    )
    mi <- match(rk, fisher_tbl$.key)
    prop_cmet_per_row <- fisher_tbl$.prop[mi]
    n_miss <- sum(is.na(mi))
    if (n_miss > 0L) {
      warning(
        n_miss,
        " inferCNV row(s) not found in CMETS Fisher table; bar uses 0. ",
        "Use cmets_infercnv_subclone_enrichment_fisher.csv (all subclones), not *_fisher_sig_with_kras_gain.csv. ",
        "Missing: ",
        paste(infercnv_colnames[is.na(mi)], collapse = ", ")
      )
      prop_cmet_per_row[is.na(mi)] <- 0
    }
    prop_cmet_per_row[is.na(prop_cmet_per_row)] <- 0
  }
} else {
  warning("CMETS Fisher CSV not found (barplot zeros): ", cmets_fisher_path)
}

chr <- as.vector(seqnames(chr_window))
chr_level <- paste0("chr", 1:22)
chr <- factor(chr, levels = chr_level)

## TCGA v3 table: cytoband track + sole gene-band marks (Candidate_gene; skip NA names)
if (!file.exists(tcga_path)) {
  stop("TCGA genes file not found: ", tcga_path, call. = FALSE)
}
tcga_genes <- read_tsv(tcga_path, show_col_types = FALSE)
if (!"Cytoband" %in% names(tcga_genes)) {
  stop("TCGA file must contain column Cytoband: ", tcga_path, call. = FALSE)
}
cn <- names(tcga_genes)
cn[cn == "chr"] <- "Chr"
cn[cn == "start"] <- "Start"
cn[cn == "end"] <- "End"
names(tcga_genes) <- cn
if (!all(c("Chr", "Start", "End") %in% names(tcga_genes))) {
  stop("TCGA file needs chr/start/end (any case): ", tcga_path, call. = FALSE)
}

tcga_genes_bed <- tcga_genes %>%
  column_to_rownames("Cytoband") %>%
  select(-dplyr::any_of("Amp_Del")) %>%
  filter(!is.na(.data$Candidate_gene), nzchar(as.character(.data$Candidate_gene)))

## Gene band: one position per gene (mid-cytoband) for window overlap
mid_gene <- as.integer((tcga_genes_bed$Start + tcga_genes_bed$End) / 2)
tcga_genes_gr <- GRanges(
  seqnames = tcga_genes_bed$Chr,
  ranges = IRanges(mid_gene, mid_gene)
)
tcga_genes_gr$gene <- tcga_genes_bed$Candidate_gene
tcga_genes_gr <- keep_std_chr(tcga_genes_gr)

mtch <- as.matrix(findOverlaps(chr_window, tcga_genes_gr))
at_gene_band <- mtch[, 1]
labels_gene_band <- as.character(tcga_genes_gr$gene[mtch[, 2]])

## TCGA cytoband amp/deletion track (numeric)
tcga_sig_bed <- tcga_genes %>%
  column_to_rownames("Cytoband") %>%
  select(-dplyr::any_of("Candidate_gene"))
tcga_sig_gr <- GRanges(
  seqnames = tcga_sig_bed$Chr,
  ranges = IRanges(as.integer(tcga_sig_bed$Start), as.integer(tcga_sig_bed$End))
)
tcga_sig_gr <- keep_std_chr(tcga_sig_gr)
tcga_sig_num_mat <- average_in_window(
  chr_window, tcga_sig_gr,
  tcga_sig_bed[, !(names(tcga_sig_bed) %in% c("Chr", "Start", "End")), drop = FALSE],
  method = "absolute"
)

ht_opt$TITLE_PADDING <- unit(c(4, 4), "points")

## HMM discrete states 1–6: 3 = neutral; blue = loss, red = amplification (averaged in windows).
col_fun <- colorRamp2(
  c(1, 2, 3, 4, 5, 6),
  c("#053061", "#4393C3", "#F7F7F7", "#F4A582", "#D6604D", "#67001D")
)

hm_infercnv <- t(num_mat)
rownames(hm_infercnv) <- infercnv_colnames
names(prop_cmet_per_row) <- infercnv_colnames
prop_cmet_per_row <- unname(prop_cmet_per_row[infercnv_colnames])
if (length(prop_cmet_per_row) != nrow(hm_infercnv)) {
  stop("CMETS bar length mismatch vs heatmap rows (internal error).")
}

## Left labels: sample + subclone (unique per row; duplicate sample-only labels confuse CH + annotations)
row_label_infercnv <- paste0(
  sample_per_infercnv_row,
  "\n",
  clone_per_infercnv_row
)

## Optional: keep rows blocked by sample (HT225 together) instead of clustering by profile similarity.
if (!heatmap_cluster_rows) {
  sample_rank <- match(sample_per_infercnv_row, sample_dirs)
  sample_rank[is.na(sample_rank)] <- 9999L
  subn <- rep(9999L, length(clone_per_infercnv_row))
  is_ts <- grepl("^Tumor_s[0-9]+$", clone_per_infercnv_row, perl = TRUE)
  subn[is_ts] <- suppressWarnings(as.integer(sub("^Tumor_s", "", clone_per_infercnv_row[is_ts])))
  subn[is.na(subn)] <- 9999L
  ord <- order(sample_rank, subn, infercnv_colnames)
  hm_infercnv <- hm_infercnv[ord, , drop = FALSE]
  infercnv_colnames <- infercnv_colnames[ord]
  sample_per_infercnv_row <- sample_per_infercnv_row[ord]
  clone_per_infercnv_row <- clone_per_infercnv_row[ord]
  prop_cmet_per_row <- prop_cmet_per_row[ord]
  row_label_infercnv <- row_label_infercnv[ord]
  rownames(hm_infercnv) <- infercnv_colnames
}

## Left strip: same color for all rows from one sample (order matches hm_infercnv rows)
sample_ids_unique <- unique(sample_per_infercnv_row)
n_samples_anno <- length(sample_ids_unique)
sample_palette <- colorspace::qualitative_hcl(
  pmax(3L, n_samples_anno),
  palette = "Dark 3"
)[seq_len(n_samples_anno)]
names(sample_palette) <- sample_ids_unique

left_sample <- rowAnnotation(
  Sample = sample_per_infercnv_row,
  col = list(Sample = sample_palette),
  show_annotation_name = FALSE,
  annotation_legend_param = list(title = "Sample", title_gp = gpar(fontsize = 8)),
  border = TRUE,
  width = unit(5, "mm")
)

## Bar height = prop_noncanonical_in_clone from CMETS_FISHER_CSV (axis_param: no labels_gp/gp — older CH rejects them)
## ylim c(0, 0.6): proportion axis; values > 0.6 clip at top of bar. use_raster=FALSE: body vs rowAnnotation align.
right_bar <- rowAnnotation(
  CMETS_prop = anno_barplot(
    prop_cmet_per_row,
    gp = gpar(fill = "#4575b4", col = NA),
    border = TRUE,
    axis = TRUE,
    ylim = c(0, 0.6),
    axis_param = list(side = "top", at = c(0, 0.2, 0.4, 0.6)),
    width = unit(2.8, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 8),
  show_annotation_name = TRUE
)

ht_list <- Heatmap(
  hm_infercnv,
  name = "HMM",
  col = col_fun,
  column_split = chr,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  cluster_rows = heatmap_cluster_rows,
  row_order = if (heatmap_cluster_rows) NULL else seq_len(nrow(hm_infercnv)),
  row_title = "inferCNV",
  row_labels = row_label_infercnv,
  row_names_side = "left",
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(12, "cm"),
  left_annotation = left_sample,
  right_annotation = right_bar,
  column_title_gp = gpar(fontsize = 8),
  border = TRUE,
  column_gap = unit(0, "points"),
  column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
  heatmap_legend_param = list(direction = "horizontal", title_position = "lefttop"),
  use_raster = FALSE
)

ht_cytoband <- Heatmap(
  t(tcga_sig_num_mat),
  name = "tcga_band",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  column_split = chr,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  row_title = "",
  column_title_gp = gpar(fontsize = 8),
  border = TRUE,
  column_gap = unit(0, "points"),
  column_title = ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n")),
  heatmap_legend_param = list(show_legend = FALSE, direction = "horizontal"),
  height = unit(10, "points")
)

## Gene band: v3 candidate genes only, ordered by genome position
mark_at <- at_gene_band
mark_lab <- labels_gene_band
if (length(mark_at)) {
  ord <- order(mark_at, mark_lab)
  mark_at <- mark_at[ord]
  mark_lab <- mark_lab[ord]
}

message("TCGA gene table: ", tcga_path)
message("CMETS Fisher table: ", cmets_fisher_path)

out_pdf_env <- trimws(Sys.getenv("OUT_PDF_FILE", unset = ""))
if (nzchar(out_pdf_env)) {
  pdf_file <- out_pdf_env
} else {
  state_tag <- basename(state_summary_dir)
  state_tag <- gsub("[^A-Za-z0-9_.-]+", "_", state_tag)
  pdf_file <- file.path(
    out_dir,
    paste0("inferCNV_states_", state_tag, "_rowCluster_TCGA_v3_gene_band_cmetsBar.pdf")
  )
}
dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
pdf(pdf_file, width = 20, height = 8)
if (length(mark_at)) {
  anno_genes <- HeatmapAnnotation(
    gene_band = anno_mark(
      at = mark_at,
      labels = mark_lab,
      side = "bottom",
      labels_gp = gpar(col = "black", fontsize = 7),
      link_gp = gpar(col = "black", lwd = 1)
    ),
    show_annotation_name = FALSE
  )
  draw(
    ht_list %v% ht_cytoband %v% anno_genes,
    heatmap_legend_side = "right",
    merge_legend = TRUE
  )
} else {
  warning("No v3 candidate gene overlaps with genome windows; drawing without gene band.")
  draw(
    ht_list %v% ht_cytoband,
    heatmap_legend_side = "right",
    merge_legend = TRUE
  )
}

dev.off()

message("Wrote ", pdf_file)
