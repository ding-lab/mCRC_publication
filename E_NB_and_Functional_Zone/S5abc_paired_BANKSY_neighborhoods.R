#!/usr/bin/env Rscript

# Extended Data Fig. 5a-c: paired primary-metastasis BANKSY neighborhoods.
# Selected Leiden resolutions: CM819C/CM579C/HT472C 0.7, CM626C/CM663C/CM798C 0.5,
# CM397C 0.4. BANKSY outputs are read from the pairwise run directory.
# Source: Revesion/Banksy/Run_all_paired_cases.ipynb

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(ggplot2)
library(forcats)

# ===== cell 1 =====
metadata_path <- Sys.getenv(
  "MCRC_XENIUM_METADATA",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv"
)
meta_df <- read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"
head(meta_df, 3)

# ===== cell 2 =====
meta_df <- meta_df %>%
  mutate(
    Case_ID = case_when(
      Case_ID == "CM394C" ~ "CM397C",
      TRUE ~ Case_ID
    )
  )

# ===== cell 3 =====
## ============================================================
## Paths
## ============================================================

banksy_dir <- Sys.getenv(
  "MCRC_BANKSY_PAIRED_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Banksy/Individiual_paired/out_pairwise/2-banksy_multi_leiden"
)
out_root <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = file.path(
    "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/E_NB_and_Functional_Zone",
    "paired_BANKSY_output"
  )
)
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## Case-specific selected resolutions
## ============================================================

case_resolution_tbl <- tribble(
  ~case_id, ~selected_res_label, ~selected_banksy_col,
  "CM819C", "res_0p7", "banksy_pair_leiden_0p7",
  "CM579C", "res_0p7", "banksy_pair_leiden_0p7",
  "CM626C", "res_0p5", "banksy_pair_leiden_0p5",
  "CM663C", "res_0p5", "banksy_pair_leiden_0p5",
  "CM798C", "res_0p5", "banksy_pair_leiden_0p5",
  "CM397C", "res_0p4", "banksy_pair_leiden_0p4",
  "HT472C", "res_0p7", "banksy_pair_leiden_0p7"
)

## ============================================================
## Cell type order
## Use Broad_cell_type1
## Doublet will be removed
## ============================================================

preferred_cell_type_order <- c(
  "Intestinal_epithelium",
  "Tumor",
  "Fibroblast",
  "Endothelial_cell",
  "Pericyte",
  "Smooth_muscle",
  "Myeloid cells",
  "T_NK_cell",
  "B_cell",
  "Plasma_cell",
  "Hepatocyte",
  "Cholangiocyte",
  "Breast_duct",
  "Aveolar_epithelium",
  "Necrosis"
)

## ============================================================
## Parameters
## ============================================================

min_cluster_pct_total <- 0.1
min_cluster_pct_PSI <- 0.1

min_max_pct_for_top_change <- 2
min_abs_delta_pct_for_top_change <- 2
top_n_changed <- 3

# ===== cell 4 =====
read_case_banksy_all_obs <- function(case_id) {
  
  all_obs_dir <- file.path(
    banksy_dir,
    case_id,
    "leiden_clusters",
    "res_0p7"
  )
  
  all_obs_files <- list.files(
    all_obs_dir,
    pattern = "_all_obs\\.csv$",
    full.names = TRUE
  )
  
  if (length(all_obs_files) == 0) {
    stop("No *_all_obs.csv found for case: ", case_id)
  }
  
  preferred_file <- all_obs_files[
    basename(all_obs_files) == paste0(case_id, "_leiden_res_0p7_all_obs.csv")
  ]
  
  if (length(preferred_file) == 1) {
    all_obs_file <- preferred_file
  } else {
    all_obs_file <- all_obs_files[1]
    warning("Using first *_all_obs.csv for ", case_id, ": ", all_obs_file)
  }
  
  message("Reading BANKSY: ", all_obs_file)
  
  banksy_raw <- read_csv(all_obs_file, show_col_types = FALSE)
  
  if (!all(c("dataset", "cell_id") %in% colnames(banksy_raw))) {
    stop("BANKSY file for ", case_id, " does not contain dataset and cell_id columns.")
  }
  
  if ("group" %in% colnames(banksy_raw) && !"leiden_res_0p7" %in% colnames(banksy_raw)) {
    banksy_raw <- banksy_raw %>%
      rename(leiden_res_0p7 = group)
  }
  
  banksy_tbl <- banksy_raw %>%
    mutate(
      barcode = paste0(dataset, "__", cell_id)
    ) %>%
    rename_with(
      ~ str_replace(.x, "^leiden_res_", "banksy_pair_leiden_"),
      starts_with("leiden_res_")
    ) %>%
    select(
      barcode,
      dataset,
      cell_id,
      starts_with("banksy_pair_leiden_")
    )
  
  return(banksy_tbl)
}


merge_case_metadata_banksy <- function(case_id) {
  
  case_meta <- meta_df %>%
    filter(Case_ID == case_id)
  
  banksy_tbl <- read_case_banksy_all_obs(case_id)
  
  case_meta_banksy <- case_meta %>%
    left_join(
      banksy_tbl %>%
        select(barcode, starts_with("banksy_pair_leiden_")),
      by = "barcode"
    )
  
  n_missing <- case_meta_banksy %>%
    summarise(n_missing = sum(is.na(banksy_pair_leiden_0p7))) %>%
    pull(n_missing)
  
  message(case_id, ": missing BANKSY assignment = ", n_missing)
  
  return(case_meta_banksy)
}

# ===== cell 5 =====
make_case_celltype_enrichment_heatmap <- function(case_meta_banksy,
                                                  case_id,
                                                  selected_banksy_col,
                                                  selected_res_label) {
  
  out_dir <- file.path(out_root, case_id, "celltype_enrichment_heatmap")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  meta_df2 <- case_meta_banksy %>%
    filter(!is.na(.data[[selected_banksy_col]])) %>%
    filter(!is.na(Broad_cell_type1), Broad_cell_type1 != "") %>%
    filter(Broad_cell_type1 != "Doublet") %>%
    mutate(
      banksy_cluster = as.character(.data[[selected_banksy_col]]),
      Broad_cell_type1 = as.character(Broad_cell_type1)
    )
  
  present_cell_types <- sort(unique(meta_df2$Broad_cell_type1))
  
  cell_types <- c(
    preferred_cell_type_order[preferred_cell_type_order %in% present_cell_types],
    setdiff(present_cell_types, preferred_cell_type_order)
  )
  
  total_non_doublet_cells <- nrow(meta_df2)
  
  cluster_size_tbl <- meta_df2 %>%
    count(banksy_cluster, name = "n_cluster_cells") %>%
    mutate(
      pct_total_non_doublet_cells = n_cluster_cells / total_non_doublet_cells * 100,
      keep_cluster = pct_total_non_doublet_cells >= min_cluster_pct_total
    ) %>%
    arrange(as.numeric(banksy_cluster))
  
  valid_clusters <- cluster_size_tbl %>%
    filter(keep_cluster) %>%
    pull(banksy_cluster)
  
  meta_df2 <- meta_df2 %>%
    filter(banksy_cluster %in% valid_clusters)
  
  neighborhoods <- cluster_size_tbl %>%
    filter(keep_cluster) %>%
    arrange(as.numeric(banksy_cluster)) %>%
    pull(banksy_cluster)
  
  ## ----------------------------
  ## Cluster composition for row ordering
  ## ----------------------------
  
  cluster_celltype_prop_tbl <- meta_df2 %>%
    count(
      banksy_cluster,
      Broad_cell_type1,
      name = "n_cells"
    ) %>%
    group_by(banksy_cluster) %>%
    mutate(
      total_cluster_cells = sum(n_cells),
      pct_within_cluster = n_cells / total_cluster_cells * 100
    ) %>%
    ungroup()
  
  cluster_celltype_wide_tbl <- cluster_celltype_prop_tbl %>%
    select(banksy_cluster, Broad_cell_type1, pct_within_cluster) %>%
    pivot_wider(
      names_from = Broad_cell_type1,
      values_from = pct_within_cluster,
      values_fill = 0
    )
  
  missing_ct_cols <- setdiff(cell_types, colnames(cluster_celltype_wide_tbl))
  for (ct in missing_ct_cols) {
    cluster_celltype_wide_tbl[[ct]] <- 0
  }
  
  cluster_order_tbl <- cluster_celltype_wide_tbl %>%
    rowwise() %>%
    mutate(
      dominant_cell_type = cell_types[which.max(c_across(all_of(cell_types)))],
      dominant_pct = max(c_across(all_of(cell_types))),
      dominant_cell_type_rank = match(dominant_cell_type, cell_types)
    ) %>%
    ungroup() %>%
    arrange(
      dominant_cell_type_rank,
      desc(dominant_pct),
      as.numeric(banksy_cluster)
    )
  
  row_order <- cluster_order_tbl$banksy_cluster
  
  ## ----------------------------
  ## Fisher test matrices
  ## ----------------------------
  
  odds_matrix <- matrix(
    NA,
    nrow = length(neighborhoods),
    ncol = length(cell_types),
    dimnames = list(neighborhoods, cell_types)
  )
  
  pval_matrix <- matrix(
    NA,
    nrow = length(neighborhoods),
    ncol = length(cell_types),
    dimnames = list(neighborhoods, cell_types)
  )
  
  count_a_matrix <- matrix(
    NA,
    nrow = length(neighborhoods),
    ncol = length(cell_types),
    dimnames = list(neighborhoods, cell_types)
  )
  
  for (nb in neighborhoods) {
    for (ct in cell_types) {
      
      a <- sum(meta_df2$banksy_cluster == nb & meta_df2$Broad_cell_type1 == ct)
      b <- sum(meta_df2$banksy_cluster == nb & meta_df2$Broad_cell_type1 != ct)
      c <- sum(meta_df2$banksy_cluster != nb & meta_df2$Broad_cell_type1 == ct)
      d <- sum(meta_df2$banksy_cluster != nb & meta_df2$Broad_cell_type1 != ct)
      
      tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      
      fisher <- fisher.test(tbl, alternative = "greater")
      
      odds_matrix[nb, ct] <- as.numeric(fisher$estimate)
      pval_matrix[nb, ct] <- fisher$p.value
      count_a_matrix[nb, ct] <- a
    }
  }
  
  ## ----------------------------
  ## FDR and log2OR
  ## ----------------------------
  
  padj_matrix <- matrix(
    p.adjust(as.vector(pval_matrix), method = "fdr"),
    nrow = nrow(pval_matrix),
    ncol = ncol(pval_matrix),
    dimnames = dimnames(pval_matrix)
  )
  
  star_matrix <- ifelse(padj_matrix < 0.05, "*", "")
  
  log2_odds_matrix <- log2(odds_matrix)
  
  max_pos_value <- max(log2_odds_matrix[is.finite(log2_odds_matrix)], na.rm = TRUE)
  max_neg_value <- min(log2_odds_matrix[is.finite(log2_odds_matrix)], na.rm = TRUE)
  
  log2_odds_matrix[
    is.infinite(log2_odds_matrix) & log2_odds_matrix > 0
  ] <- max_pos_value
  
  log2_odds_matrix[
    is.infinite(log2_odds_matrix) & log2_odds_matrix < 0
  ] <- max_neg_value
  
  log2_odds_matrix_capped <- log2_odds_matrix
  log2_odds_matrix_capped[log2_odds_matrix_capped > 3] <- 3
  log2_odds_matrix_capped[log2_odds_matrix_capped < -3] <- -3
  
  ## ----------------------------
  ## Dominant cluster annotation
  ## ----------------------------
  
  banksy_cluster_annotation_tbl <- cluster_celltype_prop_tbl %>%
    filter(banksy_cluster %in% neighborhoods) %>%
    group_by(banksy_cluster) %>%
    arrange(desc(pct_within_cluster), .by_group = TRUE) %>%
    summarise(
      dominant_cell_type = first(Broad_cell_type1),
      dominant_pct = first(pct_within_cluster),
      top3_celltype_composition = paste0(
        head(Broad_cell_type1, 3),
        " (",
        round(head(pct_within_cluster, 3), 1),
        "%)",
        collapse = "; "
      ),
      .groups = "drop"
    ) %>%
    left_join(
      cluster_size_tbl %>%
        select(banksy_cluster, n_cluster_cells, pct_total_non_doublet_cells),
      by = "banksy_cluster"
    ) %>%
    arrange(match(banksy_cluster, row_order))
  
  ## ----------------------------
  ## Heatmap
  ## ----------------------------
  
  col_fun <- colorRamp2(
    breaks = c(-3, 0, 3),
    colors = c("blue", "white", "red")
  )
  
  ht <- Heatmap(
    log2_odds_matrix_capped,
    name = "log2OR",
    col = col_fun,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(
        star_matrix[i, j],
        x,
        y,
        gp = gpar(fontsize = 10)
      )
    },
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "bottom",
    row_order = row_order,
    column_order = cell_types,
    column_names_rot = 45,
    heatmap_legend_param = list(
      title = "log2OR",
      at = c(-3, 0, 3),
      labels = c("≤ -3", "0", "≥ 3")
    )
  )
  
  pdf_file <- file.path(
    out_dir,
    paste0(case_id, "_banksy_", selected_res_label, "_Broad_cell_type1_fisher_enrichment_heatmap.pdf")
  )
  
  png_file <- file.path(
    out_dir,
    paste0(case_id, "_banksy_", selected_res_label, "_Broad_cell_type1_fisher_enrichment_heatmap.png")
  )
  
  pdf(pdf_file, width = max(8, length(cell_types) * 0.45), height = max(5, length(neighborhoods) * 0.28 + 2))
  draw(ht)
  dev.off()
  
  png(png_file, width = 2600, height = 2000, res = 300)
  draw(ht)
  dev.off()
  
  # NOTE: as.vector(matrix) is column-major; expand_grid is row-major
  # (Broad_cell_type1 varies fastest). Use as.vector(t(mat)) to align.
  fisher_result_tbl <- expand_grid(
    banksy_cluster = neighborhoods,
    Broad_cell_type1 = cell_types
  ) %>%
    mutate(
      n_cells = as.vector(t(count_a_matrix)),
      odds_ratio = as.vector(t(odds_matrix)),
      log2_odds_ratio = as.vector(t(log2_odds_matrix)),
      p_value = as.vector(t(pval_matrix)),
      p_adj = as.vector(t(padj_matrix)),
      significant = p_adj < 0.05
    ) %>%
    left_join(
      cluster_size_tbl,
      by = "banksy_cluster"
    ) %>%
    arrange(
      match(banksy_cluster, row_order),
      match(Broad_cell_type1, cell_types)
    )
  
  write_csv(
    fisher_result_tbl,
    file.path(out_dir, paste0(case_id, "_banksy_", selected_res_label, "_Broad_cell_type1_fisher_enrichment_table.csv"))
  )
  
  write_csv(
    banksy_cluster_annotation_tbl,
    file.path(out_dir, paste0(case_id, "_banksy_", selected_res_label, "_cluster_celltype_annotation.csv"))
  )
  
  write_csv(
    cluster_size_tbl,
    file.path(out_dir, paste0(case_id, "_banksy_", selected_res_label, "_cluster_size_filter_table.csv"))
  )
  
  write_csv(
    cluster_order_tbl,
    file.path(out_dir, paste0(case_id, "_banksy_", selected_res_label, "_cluster_row_order_by_celltype_composition.csv"))
  )
  
  return(
    list(
      heatmap = ht,
      fisher_result_tbl = fisher_result_tbl,
      cluster_annotation_tbl = banksy_cluster_annotation_tbl,
      cluster_size_tbl = cluster_size_tbl,
      row_order = row_order,
      cell_types = cell_types
    )
  )
}


# ===== cell 6 =====
make_case_PSI_barplot_and_top_changes <- function(case_meta_banksy,
                                                  case_id,
                                                  selected_banksy_col,
                                                  selected_res_label,
                                                  cluster_annotation_tbl = NULL) {
  
  out_dir <- file.path(out_root, case_id, "PSI_banksy_composition")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ----------------------------
  ## Define PSI
  ## PSI = tns_label > 0 and tn_label not > 0
  ## ----------------------------
  
  case_psi_tbl <- case_meta_banksy %>%
    mutate(
      tn_label_num = suppressWarnings(as.numeric(tn_label)),
      tns_label_num = suppressWarnings(as.numeric(tns_label)),
      tn_label_positive = !is.na(tn_label_num) & tn_label_num > 0,
      tns_label_positive = !is.na(tns_label_num) & tns_label_num > 0,
      is_PSI = tns_label_positive & !tn_label_positive,
      banksy_cluster = as.character(.data[[selected_banksy_col]]),
      sample_role = case_when(
        Organ == "Colon" ~ "Primary",
        Organ %in% c("Liver", "Lung") ~ "Metastasis",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(is_PSI) %>%
    filter(!is.na(banksy_cluster)) %>%
    filter(!is.na(sample_role))
  
  if (nrow(case_psi_tbl) == 0) {
    warning("No PSI cells found for ", case_id)
    return(NULL)
  }
  
  ## ----------------------------
  ## Cluster filters
  ## 1) >=0.1% of total cells
  ## 2) >=0.1% of total PSI cells
  ## ----------------------------
  
  total_cells_pair <- case_meta_banksy %>%
    filter(!is.na(.data[[selected_banksy_col]])) %>%
    nrow()
  
  total_PSI_cells_pair <- nrow(case_psi_tbl)
  
  cluster_size_total_tbl <- case_meta_banksy %>%
    filter(!is.na(.data[[selected_banksy_col]])) %>%
    mutate(
      banksy_cluster = as.character(.data[[selected_banksy_col]])
    ) %>%
    count(banksy_cluster, name = "n_total_cells_pair") %>%
    mutate(
      pct_total_cells_pair = n_total_cells_pair / total_cells_pair * 100
    )
  
  cluster_size_PSI_tbl <- case_psi_tbl %>%
    count(banksy_cluster, name = "n_total_PSI_cells_pair") %>%
    mutate(
      pct_total_PSI_cells_pair = n_total_PSI_cells_pair / total_PSI_cells_pair * 100
    )
  
  valid_banksy_clusters <- cluster_size_total_tbl %>%
    left_join(cluster_size_PSI_tbl, by = "banksy_cluster") %>%
    mutate(
      n_total_PSI_cells_pair = replace_na(n_total_PSI_cells_pair, 0L),
      pct_total_PSI_cells_pair = replace_na(pct_total_PSI_cells_pair, 0)
    ) %>%
    filter(
      pct_total_cells_pair >= min_cluster_pct_total,
      pct_total_PSI_cells_pair >= min_cluster_pct_PSI
    )
  
  ## ----------------------------
  ## PSI composition table
  ## denominator = all PSI cells in each sample
  ## ----------------------------
  
  psi_total_tbl <- case_psi_tbl %>%
    count(sample_id, Organ, sample_role, name = "total_all_PSI_cells")
  
  psi_banksy_comp_tbl <- case_psi_tbl %>%
    semi_join(valid_banksy_clusters, by = "banksy_cluster") %>%
    count(
      sample_id,
      Organ,
      sample_role,
      banksy_cluster,
      name = "n_cells"
    ) %>%
    left_join(
      psi_total_tbl,
      by = c("sample_id", "Organ", "sample_role")
    ) %>%
    mutate(
      frac = n_cells / total_all_PSI_cells,
      pct = frac * 100
    ) %>%
    ungroup()
  
  ## ----------------------------
  ## Barplot
  ## ----------------------------
  
  cluster_levels <- valid_banksy_clusters %>%
    arrange(as.numeric(banksy_cluster)) %>%
    pull(banksy_cluster)
  
  p_psi_banksy <- psi_banksy_comp_tbl %>%
    mutate(
      banksy_cluster = factor(banksy_cluster, levels = cluster_levels),
      sample_label = paste0(sample_id, "\n", Organ)
    ) %>%
    ggplot(
      aes(
        x = sample_label,
        y = pct,
        fill = banksy_cluster
      )
    ) +
    geom_col(width = 0.75, color = "black", linewidth = 0.15) +
    labs(
      title = paste0(case_id, " PSI BANKSY composition"),
      subtitle = paste0(
        selected_res_label,
        "; clusters < ", min_cluster_pct_total,
        "% total cells or < ", min_cluster_pct_PSI,
        "% PSI cells removed"
      ),
      x = NULL,
      y = "Percentage of all PSI cells",
      fill = paste0("BANKSY\n", selected_res_label)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_composition_barplot.pdf")),
    p_psi_banksy,
    width = 5.8,
    height = 4.8
  )
  
  ggsave(
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_composition_barplot.png")),
    p_psi_banksy,
    width = 5.8,
    height = 4.8,
    dpi = 300
  )
  
  ## ----------------------------
  ## Top changed clusters
  ## Collapse Colon vs Liver/Lung
  ## ----------------------------
  
  psi_total_by_role_tbl <- psi_banksy_comp_tbl %>%
    distinct(sample_role, sample_id, Organ, total_all_PSI_cells) %>%
    group_by(sample_role) %>%
    summarise(
      total_all_PSI_cells = sum(total_all_PSI_cells),
      sample_id = paste(unique(sample_id), collapse = ";"),
      Organ = paste(unique(Organ), collapse = ";"),
      .groups = "drop"
    )
  
  psi_cluster_role_count_tbl <- psi_banksy_comp_tbl %>%
    group_by(banksy_cluster, sample_role) %>%
    summarise(
      n_cells = sum(n_cells),
      .groups = "drop"
    )
  
  psi_change_summary_tbl <- expand_grid(
    banksy_cluster = sort(unique(psi_banksy_comp_tbl$banksy_cluster)),
    sample_role = c("Primary", "Metastasis")
  ) %>%
    left_join(
      psi_cluster_role_count_tbl,
      by = c("banksy_cluster", "sample_role")
    ) %>%
    left_join(
      psi_total_by_role_tbl,
      by = "sample_role"
    ) %>%
    mutate(
      n_cells = replace_na(n_cells, 0L),
      pct = n_cells / total_all_PSI_cells * 100
    )
  
  psi_banksy_change_tbl <- psi_change_summary_tbl %>%
    select(
      banksy_cluster,
      sample_role,
      n_cells,
      total_all_PSI_cells,
      pct
    ) %>%
    pivot_wider(
      id_cols = banksy_cluster,
      names_from = sample_role,
      values_from = c(n_cells, total_all_PSI_cells, pct),
      names_glue = "{.value}_{sample_role}",
      values_fill = list(
        n_cells = 0,
        pct = 0
      )
    ) %>%
    mutate(
      delta_pct_met_vs_primary = pct_Metastasis - pct_Primary,
      abs_delta_pct = abs(delta_pct_met_vs_primary),
      max_pct = pmax(pct_Primary, pct_Metastasis, na.rm = TRUE),
      direction = case_when(
        delta_pct_met_vs_primary > 0 ~ "Increased in metastasis",
        delta_pct_met_vs_primary < 0 ~ "Decreased in metastasis",
        TRUE ~ "No change"
      )
    ) %>%
    filter(
      max_pct >= min_max_pct_for_top_change,
      abs_delta_pct >= min_abs_delta_pct_for_top_change
    ) %>%
    arrange(desc(abs_delta_pct))
  
  top_changed_clusters_tbl <- psi_banksy_change_tbl %>%
    slice_head(n = top_n_changed)
  
  if (!is.null(cluster_annotation_tbl)) {
    top_changed_clusters_tbl <- top_changed_clusters_tbl %>%
      left_join(
        cluster_annotation_tbl %>%
          mutate(banksy_cluster = as.character(banksy_cluster)) %>%
          select(
            banksy_cluster,
            dominant_cell_type,
            dominant_pct,
            top3_celltype_composition
          ),
        by = "banksy_cluster"
      )
  }
  
  write_csv(
    psi_banksy_comp_tbl,
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_composition_table.csv"))
  )
  
  write_csv(
    valid_banksy_clusters,
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_valid_clusters.csv"))
  )
  
  write_csv(
    psi_banksy_change_tbl,
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_cluster_change_table.csv"))
  )
  
  write_csv(
    top_changed_clusters_tbl,
    file.path(out_dir, paste0(case_id, "_PSI_banksy_", selected_res_label, "_top3_changed_clusters.csv"))
  )
  
  return(
    list(
      psi_banksy_comp_tbl = psi_banksy_comp_tbl,
      psi_banksy_change_tbl = psi_banksy_change_tbl,
      top_changed_clusters_tbl = top_changed_clusters_tbl,
      valid_banksy_clusters = valid_banksy_clusters,
      plot = p_psi_banksy
    )
  )
}

# ===== cell 7 =====
## ============================================================
## Run all selected paired cases
## ============================================================

all_case_results <- list()

for (i in seq_len(nrow(case_resolution_tbl))) {
  
  case_id <- case_resolution_tbl$case_id[i]
  selected_res_label <- case_resolution_tbl$selected_res_label[i]
  selected_banksy_col <- case_resolution_tbl$selected_banksy_col[i]
  
  message("============================================================")
  message("Running case: ", case_id, " | ", selected_res_label)
  message("============================================================")
  
  case_meta_banksy <- merge_case_metadata_banksy(case_id)
  
  heatmap_res <- make_case_celltype_enrichment_heatmap(
    case_meta_banksy = case_meta_banksy,
    case_id = case_id,
    selected_banksy_col = selected_banksy_col,
    selected_res_label = selected_res_label
  )
  
  psi_res <- make_case_PSI_barplot_and_top_changes(
    case_meta_banksy = case_meta_banksy,
    case_id = case_id,
    selected_banksy_col = selected_banksy_col,
    selected_res_label = selected_res_label,
    cluster_annotation_tbl = heatmap_res$cluster_annotation_tbl
  )
  
  all_case_results[[case_id]] <- list(
    case_meta_banksy = case_meta_banksy,
    heatmap_res = heatmap_res,
    psi_res = psi_res
  )
}

# ===== cell 8 =====
## ============================================================
## Combine top 3 changed clusters in each direction from all cases
## Top 3 increased in metastasis + top 3 decreased in metastasis
## ============================================================

all_top3_each_direction_changed_clusters_tbl <- map_dfr(
  names(all_case_results),
  function(case_id) {
    
    psi_res <- all_case_results[[case_id]]$psi_res
    
    if (is.null(psi_res)) {
      return(NULL)
    }
    
    change_tbl <- psi_res$psi_banksy_change_tbl
    
    if (is.null(change_tbl) || nrow(change_tbl) == 0) {
      return(NULL)
    }
    
    ## Optional: add cluster annotation if available
    cluster_annotation_tbl <- all_case_results[[case_id]]$heatmap_res$cluster_annotation_tbl
    
    if (!is.null(cluster_annotation_tbl)) {
      change_tbl <- change_tbl %>%
        left_join(
          cluster_annotation_tbl %>%
            mutate(banksy_cluster = as.character(banksy_cluster)) %>%
            select(
              banksy_cluster,
              dominant_cell_type,
              dominant_pct,
              top3_celltype_composition
            ),
          by = "banksy_cluster"
        )
    }
    
    increased_tbl <- change_tbl %>%
      filter(delta_pct_met_vs_primary > 0) %>%
      arrange(desc(delta_pct_met_vs_primary)) %>%
      slice_head(n = 3) %>%
      mutate(
        change_direction = "Increased_in_metastasis",
        direction_rank = row_number()
      )
    
    decreased_tbl <- change_tbl %>%
      filter(delta_pct_met_vs_primary < 0) %>%
      arrange(delta_pct_met_vs_primary) %>%
      slice_head(n = 3) %>%
      mutate(
        change_direction = "Decreased_in_metastasis",
        direction_rank = row_number()
      )
    
    bind_rows(increased_tbl, decreased_tbl) %>%
      mutate(
        case_id = case_id,
        selected_res_label = case_resolution_tbl$selected_res_label[
          match(case_id, case_resolution_tbl$case_id)
        ],
        .before = 1
      )
  }
)

write_csv(
  all_top3_each_direction_changed_clusters_tbl,
  file.path(out_root, "all_cases_top3_each_direction_changed_PSI_banksy_clusters.csv")
)

head(all_top3_each_direction_changed_clusters_tbl)

# ===== cell 10 =====
final_cluster_list_tbl <- all_top3_each_direction_changed_clusters_tbl %>%
  mutate(
    banksy_cluster = as.character(banksy_cluster),
    case_cluster_id = paste(case_id, selected_res_label, banksy_cluster, sep = "__"),
    cluster_label = paste0(case_id, " | ", selected_res_label, " | C", banksy_cluster)
  )

# ===== cell 11 =====
## ============================================================
## Cell-type composition of final selected changed clusters
## ============================================================

cell_type_order <- preferred_cell_type_order

get_final_cluster_composition <- function(case_id) {
  
  selected_res_label <- case_resolution_tbl$selected_res_label[
    match(case_id, case_resolution_tbl$case_id)
  ]
  
  selected_banksy_col <- case_resolution_tbl$selected_banksy_col[
    match(case_id, case_resolution_tbl$case_id)
  ]
  
  case_meta_banksy <- all_case_results[[case_id]]$case_meta_banksy
  
  target_clusters <- final_cluster_list_tbl %>%
    filter(case_id == .env$case_id) %>%
    pull(banksy_cluster) %>%
    unique()
  
  comp_tbl <- case_meta_banksy %>%
    filter(!is.na(.data[[selected_banksy_col]])) %>%
    filter(!is.na(Broad_cell_type1), Broad_cell_type1 != "") %>%
    filter(Broad_cell_type1 != "Doublet") %>%
    mutate(
      case_id = .env$case_id,
      selected_res_label = .env$selected_res_label,
      banksy_cluster = as.character(.data[[selected_banksy_col]]),
      case_cluster_id = paste(case_id, selected_res_label, banksy_cluster, sep = "__")
    ) %>%
    filter(banksy_cluster %in% target_clusters) %>%
    count(
      case_id,
      selected_res_label,
      banksy_cluster,
      case_cluster_id,
      Broad_cell_type1,
      name = "n_cells"
    ) %>%
    group_by(case_id, selected_res_label, banksy_cluster, case_cluster_id) %>%
    mutate(
      total_cluster_cells = sum(n_cells),
      frac_celltype = n_cells / total_cluster_cells,
      pct_celltype = frac_celltype * 100
    ) %>%
    ungroup()
  
  return(comp_tbl)
}

final_cluster_celltype_comp_tbl <- map_dfr(
  unique(final_cluster_list_tbl$case_id),
  get_final_cluster_composition
)

final_cluster_celltype_comp_tbl

# ===== cell 12 =====


p_direction_change <- final_cluster_list_tbl %>%
  mutate(
    cluster_label = paste0(
      case_id, " | C", banksy_cluster,
      " | ", selected_res_label,
      "\n", dominant_cell_type
    ),
    cluster_label = fct_reorder(cluster_label, delta_pct_met_vs_primary)
  ) %>%
  ggplot(
    aes(
      x = delta_pct_met_vs_primary,
      y = cluster_label,
      fill = change_direction
    )
  ) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Top changed PSI-associated BANKSY clusters",
    x = "Change in PSI proportion: metastasis − primary (%)",
    y = NULL,
    fill = "Direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_direction_change

# ===== cell 13 =====
ggsave(
  file.path(out_root, "all_cases_top3_each_direction_delta_barplot.pdf"),
  p_direction_change,
  width = 8,
  height = 8
)

ggsave(
  file.path(out_root, "all_cases_top3_each_direction_delta_barplot.png"),
  p_direction_change,
  width = 8,
  height = 8,
  dpi = 300
)

# ===== cell 14 =====
## ============================================================
## Matrix: rows = changed BANKSY clusters, columns = cell types
## ============================================================

composition_wide_tbl <- final_cluster_celltype_comp_tbl %>%
  select(case_cluster_id, Broad_cell_type1, pct_celltype) %>%
  pivot_wider(
    names_from = Broad_cell_type1,
    values_from = pct_celltype,
    values_fill = 0
  )

present_cell_types <- setdiff(colnames(composition_wide_tbl), "case_cluster_id")

cell_types_for_heatmap <- c(
  cell_type_order[cell_type_order %in% present_cell_types],
  setdiff(present_cell_types, cell_type_order)
)

composition_mat <- composition_wide_tbl %>%
  select(case_cluster_id, all_of(cell_types_for_heatmap)) %>%
  column_to_rownames("case_cluster_id") %>%
  as.matrix()

## Make sure row metadata matches matrix order
row_meta_tbl <- final_cluster_list_tbl %>%
  distinct(
    case_cluster_id,
    case_id,
    selected_res_label,
    banksy_cluster,
    cluster_label,
    change_direction,
    direction_rank,
    delta_pct_met_vs_primary,
    abs_delta_pct,
    dominant_cell_type,
    top3_celltype_composition
  ) %>%
  filter(case_cluster_id %in% rownames(composition_mat)) %>%
  arrange(match(case_cluster_id, rownames(composition_mat)))

composition_mat <- composition_mat[row_meta_tbl$case_cluster_id, , drop = FALSE]

rownames(composition_mat) <- row_meta_tbl$cluster_label

# ===== cell 15 =====
## ============================================================
## Heatmap with row annotations
## ============================================================

direction_col <- c(
  "Increased_in_metastasis" = "#D73027",
  "Decreased_in_metastasis" = "#4575B4"
)

case_col <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl$case_id))),
  unique(row_meta_tbl$case_id)
)

dominant_celltype_col <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl$dominant_cell_type))),
  unique(row_meta_tbl$dominant_cell_type)
)

row_ha <- rowAnnotation(
  case = row_meta_tbl$case_id,
  direction = row_meta_tbl$change_direction,
  dominant = row_meta_tbl$dominant_cell_type,
  delta = anno_barplot(
    row_meta_tbl$delta_pct_met_vs_primary,
    baseline = 0,
    border = FALSE,
    gp = gpar(fill = ifelse(
      row_meta_tbl$delta_pct_met_vs_primary > 0,
      "#D73027",
      "#4575B4"
    )),
    width = unit(2.5, "cm")
  ),
  col = list(
    case = case_col,
    direction = direction_col,
    dominant = dominant_celltype_col
  )
)

cell_pct_col_fun <- colorRamp2(
  c(0, 25, 50, 75, 100),
  c("white", "#FEE8C8", "#FDBB84", "#E34A33", "#7F0000")
)

ht_comp <- Heatmap(
  composition_mat,
  name = "Cell type %",
  col = cell_pct_col_fun,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Cell type %",
    at = c(0, 25, 50, 75, 100)
  )
)

draw(ht_comp)

# ===== cell 16 =====
pdf(
  file.path(out_root, "all_cases_top3_each_direction_celltype_composition_heatmap.pdf"),
  width = 11,
  height = 9
)

draw(ht_comp)

dev.off()

png(
  file.path(out_root, "all_cases_top3_each_direction_celltype_composition_heatmap.png"),
  width = 3300,
  height = 2700,
  res = 300
)

draw(ht_comp)

dev.off()

# ===== cell 17 =====
## ============================================================
## Scaled cell-type composition heatmap
## Scale by cell type = column-wise z-score
## Rows = changed BANKSY clusters
## Columns = Broad_cell_type1
## ============================================================

scale_mat_by_column <- function(mat) {
  
  scaled_mat <- apply(
    mat,
    2,
    function(x) {
      
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      
      if (is.na(x_sd) || x_sd == 0) {
        return(rep(0, length(x)))
      }
      
      (x - x_mean) / x_sd
    }
  )
  
  scaled_mat <- as.matrix(scaled_mat)
  rownames(scaled_mat) <- rownames(mat)
  colnames(scaled_mat) <- colnames(mat)
  
  return(scaled_mat)
}

composition_mat_scaled <- scale_mat_by_column(composition_mat)

## Cap z-score for visualization
composition_mat_scaled_capped <- composition_mat_scaled
composition_mat_scaled_capped[composition_mat_scaled_capped > 3] <- 3
composition_mat_scaled_capped[composition_mat_scaled_capped < -3] <- -3

# ===== cell 18 =====
## ============================================================
## Heatmap color scale
## ============================================================

celltype_z_col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

ht_comp_scaled <- Heatmap(
  composition_mat_scaled_capped,
  name = "Cell type\nz-score",
  col = celltype_z_col_fun,
  right_annotation = row_ha,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Cell type\nz-score",
    at = c(-3, 0, 3),
    labels = c("≤ -3", "0", "≥ 3")
  )
)

draw(ht_comp_scaled)

# ===== cell 19 =====
pdf(
  file.path(out_root, "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_by_celltype.pdf"),
  width = 11,
  height = 9
)

draw(ht_comp_scaled)

dev.off()

png(
  file.path(out_root, "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_by_celltype.png"),
  width = 3300,
  height = 2700,
  res = 300
)

draw(ht_comp_scaled)

dev.off()

# ===== cell 20 =====
## ============================================================
## Simplified scaled heatmap:
## - exclude CM626C
## - remove Intestinal_epithelium and Enteric_neuron
## - rescale by cell type across clusters
## ============================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

exclude_case <- "CM626C"

exclude_cell_types <- c(
  "Intestinal_epithelium",
  "Enteric_neuron",
  "Enteric neuron"
)

## ----------------------------
## Filter final cluster list
## ----------------------------

final_cluster_list_tbl_simplified <- final_cluster_list_tbl %>%
  filter(case_id != exclude_case)

## ----------------------------
## Filter composition table
## Then renormalize among remaining cell types
## ----------------------------

final_cluster_celltype_comp_tbl_simplified <- final_cluster_celltype_comp_tbl %>%
  semi_join(
    final_cluster_list_tbl_simplified %>%
      select(case_cluster_id) %>%
      distinct(),
    by = "case_cluster_id"
  ) %>%
  filter(!Broad_cell_type1 %in% exclude_cell_types) %>%
  group_by(case_cluster_id) %>%
  mutate(
    total_remaining_cells = sum(n_cells),
    frac_celltype_remaining = n_cells / total_remaining_cells,
    pct_celltype_remaining = frac_celltype_remaining * 100
  ) %>%
  ungroup() %>%
  filter(total_remaining_cells > 0)

## Check remaining cell types
final_cluster_celltype_comp_tbl_simplified %>%
  count(Broad_cell_type1) %>%
  arrange(desc(n))

# ===== cell 21 =====
## ============================================================
## Matrix: rows = changed clusters, columns = remaining cell types
## ============================================================

composition_wide_tbl_simplified <- final_cluster_celltype_comp_tbl_simplified %>%
  select(case_cluster_id, Broad_cell_type1, pct_celltype_remaining) %>%
  pivot_wider(
    names_from = Broad_cell_type1,
    values_from = pct_celltype_remaining,
    values_fill = 0
  )

present_cell_types_simplified <- setdiff(
  colnames(composition_wide_tbl_simplified),
  "case_cluster_id"
)

cell_types_for_heatmap_simplified <- c(
  cell_type_order[cell_type_order %in% present_cell_types_simplified],
  setdiff(present_cell_types_simplified, cell_type_order)
)

composition_mat_simplified <- composition_wide_tbl_simplified %>%
  select(case_cluster_id, all_of(cell_types_for_heatmap_simplified)) %>%
  column_to_rownames("case_cluster_id") %>%
  as.matrix()

# ===== cell 22 =====
## ============================================================
## Row metadata
## Recompute dominant cell type after removing intestinal/enteric cells
## ============================================================

dominant_remaining_tbl <- final_cluster_celltype_comp_tbl_simplified %>%
  group_by(case_cluster_id) %>%
  arrange(desc(pct_celltype_remaining), .by_group = TRUE) %>%
  summarise(
    dominant_remaining_cell_type = first(Broad_cell_type1),
    dominant_remaining_pct = first(pct_celltype_remaining),
    top3_remaining_celltype_composition = paste0(
      head(Broad_cell_type1, 3),
      " (",
      round(head(pct_celltype_remaining, 3), 1),
      "%)",
      collapse = "; "
    ),
    .groups = "drop"
  )

row_meta_tbl_simplified <- final_cluster_list_tbl_simplified %>%
  distinct(
    case_cluster_id,
    case_id,
    selected_res_label,
    banksy_cluster,
    cluster_label,
    change_direction,
    direction_rank,
    delta_pct_met_vs_primary,
    abs_delta_pct
  ) %>%
  filter(case_cluster_id %in% rownames(composition_mat_simplified)) %>%
  left_join(
    dominant_remaining_tbl,
    by = "case_cluster_id"
  ) %>%
  arrange(match(case_cluster_id, rownames(composition_mat_simplified)))

composition_mat_simplified <- composition_mat_simplified[
  row_meta_tbl_simplified$case_cluster_id,
  ,
  drop = FALSE
]

rownames(composition_mat_simplified) <- row_meta_tbl_simplified$cluster_label

# ===== cell 23 =====
## ============================================================
## Column-wise scaling
## Each cell type is scaled across selected clusters
## ============================================================

scale_mat_by_column <- function(mat) {
  
  scaled_mat <- apply(
    mat,
    2,
    function(x) {
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      
      if (is.na(x_sd) || x_sd == 0) {
        return(rep(0, length(x)))
      }
      
      (x - x_mean) / x_sd
    }
  )
  
  scaled_mat <- as.matrix(scaled_mat)
  rownames(scaled_mat) <- rownames(mat)
  colnames(scaled_mat) <- colnames(mat)
  
  return(scaled_mat)
}

composition_mat_scaled_simplified <- scale_mat_by_column(composition_mat_simplified)

composition_mat_scaled_capped_simplified <- composition_mat_scaled_simplified
composition_mat_scaled_capped_simplified[composition_mat_scaled_capped_simplified > 3] <- 3
composition_mat_scaled_capped_simplified[composition_mat_scaled_capped_simplified < -3] <- -3

# ===== cell 24 =====
## ============================================================
## Row annotations
## ============================================================

direction_values <- unique(row_meta_tbl_simplified$change_direction)

direction_col_simplified <- setNames(
  ifelse(
    str_detect(direction_values, "Increased"),
    "#D73027",
    "#4575B4"
  ),
  direction_values
)

case_col_simplified <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_simplified$case_id))),
  unique(row_meta_tbl_simplified$case_id)
)

dominant_celltype_col_simplified <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_simplified$dominant_remaining_cell_type))),
  unique(row_meta_tbl_simplified$dominant_remaining_cell_type)
)

row_ha_simplified <- rowAnnotation(
  case = row_meta_tbl_simplified$case_id,
  direction = row_meta_tbl_simplified$change_direction,
  dominant_remaining = row_meta_tbl_simplified$dominant_remaining_cell_type,
  delta = anno_barplot(
    row_meta_tbl_simplified$delta_pct_met_vs_primary,
    baseline = 0,
    border = FALSE,
    gp = gpar(
      fill = ifelse(
        row_meta_tbl_simplified$delta_pct_met_vs_primary > 0,
        "#D73027",
        "#4575B4"
      )
    ),
    width = unit(2.5, "cm")
  ),
  col = list(
    case = case_col_simplified,
    direction = direction_col_simplified,
    dominant_remaining = dominant_celltype_col_simplified
  )
)

# ===== cell 25 =====
## ============================================================
## Simplified scaled heatmap
## ============================================================

celltype_z_col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

ht_comp_scaled_simplified <- Heatmap(
  composition_mat_scaled_capped_simplified,
  name = "Cell type\nz-score",
  col = celltype_z_col_fun,
  right_annotation = row_ha_simplified,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  column_title = "Top changed PSI-associated BANKSY clusters, excluding CM626C and intestinal/enteric cells",
  heatmap_legend_param = list(
    title = "Cell type\nz-score",
    at = c(-3, 0, 3),
    labels = c("≤ -3", "0", "≥ 3")
  )
)

draw(ht_comp_scaled_simplified)

# ===== cell 26 =====
## ============================================================
## Save simplified heatmap and matrices
## ============================================================

pdf(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_no_CM626C_no_intestinal_enteric.pdf"
  ),
  width = 11,
  height = 9
)

draw(ht_comp_scaled_simplified)

dev.off()

png(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_no_CM626C_no_intestinal_enteric.png"
  ),
  width = 3300,
  height = 2700,
  res = 300
)

draw(ht_comp_scaled_simplified)

dev.off()

write_csv(
  as.data.frame(composition_mat_simplified) %>%
    rownames_to_column("cluster_label"),
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_matrix_no_CM626C_no_intestinal_enteric.csv"
  )
)

write_csv(
  as.data.frame(composition_mat_scaled_simplified) %>%
    rownames_to_column("cluster_label"),
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_matrix_scaled_no_CM626C_no_intestinal_enteric.csv"
  )
)

write_csv(
  row_meta_tbl_simplified,
  file.path(
    out_root,
    "all_cases_top3_each_direction_row_metadata_no_CM626C_no_intestinal_enteric.csv"
  )
)

# ===== cell 28 =====
## ============================================================
## Make scaled cell-type composition heatmap excluding CM626C
## ============================================================

exclude_case <- "CM626C"

final_cluster_list_tbl_no_CM626C <- final_cluster_list_tbl %>%
  filter(case_id != exclude_case)

final_cluster_celltype_comp_tbl_no_CM626C <- final_cluster_celltype_comp_tbl %>%
  semi_join(
    final_cluster_list_tbl_no_CM626C %>%
      select(case_cluster_id) %>%
      distinct(),
    by = "case_cluster_id"
  )

# ===== cell 29 =====
## ============================================================
## Matrix: rows = changed clusters, columns = cell types
## excluding CM626C
## ============================================================

composition_wide_tbl_no_CM626C <- final_cluster_celltype_comp_tbl_no_CM626C %>%
  select(case_cluster_id, Broad_cell_type1, pct_celltype) %>%
  pivot_wider(
    names_from = Broad_cell_type1,
    values_from = pct_celltype,
    values_fill = 0
  )

present_cell_types_no_CM626C <- setdiff(
  colnames(composition_wide_tbl_no_CM626C),
  "case_cluster_id"
)

cell_types_for_heatmap_no_CM626C <- c(
  cell_type_order[cell_type_order %in% present_cell_types_no_CM626C],
  setdiff(present_cell_types_no_CM626C, cell_type_order)
)

composition_mat_no_CM626C <- composition_wide_tbl_no_CM626C %>%
  select(case_cluster_id, all_of(cell_types_for_heatmap_no_CM626C)) %>%
  column_to_rownames("case_cluster_id") %>%
  as.matrix()

row_meta_tbl_no_CM626C <- final_cluster_list_tbl_no_CM626C %>%
  distinct(
    case_cluster_id,
    case_id,
    selected_res_label,
    banksy_cluster,
    cluster_label,
    change_direction,
    direction_rank,
    delta_pct_met_vs_primary,
    abs_delta_pct,
    dominant_cell_type,
    top3_celltype_composition
  ) %>%
  filter(case_cluster_id %in% rownames(composition_mat_no_CM626C)) %>%
  arrange(match(case_cluster_id, rownames(composition_mat_no_CM626C))) %>%
  mutate(
    dominant_cell_type = replace_na(dominant_cell_type, "Unknown")
  )

composition_mat_no_CM626C <- composition_mat_no_CM626C[
  row_meta_tbl_no_CM626C$case_cluster_id,
  ,
  drop = FALSE
]

rownames(composition_mat_no_CM626C) <- row_meta_tbl_no_CM626C$cluster_label

# ===== cell 30 =====
## ============================================================
## Column-wise scaling: each cell type scaled across clusters
## ============================================================

scale_mat_by_column <- function(mat) {
  
  scaled_mat <- apply(
    mat,
    2,
    function(x) {
      x_mean <- mean(x, na.rm = TRUE)
      x_sd <- sd(x, na.rm = TRUE)
      
      if (is.na(x_sd) || x_sd == 0) {
        return(rep(0, length(x)))
      }
      
      (x - x_mean) / x_sd
    }
  )
  
  scaled_mat <- as.matrix(scaled_mat)
  rownames(scaled_mat) <- rownames(mat)
  colnames(scaled_mat) <- colnames(mat)
  
  return(scaled_mat)
}

composition_mat_scaled_no_CM626C <- scale_mat_by_column(composition_mat_no_CM626C)

## Cap z-score for visualization
composition_mat_scaled_capped_no_CM626C <- composition_mat_scaled_no_CM626C
composition_mat_scaled_capped_no_CM626C[composition_mat_scaled_capped_no_CM626C > 3] <- 3
composition_mat_scaled_capped_no_CM626C[composition_mat_scaled_capped_no_CM626C < -3] <- -3

# ===== cell 31 =====
## ============================================================
## Row annotations
## ============================================================

direction_col <- c(
  "Increased_in_metastasis" = "#D73027",
  "Decreased_in_metastasis" = "#4575B4"
)

case_col_no_CM626C <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_no_CM626C$case_id))),
  unique(row_meta_tbl_no_CM626C$case_id)
)

dominant_celltype_col_no_CM626C <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_no_CM626C$dominant_cell_type))),
  unique(row_meta_tbl_no_CM626C$dominant_cell_type)
)

row_ha_no_CM626C <- rowAnnotation(
  case = row_meta_tbl_no_CM626C$case_id,
  direction = row_meta_tbl_no_CM626C$change_direction,
  dominant = row_meta_tbl_no_CM626C$dominant_cell_type,
  delta = anno_barplot(
    row_meta_tbl_no_CM626C$delta_pct_met_vs_primary,
    baseline = 0,
    border = FALSE,
    gp = gpar(
      fill = ifelse(
        row_meta_tbl_no_CM626C$delta_pct_met_vs_primary > 0,
        "#D73027",
        "#4575B4"
      )
    ),
    width = unit(2.5, "cm")
  ),
  col = list(
    case = case_col_no_CM626C,
    direction = direction_col,
    dominant = dominant_celltype_col_no_CM626C
  )
)

# ===== cell 32 =====
## ============================================================
## Scaled heatmap excluding CM626C
## ============================================================

celltype_z_col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

ht_comp_scaled_no_CM626C <- Heatmap(
  composition_mat_scaled_capped_no_CM626C,
  name = "Cell type\nz-score",
  col = celltype_z_col_fun,
  right_annotation = row_ha_no_CM626C,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Cell type\nz-score",
    at = c(-3, 0, 3),
    labels = c("≤ -3", "0", "≥ 3")
  ),
  column_title = "Top changed PSI-associated BANKSY clusters, excluding CM626C"
)

draw(ht_comp_scaled_no_CM626C)

# ===== cell 33 =====
pdf(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_by_celltype_no_CM626C.pdf"
  ),
  width = 11,
  height = 9
)

draw(ht_comp_scaled_no_CM626C)

dev.off()

png(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_scaled_by_celltype_no_CM626C.png"
  ),
  width = 3300,
  height = 2700,
  res = 300
)

draw(ht_comp_scaled_no_CM626C)

dev.off()

# ===== cell 34 =====
write_csv(
  as.data.frame(composition_mat_scaled_no_CM626C) %>%
    rownames_to_column("cluster_label"),
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_matrix_scaled_by_celltype_no_CM626C.csv"
  )
)

# ===== cell 37 =====
## ============================================================
## Matrix: rows = changed BANKSY clusters, columns = cell types
## Exclude CM626C; remove Intestinal_epithelium and Enteric_neuron
## Heatmap scale capped at 0–50%
## ============================================================

exclude_case_id <- "CM626C"
exclude_cell_types <- c("Intestinal_epithelium", "Enteric_neuron")

## Define clusters to keep
keep_cluster_ids <- final_cluster_list_tbl %>%
  filter(case_id != exclude_case_id) %>%
  pull(case_cluster_id) %>%
  unique()

composition_wide_tbl_no626_no_epi_neuron <- final_cluster_celltype_comp_tbl %>%
  filter(
    case_cluster_id %in% keep_cluster_ids,
    !Broad_cell_type1 %in% exclude_cell_types
  ) %>%
  select(case_cluster_id, Broad_cell_type1, pct_celltype) %>%
  pivot_wider(
    names_from = Broad_cell_type1,
    values_from = pct_celltype,
    values_fill = 0
  )

present_cell_types <- setdiff(
  colnames(composition_wide_tbl_no626_no_epi_neuron),
  "case_cluster_id"
)

cell_types_for_heatmap <- c(
  cell_type_order[cell_type_order %in% present_cell_types],
  setdiff(present_cell_types, cell_type_order)
)

composition_mat_no626_no_epi_neuron <- composition_wide_tbl_no626_no_epi_neuron %>%
  select(case_cluster_id, all_of(cell_types_for_heatmap)) %>%
  column_to_rownames("case_cluster_id") %>%
  as.matrix()

## Make sure row metadata matches matrix order
row_meta_tbl_no626_no_epi_neuron <- final_cluster_list_tbl %>%
  filter(case_id != exclude_case_id) %>%
  distinct(
    case_cluster_id,
    case_id,
    selected_res_label,
    banksy_cluster,
    cluster_label,
    change_direction,
    direction_rank,
    delta_pct_met_vs_primary,
    abs_delta_pct,
    dominant_cell_type,
    top3_celltype_composition
  ) %>%
  filter(case_cluster_id %in% rownames(composition_mat_no626_no_epi_neuron)) %>%
  arrange(match(case_cluster_id, rownames(composition_mat_no626_no_epi_neuron)))

composition_mat_no626_no_epi_neuron <- composition_mat_no626_no_epi_neuron[
  row_meta_tbl_no626_no_epi_neuron$case_cluster_id,
  ,
  drop = FALSE
]

rownames(composition_mat_no626_no_epi_neuron) <- row_meta_tbl_no626_no_epi_neuron$cluster_label

## Cap values at 50 for visualization
composition_mat_no626_no_epi_neuron_plot <- pmin(
  composition_mat_no626_no_epi_neuron,
  50
)

## ============================================================
## Heatmap with row annotations
## ============================================================

direction_col <- c(
  "Increased_in_metastasis" = "#D73027",
  "Decreased_in_metastasis" = "#4575B4"
)

case_col <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_no626_no_epi_neuron$case_id))),
  unique(row_meta_tbl_no626_no_epi_neuron$case_id)
)

dominant_celltype_col <- setNames(
  grDevices::rainbow(length(unique(row_meta_tbl_no626_no_epi_neuron$dominant_cell_type))),
  unique(row_meta_tbl_no626_no_epi_neuron$dominant_cell_type)
)

row_ha_no626_no_epi_neuron <- rowAnnotation(
  case = row_meta_tbl_no626_no_epi_neuron$case_id,
  direction = row_meta_tbl_no626_no_epi_neuron$change_direction,
  dominant = row_meta_tbl_no626_no_epi_neuron$dominant_cell_type,
  delta = anno_barplot(
    row_meta_tbl_no626_no_epi_neuron$delta_pct_met_vs_primary,
    baseline = 0,
    border = FALSE,
    gp = gpar(fill = ifelse(
      row_meta_tbl_no626_no_epi_neuron$delta_pct_met_vs_primary > 0,
      "#D73027",
      "#4575B4"
    )),
    width = unit(2.5, "cm")
  ),
  col = list(
    case = case_col,
    direction = direction_col,
    dominant = dominant_celltype_col
  )
)

cell_pct_col_fun_0_50 <- colorRamp2(
  c(0, 12.5, 25, 37.5, 50),
  c("white", "#FEE8C8", "#FDBB84", "#E34A33", "#7F0000")
)

ht_comp_no626_no_epi_neuron <- Heatmap(
  composition_mat_no626_no_epi_neuron_plot,
  name = "Cell type %",
  col = cell_pct_col_fun_0_50,
  right_annotation = row_ha_no626_no_epi_neuron,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(
    title = "Cell type %",
    at = c(0, 12.5, 25, 37.5, 50),
    labels = c("0", "12.5", "25", "37.5", "≥50")
  )
)

draw(ht_comp_no626_no_epi_neuron)

# ===== cell 38 =====
pdf(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_proportion_by_celltype_no_CM626C.pdf"
  ),
  width = 11,
  height = 9
)

draw(ht_comp_no626_no_epi_neuron)

dev.off()

png(
  file.path(
    out_root,
    "all_cases_top3_each_direction_celltype_composition_heatmap_proportion_by_celltype_no_CM626C.png"
  ),
  width = 3300,
  height = 2700,
  res = 300
)

draw(ht_comp_no626_no_epi_neuron)

dev.off()
