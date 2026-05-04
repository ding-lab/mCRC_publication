# Standalone: per-sample tumor vs CAF proportions — Spearman correlation heatmaps by organ group.

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

output_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Neuroendocrine_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

updated_cell_file <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean3_metadata_cell_type_all_20260420.csv"
updated_cell <- read.table(updated_cell_file, header = TRUE, sep = ",", row.names = 1)
head(updated_cell, 3)

metadata_name <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean6_metadata.csv"
metadata <- read.csv(metadata_name, row.names = 1) %>%
  dplyr::select(-c(broad_cell_type, cell_type_all4))

idx <- match(rownames(metadata), rownames(updated_cell))
metadata$broad_cell_type <- updated_cell$broad_cell_type[idx]
metadata$cell_type_all4 <- updated_cell$cell_type_all4[idx]

tumor_fib_metadata <- metadata %>%
  dplyr::filter(cell_type_all4 %in% c("mCAF", "WNT5A_BMP", "WNT5A_infl", "iCAF") | broad_cell_type == "Tumor")

tumor_fib_count <- tumor_fib_metadata %>%
  group_by(orig.ident) %>%
  summarise(
    all_cell_count = n(),
    tumor_count = sum(broad_cell_type == "Tumor", na.rm = TRUE),
    caf_count = sum(cell_type_all4 %in% c("mCAF", "WNT5A_BMP", "WNT5A_infl", "iCAF"), na.rm = TRUE),
    CMETS_count = sum(cell_type_all4 == "CMETS", na.rm = TRUE),
    NE_count = sum(cell_type_all3 == "Neuroendocrine-like tumor", na.rm = TRUE),
    APCDD1_tumor_count = sum(cell_type_all2 == "APCDD1_CRC", na.rm = TRUE),
    can_stem_tumor_count = sum(cell_type_all2 == "Canonical_CRC_Stem", na.rm = TRUE),
    can_stem_prolif_count = sum(cell_type_all2 == "Canonical_CRC_Stem_Proliferation", na.rm = TRUE),
    can_intestine_tumor_count = sum(cell_type_all2 == "Canonical_CRC_Intestine", na.rm = TRUE),
    can_intestine_prolif_count = sum(cell_type_all2 == "Canonical_CRC_Intestine_Proliferation", na.rm = TRUE),
    mCAF_count = sum(cell_type_all2 == "mCAF", na.rm = TRUE),
    WNT5A_BMP_count = sum(cell_type_all2 == "WNT5A_BMP", na.rm = TRUE),
    WNT5A_infl_count = sum(cell_type_all2 == "WNT5A_infl", na.rm = TRUE),
    iCAF_count = sum(cell_type_all2 == "iCAF", na.rm = TRUE),
    Organ = dplyr::first(Site_of_Origin),
    Patient_ID = dplyr::first(Patient_ID),
    Primary_Side = dplyr::first(Primary_Side),
    Tx_in_6mo = dplyr::first(Tx_in_6mo),
    .groups = "drop"
  ) %>%
  filter(tumor_count > 20, caf_count > 10) %>%
  mutate(orig.ident = str_replace(orig.ident, "^HT413C1-Th1K[24]A2Nd1_2Bma1_1$", "HT413C1-Th1"))

unique(tumor_fib_metadata$cell_type_all4)

tumor_caf_prop <- tumor_fib_count %>%
  group_by(orig.ident) %>%
  summarise(
    across(
      c(
        all_cell_count, tumor_count, caf_count,
        CMETS_count, NE_count, APCDD1_tumor_count,
        can_stem_tumor_count, can_stem_prolif_count,
        can_intestine_tumor_count, can_intestine_prolif_count,
        mCAF_count, WNT5A_BMP_count, WNT5A_infl_count,
        iCAF_count
      ),
      \(x) sum(x, na.rm = TRUE)
    ),
    tumor_prop = round(100 * tumor_count / all_cell_count, 2),
    caf_prop = round(100 * caf_count / all_cell_count, 2),
    CMETS_prop = round(100 * CMETS_count / tumor_count, 2),
    NE_prop = round(100 * NE_count / tumor_count, 2),
    can_stem_tumor_prop = round(100 * can_stem_tumor_count / tumor_count, 2),
    can_stem_prolif_prop = round(100 * can_stem_prolif_count / tumor_count, 2),
    can_intestine_tumor_prop = round(100 * can_intestine_tumor_count / tumor_count, 2),
    can_intestine_prolif_prop = round(100 * can_intestine_prolif_count / tumor_count, 2),
    APCDD1_tumor_prop = round(100 * APCDD1_tumor_count / tumor_count, 2),
    can_tumor_prop = round(100 * (can_stem_tumor_count + can_intestine_tumor_count) / tumor_count, 2),
    can_prolif_prop = round(100 * (can_stem_prolif_count + can_intestine_prolif_count) / tumor_count, 2),
    mCAF_prop = round(100 * mCAF_count / caf_count, 2),
    WNT5A_BMP_prop = round(100 * WNT5A_BMP_count / caf_count, 2),
    WNT5A_infl_prop = round(100 * WNT5A_infl_count / caf_count, 2),
    iCAF_prop = round(100 * iCAF_count / caf_count, 2),
    CMETS_all_prop = round(100 * CMETS_count / all_cell_count, 2),
    NE_all_prop = round(100 * NE_count / all_cell_count, 2),
    can_stem_tumor_all_prop = round(100 * can_stem_tumor_count / all_cell_count, 2),
    can_intestine_tumor_all_prop = round(100 * can_intestine_tumor_count / all_cell_count, 2),
    can_intestine_prolif_all_prop = round(100 * can_intestine_prolif_count / all_cell_count, 2),
    mCAF_all_prop = round(100 * mCAF_count / all_cell_count, 2),
    WNT5A_BMP_all_prop = round(100 * WNT5A_BMP_count / all_cell_count, 2),
    WNT5A_infl_all_prop = round(100 * WNT5A_infl_count / all_cell_count, 2),
    iCAF_all_prop = round(100 * iCAF_count / all_cell_count, 2),
    Organ = dplyr::first(Organ),
    Patient_ID = dplyr::first(Patient_ID),
    Primary_Side = dplyr::first(Primary_Side),
    Tx_in_6mo = dplyr::first(Tx_in_6mo),
    .groups = "drop"
  )

colon_summary <- tumor_caf_prop %>% dplyr::filter(Organ == "colon")
meta_summary <- tumor_caf_prop %>% dplyr::filter(Organ %in% c("liver", "lung"))

organ_list2 <- list(
  colon = colon_summary,
  metastasis = meta_summary
)

x_axis <- c("mCAF_prop", "WNT5A_infl_prop", "WNT5A_BMP_prop", "iCAF_prop")
y_axis <- c(
  "CMETS_prop", "can_stem_tumor_prop",
  "can_intestine_prolif_prop", "can_intestine_tumor_prop"
)

compute_corr_matrix <- function(df, x_vars, y_vars) {
  corr_mat <- matrix(
    NA_real_,
    nrow = length(y_vars),
    ncol = length(x_vars),
    dimnames = list(y_vars, x_vars)
  )
  p_mat <- corr_mat

  for (y in y_vars) {
    for (x in x_vars) {
      tmp <- cor.test(df[[x]], df[[y]], method = "spearman")
      corr_mat[y, x] <- unname(tmp$estimate)
      p_mat[y, x] <- tmp$p.value
    }
  }
  list(corr = corr_mat, p = p_mat)
}

p_to_stars <- function(p) {
  ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", ""))
}

col_fun <- colorRamp2(
  c(-0.75, -0.5, 0, 0.5, 0.75),
  c("#2166ac", "#67a9cf", "white", "#f4a582", "#ca0020")
)

for (organ_name in names(organ_list2)) {
  df <- organ_list2[[organ_name]]
  results <- compute_corr_matrix(df, x_axis, y_axis)
  corr_mat <- results$corr
  p_mat <- results$p
  star_mat <- matrix(
    p_to_stars(p_mat),
    nrow = nrow(p_mat),
    dimnames = dimnames(p_mat)
  )

  ht <- Heatmap(
    corr_mat,
    name = paste0(organ_name, "_cor"),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = gpar(col = "black", lwd = 1.2),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(star_mat[i, j], x, y, gp = gpar(fontsize = 14, fontface = "bold"))
    },
    column_title = paste0(toupper(organ_name), " – CAF vs Tumor Subtype Correlations"),
    heatmap_legend_param = list(title = "Spearman Rho")
  )

  pdf(file.path(output_dir, paste0("corr_heatmap_", organ_name, "_combined.pdf")), width = 6, height = 5)
  draw(ht)
  dev.off()
}
