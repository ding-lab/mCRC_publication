library(tidyverse)
library(patchwork)
library(cowplot)
library(clinfun)    
library(scales)
library(ComplexHeatmap)
library(circlize)
library(grid)

out_dir = getwd()
data_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/'
metadata_name = '57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean6_metadata.csv'
metadata = read.csv(file.path(data_dir, metadata_name), row.name = 1)
colnames(metadata)

tumor_metadata <- metadata %>% dplyr::filter(broad_cell_type == 'Tumor')
fibroblast_metadata <- metadata %>% dplyr::filter(cell_type_all4 %in% c('mCAF', 'WNT5A_BMP', 'WNT5A_infl', 'iCAF'))
tumor_fib_metadata <- metadata %>% dplyr::filter(cell_type_all4 %in% c('mCAF', 'WNT5A_BMP', 'WNT5A_infl', 'iCAF') | broad_cell_type == 'Tumor')

tumor_fib_count <- tumor_fib_metadata %>%
  group_by(orig.ident) %>%
  summarise(
    all_cell_count = n(),
    tumor_count = sum(broad_cell_type == "Tumor", na.rm = TRUE),
    # UPDATED: Add 'iCAF' to the total CAF count sum
    caf_count = sum(cell_type_all2 %in% c("mCAF", "WNT5A_BMP", "WNT5A_infl", "iCAF"), na.rm = TRUE), 
    non_can_tumor_count = sum(cell_type_all2 == "Non_Canonical_CRC_1", na.rm = TRUE),
    APCDD1_tumor_count = sum(cell_type_all2 == "APCDD1_CRC", na.rm = TRUE),
    can_stem_tumor_count = sum(cell_type_all2 == "Canonical_CRC_Stem", na.rm = TRUE),
    can_stem_prolif_count = sum(cell_type_all2 == "Canonical_CRC_Stem_Proliferation", na.rm = TRUE),
    can_intestine_tumor_count = sum(cell_type_all2 == "Canonical_CRC_Intestine", na.rm = TRUE),
    can_intestine_prolif_count = sum(cell_type_all2 == "Canonical_CRC_Intestine_Proliferation", na.rm = TRUE),
    mCAF_count = sum(cell_type_all2 == "mCAF", na.rm = TRUE),
    WNT5A_BMP_count = sum(cell_type_all2 == "WNT5A_BMP", na.rm = TRUE),
    WNT5A_infl_count = sum(cell_type_all2 == "WNT5A_infl", na.rm = TRUE),
    # NEW: Add iCAF count
    iCAF_count = sum(cell_type_all2 == "iCAF", na.rm = TRUE), 
    Organ = first(Site_of_Origin),
    Patient_ID = first(Patient_ID),
    Primary_Side = first(Primary_Side),
    Tx_in_6mo = first(Tx_in_6mo),
    .groups = "drop"
  ) %>%
  filter(tumor_count > 20, caf_count > 10) %>%
  mutate(orig.ident = str_replace(orig.ident, "^HT413C1-Th1K[24]A2Nd1_2Bma1_1$", "HT413C1-Th1"))

  tumor_caf_prop <- tumor_fib_count %>%
  group_by(orig.ident) %>%
  summarise(
    across(c(all_cell_count, tumor_count, caf_count,
             non_can_tumor_count, APCDD1_tumor_count,
             can_stem_tumor_count, can_stem_prolif_count,
             can_intestine_tumor_count, can_intestine_prolif_count,
             mCAF_count, WNT5A_BMP_count, WNT5A_infl_count,
             iCAF_count), # UPDATED: ADD iCAF_count HERE
           \(x) sum(x, na.rm = TRUE)),
    tumor_prop = round(100 * tumor_count / all_cell_count, 2),
    caf_prop = round(100 * caf_count / all_cell_count, 2),
    non_can_tumor_prop = round(100 * non_can_tumor_count / tumor_count, 2),
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
    iCAF_prop = round(100 * iCAF_count / caf_count, 2), # NEW: iCAF prop / CAF total
    non_can_tumor_all_prop = round(100 * non_can_tumor_count / all_cell_count, 2),
    can_stem_tumor_all_prop = round(100 * can_stem_tumor_count / all_cell_count, 2),
    can_intestine_tumor_all_prop = round(100 * can_intestine_tumor_count / all_cell_count, 2),
    can_intestine_prolif_all_prop = round(100 * can_intestine_prolif_count / all_cell_count, 2),
    mCAF_all_prop = round(100 * mCAF_count / all_cell_count, 2),
    WNT5A_BMP_all_prop = round(100 * WNT5A_BMP_count / all_cell_count, 2),
    WNT5A_infl_all_prop = round(100 * WNT5A_infl_count / all_cell_count, 2),
    iCAF_all_prop = round(100 * iCAF_count / all_cell_count, 2), # NEW: iCAF prop / ALL total
    Organ = first(Organ),
    Patient_ID = first(Patient_ID),
    Primary_Side = first(Primary_Side),
    Tx_in_6mo = first(Tx_in_6mo),
    .groups = "drop"
  )

colon_summary = tumor_caf_prop %>% dplyr::filter(Organ == 'colon')
liver_summary = tumor_caf_prop %>% dplyr::filter(Organ == 'liver')
lung_summary = tumor_caf_prop %>% dplyr::filter(Organ == 'lung')

organ_list <- list(
  colon = colon_summary,
  liver = liver_summary,
  lung = lung_summary
)

x_axis = c('mCAF_prop', 'WNT5A_BMP_prop', 'WNT5A_infl_prop', 'iCAF_prop')
y_axis = c('non_can_tumor_prop', 'can_stem_tumor_prop', 'can_intestine_tumor_prop', 'can_intestine_prolif_prop')

# -----------------------------
# 1. Combine liver + lung as metastasis
# -----------------------------
colon_summary <- tumor_caf_prop %>% filter(Organ == "colon")
meta_summary  <- tumor_caf_prop %>% filter(Organ %in% c("liver", "lung"))

organ_list2 <- list(
  colon = colon_summary,
  metastasis = meta_summary
)

# -----------------------------
# 2. Variables of interest
# -----------------------------
x_axis = c("mCAF_prop", "WNT5A_BMP_prop", "WNT5A_infl_prop", "iCAF_prop")
y_axis = c("non_can_tumor_prop", "can_stem_tumor_prop",
           "can_intestine_tumor_prop", "can_intestine_prolif_prop")

# -----------------------------
# 3. Function: compute Spearman correlations
# -----------------------------
compute_corr_matrix <- function(df, x_vars, y_vars) {
  corr_mat <- matrix(NA, nrow = length(y_vars), ncol = length(x_vars),
                     dimnames = list(y_vars, x_vars))
  p_mat <- corr_mat

  for (y in y_vars) {
    for (x in x_vars) {
      tmp <- cor.test(df[[x]], df[[y]], method = "spearman")
      corr_mat[y, x] <- tmp$estimate
      p_mat[y, x]    <- tmp$p.value
    }
  }
  return(list(corr = corr_mat, p = p_mat))
}

# -----------------------------
# 4. Convert p-values into stars
# -----------------------------
p_to_stars <- function(p) {
  ifelse(p < 0.05, "**",
         ifelse(p < 0.1, "*", ""))
}

# -----------------------------
# 5. Color scale for correlations
# -----------------------------
col_fun = colorRamp2(
  c(-0.75, -0.5, 0, 0.5, 0.75),
  c("#2166ac", "#67a9cf", "white", "#f4a582", "#ca0020")
)

# -----------------------------
# 6. Generate heatmaps
# -----------------------------
for (organ_name in names(organ_list2)) {

  df <- organ_list2[[organ_name]]

  # Compute correlations
  results <- compute_corr_matrix(df, x_axis, y_axis)
  corr_mat <- results$corr
  p_mat    <- results$p
  star_mat <- matrix(
    p_to_stars(p_mat),
    nrow = nrow(p_mat),
    dimnames = dimnames(p_mat)
  )

  # Prepare heatmap
  ht <- Heatmap(
    corr_mat,
    name = paste0(organ_name, "_cor"),
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,

    rect_gp = gpar(col = "black", lwd = 1.2),   # <-- grid lines

    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(star_mat[i, j], x, y,
                gp = gpar(fontsize = 14, fontface = "bold"))
    },

    column_title = paste0(toupper(organ_name),
                          " – CAF vs Tumor Subtype Correlations"),
    heatmap_legend_param = list(title = "Spearman Rho")
  )

  pdf(file.path(out_dir, paste0("corr_heatmap_", organ_name, "_combined.pdf")),
      width = 6, height = 5)
  draw(ht)
  dev.off()
}
