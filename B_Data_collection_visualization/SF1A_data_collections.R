###############################################################################
#  mCRC Sample Collection Combined Heatmap (hp + hlim + hlum + hom + hn)
#  Author: Evan Peng
#  Description: Reads Google Sheet metadata and produces annotation-only heatmap
#               with pch annotations preserved for selected assays
###############################################################################

# --- 1. Load Libraries -------------------------------------------------------
suppressPackageStartupMessages({
  library(googlesheets4)
  library(tidyverse)
  library(ComplexHeatmap)
})

gs4_deauth()  # read-only access

# --- 2. Read Google Sheets ---------------------------------------------------
sheet_url <- "https://docs.google.com/spreadsheets/d/1d4Bn3guQTzqf4tWx5fZFc772JC2GhpKVvsg70V7EPbQ/edit?gid=83959954#gid=83959954"

sample_tbl <- read_sheet(sheet_url, sheet = "Sample_collections") %>%
  select(-c(MSI, APC_WXS, TP53_WXS, KRAS_WXS, Note))

clinical_tbl <- read_sheet(sheet_url, sheet = "Clinical") %>%
  mutate(
    Age = if_else(Age_at_diagnosis > 50, ">50", "<50"),
    CRC_site = case_when(
      Site_in_Colon %in% c("Descending", "Sigmoid", "Left") ~ "Left",
      Site_in_Colon %in% c("Ascending", "Cecum", "Transverse", "Right") ~ "Right",
      Site_in_Colon %in% c("Both") ~ "Left/Right",
      TRUE ~ Site_in_Colon
    )
  ) %>%
  select(-c(SP_ID, Age_at_diagnosis, Site_in_Colon,
            MSI, APC, TP53, KRAS, Survival, FU_days, Note)) %>%
  rename(Prior_Tx = Neoadjuvant_therapy)

sample_tbl <- left_join(sample_tbl, clinical_tbl, by = "Case_ID")

# --- 3. Subset by Tissue Type -----------------------------------------------
sample_primary_tbl <- subset(sample_tbl, Tissue == "Primary")
sample_lim_tbl     <- subset(sample_tbl, Tissue == "Metastasis" & Organ == "Liver")
sample_lum_tbl     <- subset(sample_tbl, Tissue == "Metastasis" & Organ == "Lung")
sample_om_tbl      <- subset(sample_tbl, Tissue == "Metastasis" & !(Organ %in% c("Liver", "Lung")))
sample_normal_tbl  <- subset(sample_tbl, Tissue == "Normal")

data_columns <- c("Bulk_DNA", "Bulk_RNA", "snRNA", "snATAC", "Visium",
                  "Xenium_regular", "Xenium_5K", "Phenocycler", "Visium_HD")

# --- 4. Define Color Palettes ------------------------------------------------
tissue_col <- c(Primary = "gold3", Metastasis = "maroon3", Normal = "lightblue1")
organ_col <- c(Colon = "#C2B280", Rectum = "#604E97", Liver = "brown", Lung = "steelblue1",
               Brain = "bisque1", Adrenal = "cyan4", Breast = "violet", Spleen = "lightcoral",
               LN = "olivedrab", Ovary = "tan4", Pancreas = "indianred1")
age_col <- c(">50" = "lightblue3", "<50" = "lightblue1", "NA" = "gray75")
sex_col <- c("M" = "royalblue", "F" = "plum", "NA" = "gray75")
race_col <- c("White" = "darkgreen", "African American" = "firebrick2", "Unknown" = "gray75", "NA" = "gray75")
site_col <- c("Left" = "darkorange", "Right" = "cornflowerblue", "Left/Right" = "olivedrab3",
              "Rectum" = "red3", "Unknown" = "gray75", "NA" = "gray75")
tx_col <- c("Yes" = "salmon", "No" = "antiquewhite", "Unknown" = "gray75", "NA" = "gray75")
assay_col <- c(available = "palegreen", unavailable = "white", tumor_only = "palegreen")

# --- 5. Create Legends -------------------------------------------------------
make_lgd <- function(labels, colors, title, nrow = 1) {
  Legend(labels = labels, legend_gp = gpar(fill = colors), nrow = nrow,
         title = title, title_position = "topleft")
}
lgd_list <- list(
  make_lgd(names(tissue_col), tissue_col, "Tissue"),
  make_lgd(names(organ_col), organ_col, "Organ", nrow = 2),
  make_lgd(names(age_col), age_col, "Age"),
  make_lgd(names(sex_col), sex_col, "Sex"),
  make_lgd(names(race_col), race_col, "Race", nrow = 2),
  make_lgd(names(site_col), site_col, "CRC Site", nrow = 2),
  make_lgd(names(tx_col), tx_col, "Prior Tx", nrow = 2),
  make_lgd(names(assay_col), assay_col, "Availability")
)

# --- 6. Helper: Clean Assay Columns + Generate *_pch -------------------------
clean_assays_with_pch <- function(tbl, assay_cols) {
  if (nrow(tbl) == 0) return(tbl)
  tbl <- tbl %>%
    mutate(across(all_of(assay_cols), \(x) {
      if (is.list(x))
        vapply(x, \(y) if (length(y) == 0) NA_character_ else as.character(y[1]), character(1))
      else as.character(x)
    }))
  for (a in assay_cols) {
    pch_name <- paste0(a, "_pch")
    tbl <- tbl %>%
      mutate(
        !!pch_name := case_when(.data[[a]] %in% c("available", "unavailable") ~ NA_character_, TRUE ~ .data[[a]]),
        !!a := if_else(.data[[a]] != "unavailable", "available", .data[[a]])
      )
  }
  tbl
}

# --- 7. Create a unified Heatmap generator (with optional pch) ---------------
make_heatmap <- function(tbl, prefix, has_pch = FALSE) {
  if (nrow(tbl) == 0) return(NULL)
  tbl <- tbl[!apply(tbl[, data_columns] == "unavailable", 1, all), ]

  assay_cols <- c("Bulk_RNA", "snRNA", "snATAC", "Xenium_regular", "Visium", "Phenocycler")
  if (has_pch) tbl <- clean_assays_with_pch(tbl, assay_cols)

  columns <- c("Tissue", "Organ", "Age", "Sex", "Race", "CRC_site", "Prior_Tx",
               "Bulk_DNA", "Bulk_RNA", "snRNA", "snATAC",
               "Visium", "Xenium_regular", "Xenium_5K", "Phenocycler", "Visium_HD",
               if (has_pch) paste0(assay_cols, "_pch"))

  for (col in columns)
    assign(paste0(prefix, "_", col), setNames(as.character(tbl[[col]]), tbl$Tissue_ID), envir = .GlobalEnv)

  anno_args <- list(
    Tissue = anno_simple(get(paste0(prefix, "_Tissue")), col = tissue_col),
    Organ = anno_simple(get(paste0(prefix, "_Organ")), col = organ_col, border = TRUE),
    Age = anno_simple(get(paste0(prefix, "_Age")), col = age_col, border = TRUE),
    Sex = anno_simple(get(paste0(prefix, "_Sex")), col = sex_col, border = TRUE),
    Race = anno_simple(get(paste0(prefix, "_Race")), col = race_col, border = TRUE),
    CRC_site = anno_simple(get(paste0(prefix, "_CRC_site")), col = site_col, border = TRUE),
    Prior_Tx = anno_simple(get(paste0(prefix, "_Prior_Tx")), col = tx_col, border = TRUE)
  )

  # Assays (with pch if available)
  for (a in c("Bulk_DNA", "Bulk_RNA", "snRNA", "snATAC", "Visium",
              "Visium_HD", "Xenium_regular", "Xenium_5K", "Phenocycler")) {
    assay_val <- get(paste0(prefix, "_", a))
    pch_col <- paste0(prefix, "_", a, "_pch")
    if (has_pch && exists(pch_col)) {
      pch_val <- get(pch_col)
      anno_args[[a]] <- anno_simple(assay_val, col = assay_col, gp = gpar(col = "black"), pch = pch_val)
    } else {
      anno_args[[a]] <- anno_simple(assay_val, col = assay_col, gp = gpar(col = "black"))
    }
  }

  ha <- do.call(HeatmapAnnotation, anno_args)
  dummy <- matrix(0, nrow = 0, ncol = length(tbl$Tissue_ID))
  colnames(dummy) <- tbl$Tissue_ID

  Heatmap(dummy, name = prefix, top_annotation = ha,
          show_row_names = FALSE, show_column_names = TRUE,
          cluster_rows = FALSE, cluster_columns = FALSE)
}

# --- 8. Build heatmaps for each category ------------------------------------
hp   <- make_heatmap(sample_primary_tbl, "p")
hlim <- make_heatmap(sample_lim_tbl, "lim", has_pch = TRUE)
hlum <- make_heatmap(sample_lum_tbl, "lum")
hom  <- make_heatmap(sample_om_tbl, "om", has_pch = TRUE)
hn   <- make_heatmap(sample_normal_tbl, "n")

# --- 9. Combine and Draw -----------------------------------------------------
out_dir <- getwd()
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

pdf(file.path(out_dir, "mCRC_sample_collection.pdf"), width = 20, height = 8)
ha <- hp + hlim + hlum + hom + hn
draw(ha,
     padding = unit(c(5, 5, 5, 5), "mm"),
     annotation_legend_side = "bottom",
     annotation_legend_list = lgd_list)
dev.off()

message("✅ Combined mCRC sample collection heatmap saved successfully (with pch annotations).")


