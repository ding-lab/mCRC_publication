# ---- Libraries & setup -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)      # includes dplyr, ggplot2, readr, stringr, purrr, tibble
  library(ComplexHeatmap)
  library(circlize)
  library(maftools)
  library(googlesheets4)
  library(data.table)
  library(survtype)
})
gs4_deauth()

source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')
# NOTE: assumes jupyter_support_functions.R defines KellyPalette

out_dir   <- getwd()
input_dir <- '/diskmnt/Projects/MetNet_analysis/Colorectal/data/maf_objects'
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- I/O ---------------------------------------------------------------------
merge.maf <- readRDS(file.path(input_dir, 'mCRC_pecgs_N100_maf.rds'))

# Prefer maftools' onco matrix creation (handles multi-hits and empty string nicely)
# Fallback to your maf2matrix path if needed.
if ("createOncoMatrix" %in% getNamespaceExports("maftools")) {
  onco <- maftools::createOncoMatrix(merge.maf)
  merge.mat <- as.matrix(onco$oncoMatrix)
} else {
  merge.mat <- as.matrix(maf2matrix(merge.maf@data))  # your original path
}

# No P1 somatic mutations were found in PE0068C1:
if (!"PE0068C1" %in% colnames(merge.mat)) {
  merge.mat <- cbind(merge.mat, PE0068C1 = "") 
}

# Cancer gene list
cancer_genes <- read.delim(
  "/diskmnt/Projects/PECGS_analysis/epeng/WXS/cancerGeneList.tsv",
  sep = "\t", header = TRUE, check.names = FALSE
) %>%
  filter(COSMIC_CGC_v99 == "Yes") %>%
  pull(Hugo_Symbol)

# ---- Helpers -----------------------------------------------------------------
str_count_semicolons <- function(x) stringr::str_count(x, fixed(";"))

# Reclassify multi-hit, then remove semicolons
normalize_calls <- function(mat) {
  # Ensure character matrix
  storage.mode(mat) <- "character"
  mat[] <- ifelse(str_count_semicolons(mat) >= 2, "Multi_hit", mat)
  mat[] <- gsub(";", "", mat, fixed = TRUE)
  mat
}

presence_binary <- function(mat) {
  # 1 if has any mutation string, else 0
  (mat != "" & !is.na(mat)) * 1L
}

pick_genes_by_freq <- function(mat, min_prop = 0.10) {
  bin <- presence_binary(mat)
  freq <- rowSums(bin, na.rm = TRUE) / ncol(bin)
  names(which(freq >= min_prop))
}

clean_sample_ids <- function(x) {
  x %>%
    # drop trailing run suffixes
    str_replace("(_T|Y1_T|Y2_T|D1_1-13_T|D1_1-12_T)$", "") %>%
    # collapse PE sample IDs like "PE0068C1-..." -> "PE0068C1"
    str_replace("^(PE[^-]+)-.*", "\\1") %>%
    # specific fixes from your script
    str_replace("^HT266C1-Th4$", "HT266C1-Tb1") %>%
    str_replace("^HT225C1-S2$",  "HT225C1-Th1") %>%
    str_replace("^HT413C1-Th1K4$", "HT413C1-Th1") %>%
    str_replace("^CM1519C1-S1$", "HT225C1-S1")
}

# Build an oncoPrint Heatmap for a selected set of samples
build_oncoprint <- function(mat, samples, annot_df, row_title, mutation_colors, heatmap_legend_param) {
  # Intersect to avoid subsetting errors
  samples_use <- intersect(samples, colnames(mat))
  if (length(samples_use) == 0) return(NULL)

  submat <- mat[, samples_use, drop = FALSE]

  # Annotation vectors in sample order
  ann_df <- annot_df %>%
    as.data.frame() %>%                    # ensure it's a plain data frame
    mutate(across(everything(), as.character)) %>%  # convert Rle/factor to character
    filter(Sample_ID %in% samples_use) %>%
    distinct(Sample_ID, .keep_all = TRUE)

  # Ensure order matches sample order
  ann_df <- ann_df[match(samples_use, ann_df$Sample_ID), , drop = FALSE]

  # Top annotation
  ha <- HeatmapAnnotation(
    "Organ"        = ann_df$Organ,
    "Primary site" = ann_df$Primary_Side,
    col = list(
      "Organ"        = organ_col,
      "Primary site" = site_col
    ),
    show_legend = TRUE
  )

  # Alter function (only include mutation types we actually have)
  muts_present <- sort(unique(na.omit(as.vector(submat[submat != ""]))))
  muts_present <- muts_present[muts_present %in% names(mutation_colors)]
  if (length(muts_present) == 0) return(NULL)

  alter_fun <- c(
    list(background = alter_graphic("rect", fill = "#CCCCCC")),
    setNames(
      lapply(muts_present, \(mt) alter_graphic("rect", fill = mutation_colors[[mt]])),
      muts_present
    )
  )

  oncoPrint(
    submat,
    alter_fun = alter_fun,
    col = mutation_colors,
    alter_fun_is_vectorized = FALSE,
    row_title = row_title,
    row_names_gp = gpar(fontsize = 8),
    top_annotation = ha,
    pct_side = "right",
    row_names_side = "left",
    heatmap_legend_param = heatmap_legend_param,
    show_column_names = TRUE
  )
}

# ---- Normalize alteration calls ----------------------------------------------
merge.mat <- normalize_calls(merge.mat)

# Focus on cancer genes
gene_interested <- intersect(rownames(merge.mat), cancer_genes)
merge.mat1 <- merge.mat[gene_interested, , drop = FALSE]

# Frequency filter (>=10%)
genes_10p <- pick_genes_by_freq(merge.mat1, min_prop = 0.10)
message(sprintf(">10%% genes: N = %d", length(genes_10p)))
merge.mat1 <- merge.mat1[genes_10p, , drop = FALSE]

# Clean column names
colnames(merge.mat1) <- clean_sample_ids(colnames(merge.mat1))

# ---- Clinical table & annotations --------------------------------------------
sheet_url <- "https://docs.google.com/spreadsheets/d/1d4Bn3guQTzqf4tWx5fZFc772JC2GhpKVvsg70V7EPbQ/edit?gid=83959954#gid=83959954"

clinical_tbl <- read_sheet(sheet_url, sheet = 'WXS_clinical') %>%
  as.data.frame() %>%
  select(any_of(c(
    "Sample_ID","Case_ID","Age","Sex","Race","Organ","Tissue_Type",
    "Primary_Side","Dx_Met","collections","mutation_calling","Storage"
  ))) %>%
  mutate(
    source_stage = case_when(
      Tissue_Type == "primary"    & Dx_Met == "No"  ~ "CRC_no_Met",
      Tissue_Type == "primary"    & Dx_Met == "Yes" ~ "CRC_with_Met",
      Tissue_Type == "metastasis"                    ~ "mCRC",
      TRUE ~ NA_character_
    ),
    source_stage2 = case_when(
      source_stage == "mCRC" & Organ == "liver" ~ "mCRC_Liver",
      source_stage == "mCRC" & Organ == "lung"  ~ "mCRC_Lung",
      source_stage == "mCRC"                    ~ "mCRC_Other",
      TRUE                                      ~ source_stage
    )
  ) %>%
  # Keep only samples present in matrix
  filter(Sample_ID %in% colnames(merge.mat1))

# Annotation palettes (kept from your script; adjust as needed)
methx_col <- c(
  CRC_no_Met = 'gold3',
  CRC_with_Met = 'gold1',
  mCRC_Liver = 'brown',
  mCRC_Lung  = 'steelblue1',
  mCRC_Other = 'tan3'
)

organ_col <- c(
  colon = '#C2B280', rectum = '#604E97', liver = 'brown', lung = 'steelblue1',
  brain = 'bisque1', adrenal = 'cyan4', breast = 'violet', spleen = 'lightcoral',
  LN ='seagreen' , ovary = 'tan', pancreas = 'sandybrown', 'NA' = 'ivory4'
)

site_col <- c('L' = 'darkorange', 'R' = 'cornflowerblue', 'Rectum' = 'red3')

# ---- Mutation color map (KellyPalette expected from your sourced script) ------
KellyPalette <- unlist(KellyPalette)
mutation_colors <- c(
  Missense_Mutation = KellyPalette["blue"],
  Nonsense_Mutation = KellyPalette["red"],
  Frame_Shift_Del   = KellyPalette["orange"],
  Multi_hit         = KellyPalette["black"],
  Frame_Shift_Ins   = KellyPalette["yellow"],
  Splice_Site       = KellyPalette["purpleRed"],
  In_Frame_Ins      = KellyPalette["green"],
  In_Frame_Del      = KellyPalette["violet"],
  Nonstop_Mutation  = KellyPalette["grey"] %||% "#7f7f7f"  # fallback if absent
)
# Ensure names have no trailing palette suffixes
names(mutation_colors) <- sub("\\..*$", "", names(mutation_colors))

# Legend labels: only those present in your matrix
mutation_types_present <- sort(unique(na.omit(as.vector(merge.mat1[merge.mat1 != ""]))))
mutation_types_present <- mutation_types_present[mutation_types_present %in% names(mutation_colors)]
heatmap_legend_param1 <- list(
  title = "Somatic Alterations",
  at = mutation_types_present,
  labels = mutation_types_present
)

# ---- Group samples & build panels --------------------------------------------
# Sample groups
grp <- clinical_tbl %>%
  transmute(Sample_ID, source_stage2)

get_samples <- function(stage_label) {
  clinical_tbl %>% filter(source_stage2 == stage_label) %>% pull(Sample_ID)
}

# One annotation df to reuse
annot_df <- clinical_tbl %>%
  select(Sample_ID, Organ, Primary_Side, source_stage2)

p1 <- build_oncoprint(merge.mat1, get_samples("CRC_no_Met"),   annot_df,
                      "Somatic mutations (>=10% samples)", mutation_colors, heatmap_legend_param1)
p2 <- build_oncoprint(merge.mat1, get_samples("CRC_with_Met"), annot_df,
                      "Somatic mutations (>=10% samples)", mutation_colors, heatmap_legend_param1)
p3 <- build_oncoprint(merge.mat1, get_samples("mCRC_Liver"),   annot_df,
                      "Somatic mutations (>=10% samples)", mutation_colors, heatmap_legend_param1)
p4 <- build_oncoprint(merge.mat1, get_samples("mCRC_Lung"),    annot_df,
                      "Somatic mutations (>=10% samples)", mutation_colors, heatmap_legend_param1)
p5 <- build_oncoprint(merge.mat1, get_samples("mCRC_Other"),   annot_df,
                      "Somatic mutations (>=10% samples)", mutation_colors, heatmap_legend_param1)

# Drop NULL panels (in case any group is empty)
panels <- list(p1,p2,p3,p4,p5) %>% keep(~ !is.null(.x))

# ---- Draw & save -------------------------------------------------------------
pdf(file = file.path(out_dir, 'mCRC_N100_organ_oncoprint.pdf'), width = 30, height = 4)
if (length(panels) == 0) {
  grid::grid.text("No samples to display after filtering.", gp = grid::gpar(cex = 1.2))
} else {
  ht <- Reduce(`+`, panels)
  draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
}
dev.off()

message("Saved: ", file.path(out_dir, 'mCRC_N100_organ_oncoprint.pdf'))
