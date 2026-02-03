# ---- Libraries & setup -------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(maftools)
  library(googlesheets4)
})
gs4_deauth()

# Optional default theme if not loaded elsewhere
if (!exists("theme_mydefault")) {
  theme_mydefault <- function(base_size = 12) {
    theme_minimal(base_size = base_size) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
      )
  }
}

source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')

input_dir  <- '/diskmnt/Projects/MetNet_analysis/Colorectal/data/maf_objects'
output_dir <- getwd()
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Read mutation matrix ----------------------------------------------------
merge.mat1 <- read.csv(
  file = file.path(input_dir, 'mCRC_pecgs_N100_name_cleaned.maf.mat.csv'),
  row.names = 1,
  check.names = FALSE
)
merge.mat1 <- as.matrix(merge.mat1)

# Determine KRAS mutation presence
stopifnot("KRAS" %in% rownames(merge.mat1))
kras_status <- as.integer(merge.mat1["KRAS", ] != "" & !is.na(merge.mat1["KRAS", ]))
names(kras_status) <- colnames(merge.mat1)

# ---- Read sample + case clinical tables --------------------------------------
sheet_url <- "https://docs.google.com/spreadsheets/d/1d4Bn3guQTzqf4tWx5fZFc772JC2GhpKVvsg70V7EPbQ/edit?gid=83959954#gid=83959954"

# Sample-level WXS metadata
sample_tbl <- read_sheet(sheet_url, sheet = 'WXS_clinical') %>%
  as.data.frame()
sample_col <- names(sample_tbl)[grepl("sample", names(sample_tbl), ignore.case = TRUE)][1]
stopifnot(!is.na(sample_col))

sample_tbl <- sample_tbl %>%
  rename(Sample_ID = all_of(sample_col)) %>%
  select(any_of(c("Sample_ID", "Case_ID"))) %>%
  mutate(across(everything(), as.character)) %>%
  distinct()

# Case-level clinical table
case_tbl <- read_sheet(sheet_url, sheet = 'Clinical') %>%
  mutate(
    Age = if_else(Age_at_diagnosis > 50, '>50', '<50'),
    CRC_site = case_when(
      Site_in_Colon %in% c('Descending', 'Sigmoid', 'Left') ~ 'Left',
      Site_in_Colon %in% c('Ascending', 'Cecum', 'Transverse', 'Right') ~ 'Right',
      Site_in_Colon %in% c('Both') ~ 'Left/Right',
      TRUE ~ as.character(Site_in_Colon)
    )
  ) %>%
  select(-c(SP_ID, Age_at_diagnosis, Site_in_Colon,
            MSI, APC, TP53, KRAS, KRAS_detail, Survival, FU_days, Note)) %>%
  rename(Prior_Tx = Neoadjuvant_therapy) %>%
  mutate(
    Met_status = case_when(
      Dx_with_Metastasis == 'No' ~ 'No_Met',
      Dx_with_Metastasis == 'Yes' & Lung_Met == 'No' ~ 'Met_wo_Lung',
      Lung_Met == 'Yes' ~ 'Lung_Met'
    )
  ) %>%
  mutate(across(everything(), ~ if (is.list(.)) as.character(.) else .)) %>%
  as.data.frame()

# ---- Map sample → case and aggregate KRAS presence ---------------------------
sample_kras <- tibble(Sample_ID = names(kras_status), KRAS_mut = kras_status) %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  inner_join(sample_tbl, by = "Sample_ID") %>%
  distinct(Case_ID, Sample_ID, KRAS_mut)

case_kras <- sample_kras %>%
  group_by(Case_ID) %>%
  summarise(
    KRAS = as.integer(any(KRAS_mut == 1)),
    .groups = "drop"
  )

case_tbl <- left_join(case_tbl, case_kras, by = "Case_ID")

# ---- Define metastasis status (No_Met / Liver_Met_only / Lung_Met) ----------
case_tbl2 <- case_tbl %>%
  filter(!is.na(Dx_with_Metastasis), !is.na(Liver_Met), !is.na(Lung_Met)) %>%
  mutate(LiLu_Met_status = case_when(
    Dx_with_Metastasis == 'No'            ~ 'No_Met',
    Lung_Met == 'Yes'                     ~ 'Lung_Met',
    Liver_Met == 'Yes' & Lung_Met == 'No' ~ 'Liver_Met_only',
    TRUE                                  ~ NA_character_
  )) %>%
  filter(!is.na(LiLu_Met_status), !is.na(KRAS))

# ---- Build per-sample KRAS variant summary ----------------------------------
merge.maf <- readRDS(file.path(input_dir, 'mCRC_pecgs_N100_maf.rds'))

kras <- subsetMaf(
  maf = merge.maf,
  genes = "KRAS",
  includeSyn = FALSE,
  fields = c("Tumor_Sample_Barcode","HGVSp_Short","HGVSc",
             "Variant_Classification","Variant_Type")
)

dt <- as.data.table(kras@data)
dt[, AA := sub("^p\\.", "", HGVSp_Short)]
dt[Variant_Classification != "Missense_Mutation",
   AA := paste0(Variant_Classification, " (", sub("^c\\.", "", HGVSc), ")")]

per_sample_summary <- dt[, .(
  n_kras_mut = .N,
  KRAS_detail = paste(sort(unique(AA)), collapse = ";")
), by = Tumor_Sample_Barcode][order(Tumor_Sample_Barcode)]

# Derive Case_ID from sample prefix (keep trailing C)
per_sample_summary <- per_sample_summary %>%
  mutate(Case_ID = str_extract(Tumor_Sample_Barcode, "^[A-Za-z0-9]+C"))

# Collapse duplicate Case_IDs
case_detail <- per_sample_summary %>%
  group_by(Case_ID) %>%
  summarise(
    KRAS_detail = paste(sort(unique(KRAS_detail[!is.na(KRAS_detail) & KRAS_detail != ""])), collapse = ";"),
    .groups = "drop"
  )

# # ---- Add manual case HT347C (No WXS but Medical record showed KRAS G12D)-----
# case_detail <- case_detail %>%
#   bind_rows(tibble(Case_ID = "HT347C", KRAS_detail = "G12D")) %>%
#   distinct(Case_ID, .keep_all = TRUE)

# ---- Merge back to case-level table -----------------------------------------
case_tbl3 <- left_join(case_tbl2, case_detail, by = "Case_ID") %>%
  mutate(
    KRAS_details = case_when(
      grepl("G12D", KRAS_detail) ~ "G12D",
      grepl("G12V", KRAS_detail) ~ "G12V",
      grepl("G13D", KRAS_detail) ~ "G13D",
      !is.na(KRAS_detail)        ~ "Other",
      TRUE                       ~ NA_character_
    )
  ) %>%
  mutate(KRAS_details = factor(KRAS_details, levels = c("G12D","G12V","G13D","Other")))

# ---- Calculate proportions ---------------------------------------------------
plot_df2 <- case_tbl3 %>%
  group_by(LiLu_Met_status) %>%
  summarise(
    total_cases = n(),
    kras_mut_cases = sum(KRAS == 1, na.rm = TRUE),
    prop_kras = kras_mut_cases / total_cases,
    .groups = "drop"
  ) %>%
  mutate(
    LiLu_Met_status = factor(LiLu_Met_status, levels = c("No_Met","Liver_Met_only","Lung_Met")),
    label = paste0(kras_mut_cases, "/", total_cases)
  )

denom_df <- case_tbl3 %>%
  group_by(LiLu_Met_status) %>%
  summarise(total_cases = n(), .groups = "drop")

stack_df <- case_tbl3 %>%
  filter(KRAS == 1, !is.na(KRAS_details), !is.na(LiLu_Met_status)) %>%
  group_by(LiLu_Met_status, KRAS_details) %>%
  summarise(kras_detail_cases = n(), .groups = "drop") %>%
  left_join(denom_df, by = "LiLu_Met_status") %>%
  mutate(
    LiLu_Met_status = factor(LiLu_Met_status, levels = c("No_Met","Liver_Met_only","Lung_Met")),
    prop = kras_detail_cases / total_cases
  )

# ---- Plot --------------------------------------------------------------------
kras_detail_cols <- c(G12D="#4277b6", G12V="#db6917", G13D="#5fa641", Other="#91218c")

p4 <- ggplot(stack_df, aes(x = LiLu_Met_status, y = prop, fill = KRAS_details)) +
  geom_col(width = 0.6) +
  geom_text(
    data = plot_df2,
    aes(x = LiLu_Met_status, y = prop_kras, label = label),
    vjust = -0.5, size = 3, inherit.aes = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = kras_detail_cols) +
  labs(x = "Metastasis Status",
       y = "KRAS-mutant Case Proportion",
       fill = "KRAS detail") +
  theme_mydefault() +
  theme(legend.position = "right")


pdf(file.path(output_dir, "mCRC_LiLu_KRAS_mut_detail_frequency_barplot.pdf"), width = 4, height = 4)
p4
dev.off()
message("Saved: ", file.path(output_dir, "mCRC_KRAS_mut_detail_frequency_barplot.pdf"))

