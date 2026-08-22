# ======================
# Libraries
# ======================
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DESeq2)
library(CMScaller)
library(ggplot2)
library(ggalluvial)
library(ggpubr)

source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/mCRC_colors.R")

# ======================
# Output directory
# ======================
outdir = getwd()
if (!dir.exists(outdir)) dir.create(outdir)

# ======================
# Load featureCounts output
# ======================
base_dir <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/RNAseq/Featurecount"

files <- list.files(
  base_dir,
  pattern = "featurecounts_unstranded_readcount.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# read_fc <- function(path) {
#   df <- read_tsv(path, comment = "#")
#   sample_id <- basename(dirname(path))
#   df <- df %>% select(Geneid, !!sample_id := 7)
#   return(df)
# }

read_fc <- function(path) {
  df <- read_tsv(path, comment = "#", show_col_types = FALSE)

  sample_id <- basename(dirname(path))
  count_col <- colnames(df)[ncol(df)]

  # Extract the 2 required columns
  df2 <- df[, c("Geneid", count_col)]

  # Rename the count column to sample_id
  colnames(df2)[2] <- sample_id

  return(df2)
}

count_list <- lapply(files, read_fc)

count_mat <- purrr::reduce(count_list, full_join, by = "Geneid") %>%
  column_to_rownames("Geneid") %>%
  replace(is.na(.), 0)

# ======================
# Convert Ensembl → Symbol
# ======================
count_mat <- count_mat %>% rownames_to_column("ensembl_id")
count_mat$ensembl_id <- sub("\\..*", "", count_mat$ensembl_id)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = count_mat$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

count_mat$symbol <- gene_symbols
count_mat <- count_mat[!is.na(count_mat$symbol), ]

count_mat <- count_mat %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

count_mat <- as.data.frame(count_mat)
rownames(count_mat) <- count_mat$symbol
count_mat$symbol <- NULL

# ======================
# VST normalization via DESeq2
# ======================
col_data <- data.frame(
  condition = rep("unknown", ncol(count_mat)),
  row.names = colnames(count_mat)
)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_mat),
  colData = col_data,
  design = ~1
)

dds <- estimateSizeFactors(dds)
vst_mat <- vst(dds, blind = TRUE)
vst_expr <- assay(vst_mat)

# ======================
# Prepare matrix for CMS caller
# ======================
symbol2entrez <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vst_expr),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

keep <- !is.na(symbol2entrez)

emat <- vst_expr[keep, , drop = FALSE]
rownames(emat) <- as.character(symbol2entrez[rownames(emat)])

# Revision #1: CMScaller for primary CRC, lmCMScaller for liver metastases.
# Paired primary-to-metastasis transitions were dropped because only one
# matched pair remained callable after the liver-tailored classifier.
colnames(count_mat)[colnames(count_mat) == "CM1519C1-S1Y2E1"] <- "HT225C1-S1Y2E1"
colnames(count_mat)[colnames(count_mat) == "CM1519C1-T1Y2E1"] <- "HT225C1-T1Y2E1"
count_mat <- count_mat[, colnames(count_mat) != "HT225C1-S2T2A4E1", drop = FALSE]

colon_ids <- c(
  "HT225C1-S1Y2E1", "CM1563C1-S1Y1E1", "CM268C1-S1Y2E1", "CM318C1-S1Y1E1",
  "CM329C1", "CM349C1", "CM350C1-S1Y1E1", "CM357C1-S1Y1E1", "CM374C1-S1Y1E1",
  "CM376C1-S1Y1E1", "CM392C1-S1Y3E1", "CM397C1-S1Y1E1", "CM618C1-S1Y2E1",
  "CM655C1-S1E2", "CM663C1-S1Y2E1", "CM724C1-S1E1", "CM873C1-S1E1",
  "HT472C1-S1H4FC2E1"
)
rectum_ids <- c("CM492C1-S1Y2E1", "CM743C1-S1E1")
liver_ids <- c(
  "CM1563C1-T1Y1E1", "CM354C1-T3Y1E1", "CM392C2-Th1Y2E1",
  "CM426C1-Th1A2E1", "CM492C2-T1Y2E1", "CM492C2-T2Y2E1", "CM556C1-T1Y2E1",
  "CM618C2-T1Y2E1", "CM663C1-T1Y2E1", "HT112C1", "HT165C1-T1A3Y1E2",
  "HT179C1-T1A2K2E3", "HT186C1-T1A3Y1E1", "HT225C1-T1Y2E1", "HT230C1-T2A3Y1E1",
  "HT254C1-Th1K4A2E1", "HT260C1-Th1K1A3E1", "HT307C1-Th1K1A3E4",
  "HT413C1-Th1K4A1E1", "HT472C1-Th1K2E1", "HT525C1-Th1K2E1", "HT342C1-Th1K3A3E1"
)

primary_ids <- intersect(c(colon_ids, rectum_ids), colnames(count_mat))
liver_keep <- intersect(liver_ids, colnames(count_mat))

cms_primary <- CMScaller(
  emat = as.matrix(count_mat[, primary_ids, drop = FALSE]),
  RNAseq = TRUE,
  rowNames = "symbol",
  doPlot = FALSE
)
cms_liver <- lmCMScaller(
  emat = as.matrix(count_mat[, liver_keep, drop = FALSE]),
  RNAseq = TRUE,
  rowNames = "symbol"
)

extract_call <- function(res_obj, classifier_name) {
  out <- as.data.frame(res_obj)
  pred_col <- intersect(c("prediction", "class", "CMS", "subtype", "predictedCMS"), colnames(out))[1]
  tibble(
    sample_id = rownames(out),
    classifier = classifier_name,
    CMS = as.character(out[[pred_col]])
  )
}

cms_results <- bind_rows(
  extract_call(cms_primary, "CMScaller"),
  extract_call(cms_liver, "lmCMScaller")
) %>%
  column_to_rownames("sample_id")  

# ======================
# Annotate Organ & CMS
# ======================
cms_results <- as.data.frame(cms_results)
cms_results$CMS[is.na(cms_results$CMS) | cms_results$CMS %in% c("", "<NA>")] <- "NA"

cms_results <- cms_results %>%
  mutate(
    Organ = case_when(
      # Colon
      rownames(cms_results) %in% c(
        'HT225C1-S1Y2E1', 'CM1563C1-S1Y1E1', 'CM268C1-S1Y2E1', 'CM318C1-S1Y1E1',
        'CM329C1', 'CM349C1', 'CM350C1-S1Y1E1', 'CM357C1-S1Y1E1', 'CM374C1-S1Y1E1',
        'CM376C1-S1Y1E1', 'CM392C1-S1Y3E1', 'CM397C1-S1Y1E1', 'CM618C1-S1Y2E1',
        'CM655C1-S1E2', 'CM663C1-S1Y2E1', 'CM724C1-S1E1', 'CM873C1-S1E1',
        'HT472C1-S1H4FC2E1'
      ) ~ 'Colon',

      # Rectum
      rownames(cms_results) %in% c(
        'CM492C1-S1Y2E1', 'CM743C1-S1E1'
      ) ~ 'Rectum',

      # Liver
      rownames(cms_results) %in% c(
        'CM1563C1-T1Y1E1', 'CM354C1-T3Y1E1', 'CM392C2-Th1Y2E1',
        'CM426C1-Th1A2E1', 'CM492C2-T1Y2E1', 'CM492C2-T2Y2E1',
        'CM556C1-T1Y2E1', 'CM618C2-T1Y2E1', 'CM663C1-T1Y2E1',
        'HT112C1', 'HT165C1-T1A3Y1E2', 'HT179C1-T1A2K2E3', 'HT186C1-T1A3Y1E1',
        'HT225C1-S2T2A4E1', 'HT230C1-T2A3Y1E1', 'HT254C1-Th1K4A2E1',
        'HT260C1-Th1K1A3E1', 'HT307C1-Th1K1A3E4', 'HT413C1-Th1K4A1E1',
        'HT472C1-Th1K2E1', 'HT525C1-Th1K2E1', 'HT342C1-Th1K3A3E1', "HT225C1-T1Y2E1"
      ) ~ 'Liver',

      # Lung
      rownames(cms_results) %in% c(
        'CM354C2-T1Y2E1', 'CM426C2-Tp1E1', 'CM478C1-T1Y2E1', 'CM478C2-T1Y2E1'
      ) ~ 'Lung',

      # Other metastasis sites
      rownames(cms_results) %in% c('CM556C1-T2Y2E1', 'CM556C1-T3Y2E1') ~ 'Adrenal',
      rownames(cms_results) %in% c('CM556C2-T1Y2E1', 'CM556C3-T1Y2E1') ~ 'Brain',
      rownames(cms_results) %in% c('HT213C1-T1A3E1') ~ 'Spleen',
      rownames(cms_results) %in% c('HT226C1-T2A2E1') ~ 'CHOL',
      rownames(cms_results) %in% c('HT266C1-Th4K2E1') ~ 'Breast',

      TRUE ~ NA_character_
    ),

    case_id = str_extract(rownames(cms_results), "^[A-Za-z0-9]+(?=C[0-9])"),
    CMS = factor(CMS, levels = c("CMS1", "CMS2", "CMS3", "CMS4", "NA"))
  ) %>%
  filter(Organ != "CHOL") %>%
  mutate(
    Organ = factor(Organ),
    Organ2 = if_else(Organ %in% c('Colon', 'Rectum'), 'Colorectal', Organ),
    Organ2 = factor(Organ2, 
                    levels = c('Colorectal', 'Liver', 'Lung', 
                               'Breast', 'Adrenal', 'Brain', 'Spleen')),
    PM_status = if_else(Organ2 == "Colorectal", "Primary", "Metastasis")
  )

# ======================
# Save table
# ======================
write.csv(cms_results, file.path(outdir, "cms_results.csv"))

# ======================
# Barplot: all organs
# ======================
p1 <- ggplot(cms_results, aes(x = Organ2, fill = CMS)) +
  geom_bar(position = "fill", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c(
    CMS1="#1F77B4", CMS2="#FF7F0E",
    CMS3="#2CA02C", CMS4="#D62728", "NA"="grey70"
  )) +
  labs(
    title = "CMS subtype distribution by organ",
    x = "Organ", y = "Proportion of samples", fill = "CMS"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf(file.path(outdir, "CMS_barplot.pdf"), width=4, height=4)
p1
dev.off()

# Paired alluvium plot was removed in Revision #1 (only one callable pair
# remained after lmCMScaller). Keep colorectal vs liver comparison only.

# ======================
# Colon vs Liver barplot + Fisher tests
# ======================
cms_bar_df <- cms_results %>%
  filter(Organ2 %in% c("Colorectal", "Liver"), CMS != "NA")

cms_counts <- cms_bar_df %>%
  dplyr::count(Organ2, CMS) %>%
  group_by(Organ2) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

pvals <- cms_counts %>%
  group_by(CMS) %>%
  summarise(
    p_value = {
      col_n <- n[Organ2=="Colorectal"]
      liv_n <- n[Organ2=="Liver"]
      other_col <- sum(cms_counts$n[cms_counts$Organ2=="Colorectal" & cms_counts$CMS!=unique(CMS)])
      other_liv <- sum(cms_counts$n[cms_counts$Organ2=="Liver" & cms_counts$CMS!=unique(CMS)])
      
      tab <- matrix(c(col_n, other_col, liv_n, other_liv),
                    nrow=2, byrow=TRUE)
      
      if (any(is.na(tab)) || any(tab==0)) NA_real_ else fisher.test(tab)$p.value
    },
    .groups="drop"
  ) %>%
  mutate(
    p_label = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

cms_counts <- left_join(cms_counts, pvals, by="CMS")

p3 <- ggplot(cms_counts, aes(x = CMS, y = prop, fill = Organ2)) +
  geom_col(position=position_dodge(width=0.8), color="black", width=0.7) +
  geom_text(
    aes(label = paste0(scales::percent(prop, accuracy=1), " (", n, ")")),
    position = position_dodge(width=0.8),
    vjust = -0.3, size = 4
  ) +
  geom_text(
    data = cms_counts %>% distinct(CMS, p_label),
    aes(x = CMS, y = 0.48, label = p_label),
    inherit.aes = FALSE, size = 5
  ) +
  scale_fill_manual(values=c("Colorectal"="#C2B280", "Liver"="brown")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.5)) +
  labs(
    title = "CMS subtype distribution: Colorectal vs Liver",
    x = "CMS subtype", y = "Proportion of samples", fill = "Organ"
  ) +
  theme_bw(base_size = 14) +
  theme_mydefault()

pdf(file.path(outdir, "CMS_barplot_colon_vs_liver.pdf"), width=6, height=4)
p3
dev.off()
