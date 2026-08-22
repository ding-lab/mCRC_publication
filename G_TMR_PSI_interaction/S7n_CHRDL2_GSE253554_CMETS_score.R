#!/usr/bin/env Rscript

# Extended Data Fig. 7m / 8m: CMETS score after CHRDL2-mediated BMP-SMAD
# attenuation in Caco-2 cells (GSE253554). Revision #2 plot uses day on the
# x-axis and doxycycline dose as color.
# Source: Revesion/CMETS/CMETS_interaction_functional.ipynb

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(edgeR)
  library(GSVA)
})

output_dir <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/G_TMR_PSI_interaction"
)
base_dir <- Sys.getenv(
  "MCRC_GSE253554_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/CMETS/Functional"
)
geneset_path <- Sys.getenv(
  "CMETS_GENESET_RDS",
  unset = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets/WASHU_snRNA_CRC_genesets_updated.rds"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

geneset <- readRDS(geneset_path)
geneset_use <- list(
  Canonical_CRC_Stem = geneset$Canonical_CRC_Stem,
  CMETS = geneset$CMETS
)

sample_info <- tibble(
  sample_id = c(
    "GSM8023145", "GSM8023146", "GSM8023147", "GSM8023148",
    "GSM8023149", "GSM8023150", "GSM8023151", "GSM8023152",
    "GSM8023153", "GSM8023154", "GSM8023155", "GSM8023156"
  ),
  sample_name = c(
    "ChrDMSO13", "ChrLOW3", "ChrDMSO15", "ChrLOW5",
    "ChrDMSO17", "ChrLOW7", "ChrMED3", "ChrMED5",
    "ChrMED7", "ChrHIGH3", "ChrHIGH5", "ChrHIGH7"
  ),
  treatment = factor(
    c("DMSO", "Dox_low", "DMSO", "Dox_low", "DMSO", "Dox_low",
      "Dox_medium", "Dox_medium", "Dox_medium", "Dox_high", "Dox_high", "Dox_high"),
    levels = c("DMSO", "Dox_low", "Dox_medium", "Dox_high")
  ),
  day = factor(c(3, 3, 5, 5, 7, 7, 3, 5, 7, 3, 5, 7), levels = c(3, 5, 7))
)

cts <- fread(file.path(base_dir, "GSE253554_raw_counts_GRCh38.p13_NCBI.tsv"))
ann <- fread(file.path(base_dir, "Human.GRCh38.p13.annot.tsv"))
present_samples <- intersect(sample_info$sample_id, colnames(cts))
sample_info2 <- sample_info %>% filter(sample_id %in% present_samples)

cts2 <- cts %>%
  left_join(ann %>% select(GeneID, Symbol), by = "GeneID") %>%
  filter(!is.na(Symbol), Symbol != "") %>%
  mutate(total_count = rowSums(across(all_of(present_samples)))) %>%
  arrange(desc(total_count)) %>%
  distinct(Symbol, .keep_all = TRUE)

mat_raw <- as.data.frame(cts2 %>% select(Symbol, all_of(present_samples)))
rownames(mat_raw) <- mat_raw$Symbol
mat_raw <- as.matrix(mat_raw[, -1, drop = FALSE])[, sample_info2$sample_id, drop = FALSE]

dge <- DGEList(counts = mat_raw)
dge <- calcNormFactors(dge)
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
expr_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
rownames(expr_logcpm) <- toupper(rownames(expr_logcpm))

geneset_use2 <- lapply(geneset_use, function(gs) intersect(unique(toupper(gs)), rownames(expr_logcpm)))
if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  ssgsea_scores <- GSVA::gsva(GSVA::ssgseaParam(exprData = expr_logcpm, geneSets = geneset_use2), verbose = FALSE)
} else {
  ssgsea_scores <- GSVA::gsva(expr = expr_logcpm, gset.idx.list = geneset_use2, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
}

score_long <- as.data.frame(t(ssgsea_scores)) %>%
  rownames_to_column("sample_id") %>%
  left_join(sample_info2, by = "sample_id") %>%
  pivot_longer(cols = names(geneset_use2), names_to = "Geneset", values_to = "Score")

p <- ggplot(score_long, aes(x = day, y = Score, color = treatment, group = treatment)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Geneset, scales = "free_y") +
  theme_classic() +
  labs(x = "Day", y = "ssGSEA score", color = "Treatment")

ggsave(file.path(output_dir, "CHRDL2_GSE253554_CMETS_score_by_day.pdf"), p, width = 8, height = 4)
write.csv(score_long, file.path(output_dir, "CHRDL2_GSE253554_CMETS_scores.csv"), row.names = FALSE)
message("Wrote CHRDL2 / GSE253554 scores to ", output_dir)
