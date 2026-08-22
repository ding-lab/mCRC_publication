#!/usr/bin/env Rscript

# Extended Data Fig. 7m-n: Tabula Sapiens organ- and cell-type-level expression
# of HGF, HBEGF, and SLIT2.
# Source: Revesion/Organ_specific_interactions/TSC_HBEGF_HGF_validation.ipynb
# Input CSVs are CELLxGENE gene-expression exports.

suppressPackageStartupMessages({
  library(tidyverse)
})

output_dir <- Sys.getenv(
  "MCRC_OUTPUT_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/H_Organ_specific_adaption"
)
input_dir <- Sys.getenv(
  "MCRC_TABULA_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/HGF_HBEGF"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(
  file.path(input_dir, "CELLxGENE_gene_expression_HGF_SLIT2_HBEGF_Tabula_Sapiens.csv"),
  comment.char = "#",
  check.names = FALSE
)

gene_levels <- c("HGF", "HBEGF", "SLIT2")
tissue_levels <- c("colon", "liver", "lung")

z_by_gene <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(x)) else as.numeric(scale(x))
}

plot_organs <- df %>%
  filter(
    Tissue %in% tissue_levels,
    `Cell Type` == "aggregated",
    `Gene Symbol` %in% gene_levels
  ) %>%
  mutate(Percent_Expressed = (`Number of Cells Expressing Genes` / `Cell Count`) * 100) %>%
  group_by(`Gene Symbol`) %>%
  mutate(Z_by_Gene = z_by_gene(Expression)) %>%
  ungroup() %>%
  mutate(
    Tissue = factor(Tissue, levels = tissue_levels, labels = c("Colon", "Liver", "Lung")),
    `Gene Symbol` = factor(`Gene Symbol`, levels = rev(gene_levels))
  )

p_organs <- ggplot(plot_organs, aes(x = Tissue, y = `Gene Symbol`)) +
  geom_point(aes(size = Percent_Expressed, color = Z_by_Gene)) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Z-score") +
  scale_size_continuous(name = "Percent expressed", range = c(3, 12)) +
  theme_bw() +
  labs(x = "Organ", y = "Gene")

plot_celltypes <- df %>%
  filter(
    Tissue %in% tissue_levels,
    `Cell Type` != "aggregated",
    `Gene Symbol` %in% gene_levels
  ) %>%
  mutate(
    Percent_Expressed = (`Number of Cells Expressing Genes` / `Cell Count`) * 100,
    Cell_Group = paste0(str_to_title(Tissue), ": ", `Cell Type`)
  ) %>%
  group_by(`Gene Symbol`) %>%
  mutate(Z_by_Gene = z_by_gene(Expression)) %>%
  ungroup() %>%
  mutate(`Gene Symbol` = factor(`Gene Symbol`, levels = rev(gene_levels)))

p_cells <- ggplot(plot_celltypes, aes(x = Cell_Group, y = `Gene Symbol`)) +
  geom_point(aes(size = Percent_Expressed, color = Z_by_Gene)) +
  scale_color_gradient2(low = "blue", mid = "grey90", high = "red", midpoint = 0, name = "Z-score") +
  scale_size_continuous(name = "Percent expressed", range = c(2, 10)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8)) +
  labs(x = "Cell type", y = "Gene")

ggsave(file.path(output_dir, "S7m_Tabula_Sapiens_HGF_HBEGF_SLIT2_organs.pdf"), p_organs, width = 5, height = 4)
ggsave(file.path(output_dir, "S7m_Tabula_Sapiens_HGF_HBEGF_SLIT2_celltypes.pdf"), p_cells, width = 12, height = 4)
message("Wrote Tabula Sapiens ligand plots to ", output_dir)
