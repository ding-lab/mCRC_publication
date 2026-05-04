# Standalone: AP1 / NCC-related target gene dotplot on tumor epithelial subset.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

support_r <- "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/support_scripts/jupyter_support_functions.R"
if (file.exists(support_r)) {
  source(support_r)
}

output_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/Revision/Neuroendocrine_version"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tumor_rds_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4_tumor_subset_NE_included.rds"
tumor_obj <- readRDS(tumor_rds_path)

tumor_epi_types <- c(
  "CMETS",
  "Neuroendocrine-like tumor",
  "APCDD1_CRC",
  "Canonical_CRC_Stem",
  "Canonical_CRC_Stem_Proliferation",
  "Canonical_CRC_Intestine_Proliferation",
  "Canonical_CRC_Intestine"
)

tumor_obj2 <- subset(tumor_obj, subset = cell_type_all3 %in% tumor_epi_types)
tumor_obj2$cell_type_all3 <- factor(tumor_obj2$cell_type_all3, levels = rev(tumor_epi_types))
DefaultAssay(tumor_obj2) <- "SCT"

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

NCC_AP1_sorted_genes <- c(
  "SLC2A1", "LDHA", "PDK1", "HSPB1", "HSPA1B",
  "LAMC2", "LAMB3", "LAMA3", "ITGA3", "LTBP1", "SDC4",
  "TGFA", "JAG1", "MAPKAPK2", "BIRC2",
  "PLAUR", "IL13RA1", "CD55",
  "MCL1", "PPP1R13L", "EMP1"
)

p7 <- DotPlot(tumor_obj2, features = NCC_AP1_sorted_genes, group.by = "cell_type_all3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_fixed() +
  scale_color_gradientn(
    colors = rdylbu_colors,
    limits = c(-3, 3),
    breaks = color_breaks
  ) +
  scale_size_area(limits = c(0, 100), oob = scales::squish)

if (exists("set_size", mode = "function")) {
  set_size(6, 6)
}
p7

pdf(file.path(output_dir, "Dotplot_mCRC_CMETS_AP1_target.pdf"), width = 10, height = 5)
print(p7)
dev.off()
