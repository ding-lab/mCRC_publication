# Standalone: KRAS pathway / related target gene dotplot on tumor epithelial subset.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
})

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

MAPK_target_sig <- c("VEGFA", "MCL1", "FOS", "FOSB")
PI3K_AKT_mTOR_target_sig <- c("SLC2A1", "HK2", "LDHA", "XIAP", "NFKB1", "NFKB2", "HIF1A")
RalGEF_RalA_RalB_target_sig <- c("ACTB", "CDC42", "EXOC2")
Hippo_YAP_TAZ_target_sig <- c("CCN1", "CCN2", "AXL", "ITGB2", "ZEB1")

KRAS_target_sig <- c(
  MAPK_target_sig,
  PI3K_AKT_mTOR_target_sig,
  RalGEF_RalA_RalB_target_sig,
  Hippo_YAP_TAZ_target_sig
)

p6 <- DotPlot(tumor_obj2, features = KRAS_target_sig, group.by = "cell_type_all3") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  coord_fixed() +
  scale_color_gradientn(
    colors = rdylbu_colors,
    limits = c(-3, 3),
    breaks = color_breaks
  ) +
  scale_size_area(limits = c(0, 100), oob = scales::squish)

options(repr.plot.width = 12, repr.plot.height = 12)
p6

pdf(file.path(output_dir, "Dotplot_mCRC_CMETS_KRAS_target.pdf"), width = 10, height = 5)
print(p6)
dev.off()
