#!/usr/bin/env Rscript
# Single-sample spatial exemplar of local KRAS VAF at 100/200/300/400 um.
#
# S15-1909-A5U1 is the default: largest tumor cell count and best transcript
# coverage, so it is the cleanest illustration of a spatially uniform section.
#
# Writes a wide 1x4 strip (to sit as one row of a multi-panel figure) and a
# square 2x2 version.

suppressPackageStartupMessages({
  library(tidyverse)
})

vaf_root <- Sys.getenv(
  "MCRC_KRAS_VAF_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_local_VAF_output"
)
vaf_dir <- file.path(vaf_root, "Xenium_local_VAF")
input_csv <- file.path(vaf_dir, "Xenium_KRAS_local_VAF_downsampled_100-400um.csv")
output_dir <- file.path(vaf_dir, "spatial_exemplar")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

exemplar <- "S15-1909-A5U1"
radii <- c(100, 200, 300, 400)
pt_size <- 0.9
spatial_bg <- "white"
vaf_na_color <- "grey55"

df <- read_csv(input_csv, show_col_types = FALSE) %>%
  filter(Sample == exemplar, radius_um %in% radii) %>%
  mutate(radius_lab = factor(paste0(radius_um, " um"), levels = paste0(radii, " um"))) %>%
  arrange(desc(is.na(local_VAF)))

if (nrow(df) == 0) stop("No rows for ", exemplar)

kras <- unique(df$KRAS_call)
svaf <- unique(df$sample_VAF)[1]
n_low <- df %>%
  group_by(radius_um) %>%
  summarise(low = mean(is.na(local_VAF)), .groups = "drop")
message("Low-transcript fraction in plotted downsample:")
print(as.data.frame(n_low), row.names = FALSE, digits = 3)

build <- function(ncol_facet) {
  ggplot(df, aes(x = x, y = y, color = local_VAF)) +
    geom_point(size = pt_size, stroke = 0) +
    facet_wrap(~ radius_lab, ncol = ncol_facet) +
    scale_color_viridis_c(
      option = "viridis",
      begin = 0.15,
      end = 1,
      limits = c(0, 1),
      na.value = vaf_na_color,
      name = "Local\nKRAS VAF"
    ) +
    labs(
      title = paste0(exemplar, "  (", kras, ", section VAF = ", sprintf("%.2f", svaf), ")"),
      subtitle = "Each cell colored by the KRAS VAF of all tumor cells within the stated radius. Grey = <10 KRAS transcripts nearby."
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey30"),
      legend.position = "right",
      strip.text = element_text(size = 11, face = "bold", color = "black"),
      strip.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = spatial_bg, color = NA),
      panel.background = element_rect(fill = spatial_bg, color = NA),
      aspect.ratio = 1
    )
}

p_wide <- build(4)
ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_exemplar_", exemplar, "_1x4.pdf")),
       p_wide, width = 15, height = 4.6, bg = spatial_bg)
ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_exemplar_", exemplar, "_1x4.png")),
       p_wide, width = 15, height = 4.6, dpi = 300, bg = spatial_bg)

p_sq <- build(2)
ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_exemplar_", exemplar, "_2x2.pdf")),
       p_sq, width = 9.5, height = 8.8, bg = spatial_bg)
ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_exemplar_", exemplar, "_2x2.png")),
       p_sq, width = 9.5, height = 8.8, dpi = 300, bg = spatial_bg)

message("Wrote exemplar spatial panels to ", output_dir)
