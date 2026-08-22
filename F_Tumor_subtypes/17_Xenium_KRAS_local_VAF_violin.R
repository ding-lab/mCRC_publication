#!/usr/bin/env Rscript
# Reviewer-facing version of the local KRAS VAF violin figure.
#
# Differences from the internal version:
#   - the low-coverage label is spelled out instead of abbreviated "LD"
#   - radii where most cells lack 10 KRAS transcripts are drawn faded and
#     dashed, and called out as not interpretable rather than left to the reader
#   - the "% below 50% VAF" annotation is replaced by one sentence, since it is
#     near zero everywhere and only cluttered the panels
#
# Exact per-radius coverage comes from the neighborhood run (script 15)
# sample-summary and stats tables, which count every tumor cell rather than
# the plotting downsample.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

vaf_root <- Sys.getenv(
  "MCRC_KRAS_VAF_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_local_VAF_output"
)
vaf_dir <- file.path(vaf_root, "Xenium_local_VAF")
sample_summary_csv <- file.path(vaf_dir, "Xenium_KRAS_local_VAF_sample_summary.csv")
stats_csv <- file.path(vaf_dir, "Xenium_KRAS_local_VAF_stats.csv")
output_dir <- file.path(vaf_dir, "violin_reviewer")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_levels <- c(
  "CM1799-A1-Th1Fp1U1",
  "PE0319C1-A4-S1Fp1U1",
  "S15-1909-A5U1",
  "S15-1909-D2U1",
  "S16-38794-A4U2",
  "S16-38794-E2U1",
  "S17-16442-AFR1U1",
  "S20-46186-A16U1"
)

# Above this fraction of cells without 10 nearby KRAS transcripts, the violin
# describes a small unrepresentative subset and is not interpreted.
underpowered_cutoff <- 0.50

sample_n <- read_csv(sample_summary_csv, show_col_types = FALSE) %>%
  select(Sample, KRAS_call, n_tumor_cells)
coverage <- read_csv(stats_csv, show_col_types = FALSE) %>%
  left_join(sample_n, by = c("Sample", "KRAS_call")) %>%
  transmute(
    Sample,
    radius_um,
    KRAS_call,
    n_tumor_cells,
    n_informative = n_cells,
    low_count_fraction = 1 - n_cells / n_tumor_cells,
    underpowered = low_count_fraction > underpowered_cutoff
  )

build_figure <- function(csv_path, radii, tag, radius_palette) {
  vaf_df <- read_csv(csv_path, show_col_types = FALSE) %>%
    filter(Sample %in% sample_levels, radius_um %in% radii, !is.na(local_VAF)) %>%
    mutate(
      Sample = factor(Sample, levels = sample_levels),
      radius = factor(radius_um, levels = radii)
    )

  cov <- coverage %>%
    filter(radius_um %in% radii) %>%
    mutate(
      Sample = factor(Sample, levels = sample_levels),
      radius = factor(radius_um, levels = radii),
      label = percent(low_count_fraction, accuracy = 0.1)
    )

  max_below50 <- vaf_df %>%
    group_by(Sample, radius_um) %>%
    summarise(f = mean(local_VAF < 0.5), .groups = "drop") %>%
    pull(f) %>%
    max(na.rm = TRUE)

  violin_df <- vaf_df %>%
    filter(local_VAF >= 0.50, local_VAF <= 1) %>%
    left_join(cov %>% select(Sample, radius, underpowered), by = c("Sample", "radius"))

  vd_ok <- violin_df %>% filter(!underpowered)
  vd_low <- violin_df %>% filter(underpowered)

  p <- ggplot(mapping = aes(x = radius, y = local_VAF, fill = radius)) +
    geom_violin(
      data = vd_ok,
      scale = "width", trim = TRUE, linewidth = 0.3, color = "grey25"
    ) +
    geom_violin(
      data = vd_low,
      scale = "width", trim = TRUE, linewidth = 0.3,
      color = "grey55", linetype = "22", alpha = 0.30
    ) +
    geom_boxplot(
      data = violin_df,
      width = 0.12, outlier.shape = NA, fill = "white",
      alpha = 0.85, linewidth = 0.3
    ) +
    geom_text(
      data = cov,
      aes(x = radius, y = 1.035, label = label, fontface = ifelse(underpowered, "italic", "plain")),
      inherit.aes = FALSE,
      size = 2.6,
      color = "grey15"
    ) +
    facet_wrap(~ Sample, ncol = 4) +
    scale_x_discrete(limits = as.character(radii)) +
    scale_fill_manual(values = radius_palette, guide = "none") +
    scale_y_continuous(
      breaks = seq(0.5, 1, 0.1),
      labels = percent_format(accuracy = 1),
      limits = c(0.5, 1.06),
      expand = expansion(mult = c(0.01, 0))
    ) +
    labs(
      title = "Local KRAS VAF is stable across spatial scales within each section",
      subtitle = paste0(
        "Each violin pools KRAS mutant + wild-type transcripts from all tumor cells within the stated radius of each cell ",
        "(neighborhoods with >=10 transcripts); white box = median and IQR.\n",
        "Number above each violin = % of tumor cells with <10 KRAS transcripts nearby, i.e. below the detection limit of the mutation probes.\n",
        "Faded dashed violins (italic %): >", percent(underpowered_cutoff, accuracy = 1),
        " of cells below that limit, reflecting probe sensitivity rather than absent signal; shown for completeness, not interpreted.\n",
        "Across all samples and radii, fewer than ", percent(max_below50, accuracy = 0.1),
        " of measurable neighborhoods fell below 50% VAF."
      ),
      x = "Neighborhood radius (um)",
      y = "Local KRAS VAF"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8.2, lineheight = 1.15),
      axis.text.x = element_text(face = "bold")
    )

  ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_violin_reviewer_", tag, ".pdf")),
         p, width = 14, height = 7.8)
  ggsave(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_violin_reviewer_", tag, ".png")),
         p, width = 14, height = 7.8, dpi = 300)

  cov %>%
    select(Sample, radius_um, n_tumor_cells, n_informative, low_count_fraction, underpowered) %>%
    arrange(Sample, radius_um) %>%
    write_csv(file.path(output_dir, paste0("Xenium_KRAS_local_VAF_violin_reviewer_", tag, "_coverage.csv")))

  message("max fraction below 50% VAF (", tag, "): ", signif(max_below50, 3))
}

build_figure(
  file.path(vaf_dir, "Xenium_KRAS_local_VAF_downsampled_100-250um.csv"),
  c(100, 150, 200, 250),
  "100-250um",
  c("100" = "#4575b4", "150" = "#74add1", "200" = "#fdae61", "250" = "#d73027")
)

build_figure(
  file.path(vaf_dir, "Xenium_KRAS_local_VAF_downsampled_100-400um.csv"),
  c(100, 200, 300, 400),
  "100-400um",
  c("100" = "#4575b4", "200" = "#74add1", "300" = "#fdae61", "400" = "#d73027")
)

message("Wrote reviewer violin figures to ", output_dir)
