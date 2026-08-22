#!/usr/bin/env Rscript
# Compact CMETS vs local KRAS VAF summaries.
# Bins: low_depth, <75%, 75-90%, >90%  (no quartiles; <50% is too sparse).
# Reads count tables from Xenium_KRAS_local_VAF_neighborhood.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

cmets_bin_dir <- Sys.getenv(
  "MCRC_KRAS_CMETS_VAF_DIR",
  unset = file.path(
    Sys.getenv(
      "MCRC_KRAS_VAF_DIR",
      unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_local_VAF_output"
    ),
    "local_VAF_CMETS"
  )
)
dir.create(cmets_bin_dir, recursive = TRUE, showWarnings = FALSE)

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
sample_short <- c(
  "CM1799-A1-Th1Fp1U1" = "CM1799",
  "PE0319C1-A4-S1Fp1U1" = "PE0319",
  "S15-1909-A5U1" = "S15-A5U1",
  "S15-1909-D2U1" = "S15-D2U1",
  "S16-38794-A4U2" = "S16-A4U2",
  "S16-38794-E2U1" = "S16-E2U1",
  "S17-16442-AFR1U1" = "S17-AFR1",
  "S20-46186-A16U1" = "S20-A16"
)
vaf_levels <- c("low_depth", "<75%", "75-90%", ">90%")
vaf_pal <- c(
  "low_depth" = "grey65",
  "<75%" = "#4575b4",
  "75-90%" = "#fdae61",
  ">90%" = "#d73027"
)
cmets_pal <- c("CMETS" = "firebrick2", "Non-CMETS" = "palegreen4")

label_sample <- function(x) {
  factor(unname(sample_short[as.character(x)]), levels = unname(sample_short[sample_levels]))
}

wilson_ci <- function(x, n, z = 1.96) {
  prop <- ifelse(n > 0, x / n, NA_real_)
  lo <- hi <- rep(NA_real_, length(n))
  ok <- !is.na(n) & n > 0
  if (any(ok)) {
    p <- x[ok] / n[ok]
    nn <- n[ok]
    denom <- 1 + z^2 / nn
    center <- (p + z^2 / (2 * nn)) / denom
    half <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2)) / denom
    lo[ok] <- pmax(0, center - half)
    hi[ok] <- pmin(1, center + half)
  }
  tibble(prop = prop, prop_lo = lo, prop_hi = hi)
}

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) suppressPackageStartupMessages(library(patchwork))

raw <- read_csv(
  file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_counts.csv"),
  show_col_types = FALSE
)

bin_counts <- raw %>%
  mutate(
    vaf_cat = factor(
      case_when(
        as.character(vaf_cat) %in% c("<50%", "50-75%", "<75%") ~ "<75%",
        as.character(vaf_cat) %in% c("75-90%") ~ "75-90%",
        as.character(vaf_cat) %in% c(">90%") ~ ">90%",
        TRUE ~ as.character(vaf_cat)
      ),
      levels = vaf_levels
    )
  ) %>%
  group_by(Sample, radius_um, vaf_cat, CMETS_group) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(Sample, radius_um, vaf_cat) %>%
  mutate(n_bin = sum(n), prop_in_bin = ifelse(n_bin > 0, n / n_bin, NA_real_)) %>%
  group_by(Sample, radius_um, CMETS_group) %>%
  mutate(n_group = sum(n), prop_in_group = ifelse(n_group > 0, n / n_group, NA_real_)) %>%
  ungroup() %>%
  mutate(
    Sample = factor(Sample, levels = sample_levels),
    Sample_lab = label_sample(Sample),
    radius_lab = factor(paste0(radius_um, " um"), levels = paste0(c(100, 200, 300, 400, 500), " um"))
  )

overall <- bin_counts %>%
  group_by(Sample, Sample_lab, radius_um, radius_lab) %>%
  summarise(
    n_cells = sum(n),
    n_cmets = sum(n[CMETS_group == "CMETS"]),
    overall_cmets_frac = n_cmets / n_cells,
    .groups = "drop"
  )

cmets_in_bin <- bin_counts %>%
  filter(CMETS_group == "CMETS") %>%
  transmute(Sample, Sample_lab, radius_um, radius_lab, vaf_cat, n_cmets = n, n_bin) %>%
  mutate(wilson_ci(n_cmets, n_bin)) %>%
  left_join(overall, by = c("Sample", "Sample_lab", "radius_um", "radius_lab")) %>%
  mutate(
    enrich = prop - overall_cmets_frac,
    keep = n_bin >= 50
  )

write.csv(
  cmets_in_bin,
  file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_bins.csv"),
  row.names = FALSE
)

heat_df <- cmets_in_bin %>%
  filter(radius_um %in% c(200, 500)) %>%
  mutate(
    fill_enrich = ifelse(keep, enrich, NA_real_),
    txt = ifelse(keep, percent(enrich, accuracy = 0.1), "n<50")
  )

p_heat <- ggplot(heat_df, aes(x = vaf_cat, y = Sample_lab, fill = fill_enrich)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = txt), size = 2.7) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    na.value = "grey90",
    labels = percent_format(accuracy = 1),
    name = "CMETS%\nvs sample"
  ) +
  facet_wrap(~ radius_lab, ncol = 2) +
  labs(
    title = "CMETS enrichment in each local VAF category",
    subtitle = "Color / text = CMETS% in that bin minus the sample-wide CMETS%. Grey / n<50 = too few cells. 200 um vs 500 um: if the pattern holds, it is not sparse-count noise.",
    x = "Local KRAS VAF",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "grey95", color = NA)
  )

p_stack <- ggplot(
  bin_counts %>% filter(radius_um == 500),
  aes(x = vaf_cat, y = n, fill = CMETS_group)
) +
  geom_col(position = "fill", width = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = cmets_pal, name = NULL) +
  facet_wrap(~ Sample_lab, ncol = 4) +
  labs(
    title = "500 um: CMETS vs Non-CMETS mix inside each VAF category",
    subtitle = "Each bar is 100% of tumor cells in that sample x bin. This is the comparison requested: CMETS proportion vs Non-CMETS proportion per category.",
    x = "Local KRAS VAF",
    y = "Proportion of cells in bin"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

p_frac <- ggplot(
  cmets_in_bin %>% filter(radius_um == 500),
  aes(x = vaf_cat, y = prop)
) +
  geom_hline(
    data = overall %>% filter(radius_um == 500),
    aes(yintercept = overall_cmets_frac),
    linetype = "dashed",
    color = "grey30",
    linewidth = 0.4
  ) +
  geom_col(aes(fill = vaf_cat), width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin = prop_lo, ymax = prop_hi), width = 0.2, linewidth = 0.35) +
  geom_text(
    data = ~ dplyr::filter(.x, n_bin > 0),
    aes(label = n_bin, y = pmax(dplyr::coalesce(prop_hi, prop), 0)),
    vjust = -0.3,
    size = 2.1
  ) +
  scale_fill_manual(values = vaf_pal, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.16))) +
  facet_wrap(~ Sample_lab, ncol = 4) +
  labs(
    title = "500 um: CMETS fraction in each VAF category",
    subtitle = "Dashed line = overall CMETS fraction in that sample. Error bars = Wilson 95% CI. Numbers = cells in bin.",
    x = "Local KRAS VAF",
    y = "CMETS proportion"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "grey95")
  )

p_line <- ggplot(
  cmets_in_bin %>% filter(keep),
  aes(x = factor(radius_um), y = enrich, group = vaf_cat, color = vaf_cat)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = vaf_pal, name = NULL, breaks = vaf_levels) +
  facet_wrap(~ Sample_lab, ncol = 4) +
  labs(
    title = "Does CMETS enrichment persist at larger radius?",
    subtitle = "Y = CMETS% in bin minus sample-wide CMETS%. Includes low_depth. Bins with < 50 cells omitted.",
    x = "Neighborhood radius (um)",
    y = "CMETS enrichment"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_heatmap.pdf"), p_heat, width = 10, height = 5.2)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_heatmap.png"), p_heat, width = 10, height = 5.2, dpi = 300)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_composition_500um.pdf"), p_stack, width = 11, height = 6.2)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_composition_500um.png"), p_stack, width = 11, height = 6.2, dpi = 300)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_fraction_500um.pdf"), p_frac, width = 11, height = 6.4)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_fraction_500um.png"), p_frac, width = 11, height = 6.4, dpi = 300)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_vs_radius.pdf"), p_line, width = 11, height = 6.4)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_vs_radius.png"), p_line, width = 11, height = 6.4, dpi = 300)

p_line_vafx <- ggplot(
  cmets_in_bin %>% filter(keep),
  aes(x = vaf_cat, y = enrich, group = radius_um, color = factor(radius_um))
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(
    values = c(
      "100" = "grey65",
      "200" = "#4575b4",
      "300" = "#74add1",
      "400" = "#fdae61",
      "500" = "#d73027"
    ),
    name = "Radius (um)"
  ) +
  facet_wrap(~ Sample_lab, ncol = 4) +
  labs(
    title = "CMETS enrichment across VAF categories, by neighborhood radius",
    subtitle = "Y = CMETS% in bin minus sample-wide CMETS%. Each line is one radius. Bins with < 50 cells omitted.",
    x = "Local KRAS VAF category",
    y = "CMETS enrichment"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_vs_radius_vafcat.pdf"), p_line_vafx, width = 11, height = 6.4)
ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary_enrich_vs_radius_vafcat.png"), p_line_vafx, width = 11, height = 6.4, dpi = 300)

if (has_patchwork) {
  p_page <- p_heat / p_stack / p_line +
    plot_annotation(
      title = "CMETS vs local KRAS VAF — low_depth / <75% / 75-90% / >90%",
      subtitle = "<50% was too sparse, so it is pooled with 50-75% into <75%. 75-90% vs >90% splits the high-VAF neighborhoods where most cells sit.",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30")
      )
    )
  ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary.pdf"), p_page, width = 12, height = 14)
  ggsave(file.path(cmets_bin_dir, "Xenium_KRAS_local_VAF_CMETS_summary.png"), p_page, width = 12, height = 14, dpi = 300)
}

cat("Summary figures written to", cmets_bin_dir, "\n")
