#!/usr/bin/env Rscript
# Reviewer-facing panel: are CMETS cells enriched in high local KRAS VAF
# neighborhoods? Restricted to 100-400 um so the radii match the violin
# (the internal version also carries 500 um).
#
# Bins: low_depth (<10 KRAS transcripts nearby), then <75% / 75-90% / >90% VAF.
# y = CMETS% within a bin minus the sample-wide CMETS%, so 0 means "this bin has
# exactly the sample's baseline CMETS content".
# Cramer's V is computed on the three measurable VAF bins only, so it describes
# association with VAF rather than with detection depth.

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

vaf_root <- Sys.getenv(
  "MCRC_KRAS_VAF_DIR",
  unset = "/diskmnt/Projects/MetNet_analysis/Colorectal/mCRC_Manuscript_Script/F_Tumor_subtypes/Xenium_local_VAF_output"
)
cmets_dir <- file.path(vaf_root, "local_VAF_CMETS")
output_dir <- file.path(cmets_dir, "reviewer_100-400um")
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
radii <- c(100, 200, 300, 400)
vaf_levels <- c("low_depth", "<75%", "75-90%", ">90%")
vaf_pal <- c(
  "low_depth" = "grey60",
  "<75%" = "#4575b4",
  "75-90%" = "#fdae61",
  ">90%" = "#d73027"
)
min_bin <- 50

wilson_ci <- function(x, n, z = 1.96) {
  p <- ifelse(n > 0, x / n, NA_real_)
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  tibble(prop = p, lo = pmax(0, center - half), hi = pmin(1, center + half))
}

counts <- read_csv(file.path(cmets_dir, "Xenium_KRAS_local_VAF_CMETS_counts.csv"), show_col_types = FALSE) %>%
  filter(scheme == "absolute", radius_um %in% radii, Sample %in% sample_levels) %>%
  mutate(
    vaf_cat = factor(
      case_when(
        vaf_cat %in% c("<50%", "50-75%") ~ "<75%",
        TRUE ~ vaf_cat
      ),
      levels = vaf_levels
    )
  ) %>%
  group_by(Sample, radius_um, vaf_cat, CMETS_group) %>%
  summarise(n = sum(n), .groups = "drop")

baseline <- counts %>%
  group_by(Sample, radius_um) %>%
  summarise(
    n_total = sum(n),
    baseline_cmets = sum(n[CMETS_group == "CMETS"]) / sum(n),
    .groups = "drop"
  )

enrich_df <- counts %>%
  group_by(Sample, radius_um, vaf_cat) %>%
  summarise(n_cmets = sum(n[CMETS_group == "CMETS"]), n_bin = sum(n), .groups = "drop") %>%
  bind_cols(., wilson_ci(.$n_cmets, .$n_bin)) %>%
  left_join(baseline, by = c("Sample", "radius_um")) %>%
  mutate(
    enrich = prop - baseline_cmets,
    enrich_lo = lo - baseline_cmets,
    enrich_hi = hi - baseline_cmets,
    keep = n_bin >= min_bin,
    Sample = factor(Sample, levels = sample_levels)
  )

# Cramer's V on the measurable VAF bins only (low_depth excluded).
cramers <- counts %>%
  filter(vaf_cat != "low_depth") %>%
  group_by(Sample, radius_um) %>%
  group_modify(function(d, key) {
    tab <- d %>%
      pivot_wider(names_from = CMETS_group, values_from = n, values_fill = 0) %>%
      column_to_rownames("vaf_cat") %>%
      as.matrix()
    tab <- tab[rowSums(tab) >= min_bin, , drop = FALSE]
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(cramers_v = NA_real_, n = sum(tab)))
    }
    ct <- suppressWarnings(chisq.test(tab))
    tibble(
      cramers_v = sqrt(as.numeric(ct$statistic) / (sum(tab) * (min(dim(tab)) - 1))),
      n = sum(tab)
    )
  }) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = sample_levels))

v_labels <- cramers %>%
  group_by(Sample) %>%
  summarise(
    label = paste0(
      "Cramer's V: ",
      paste(ifelse(is.na(cramers_v), "-", sprintf("%.02f", cramers_v)), collapse = " / ")
    ),
    .groups = "drop"
  )

write_csv(enrich_df, file.path(output_dir, "Xenium_KRAS_local_VAF_CMETS_reviewer_enrichment.csv"))
write_csv(cramers, file.path(output_dir, "Xenium_KRAS_local_VAF_CMETS_reviewer_cramersv.csv"))

plot_df <- enrich_df %>% filter(keep)

p <- ggplot(plot_df, aes(x = radius_um, y = enrich, color = vaf_cat, group = vaf_cat)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.4) +
  geom_line(aes(linetype = vaf_cat == "low_depth"), linewidth = 0.7) +
  geom_point(size = 1.9) +
  geom_errorbar(aes(ymin = enrich_lo, ymax = enrich_hi), width = 12, linewidth = 0.35) +
  geom_text(
    data = v_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.4, size = 2.6, color = "grey25"
  ) +
  facet_wrap(~ Sample, ncol = 4) +
  scale_x_continuous(breaks = radii) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = vaf_pal, name = "Local KRAS VAF") +
  scale_linetype_manual(values = c("TRUE" = "22", "FALSE" = "solid"), guide = "none") +
  labs(
    title = "CMETS cells are not concentrated in high-KRAS-VAF neighborhoods",
    subtitle = paste0(
      "Each point: CMETS% among tumor cells whose neighborhood falls in that VAF bin, minus the sample-wide CMETS%. ",
      "Zero = baseline CMETS content; bars are 95% CIs.\n",
      "Grey dashed = cells with <10 KRAS transcripts nearby (below the probe detection limit), shown for completeness. ",
      "Bins with <", min_bin, " cells omitted.\n",
      "Cramer's V (measurable VAF bins only, listed 100/200/300/400 um) quantifies association between local VAF and CMETS status; ",
      "<0.1 is a negligible association."
    ),
    x = "Neighborhood radius (um)",
    y = "CMETS enrichment (percentage points)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 8.5),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8.2, lineheight = 1.15),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_CMETS_reviewer_100-400um.pdf"), p, width = 13, height = 7.2)
ggsave(file.path(output_dir, "Xenium_KRAS_local_VAF_CMETS_reviewer_100-400um.png"), p, width = 13, height = 7.2, dpi = 300)

print(as.data.frame(cramers), row.names = FALSE, digits = 3)
message("Wrote CMETS reviewer panel to ", output_dir)
