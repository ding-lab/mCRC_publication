library(qs)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

out_dir = getwd()
obj_path = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tmr_association_fibroblast_obj.qs"
obj <- qread(obj_path)

obj@meta.data <- obj@meta.data %>% 
                 mutate(near_tmr_region = case_when(
                     Organ == 'Liver' ~ 'metastasis',
                     Organ == 'Colon' & (near_mucosa_TMR == TRUE | tmr_region == 'mucosa') ~ 'mucosa',
                     Organ == 'Colon' & (near_submucosa_TMR == TRUE | tmr_region == 'submucosa')  ~ 'submucosa',
                     Organ == 'Colon' & (near_muscle_TMR == TRUE | tmr_region == 'muscle')  ~ 'muscle'))

obj@meta.data <- obj@meta.data %>% 
                        mutate(sample_region = paste0(Tissue_ID, '_', near_tmr_region))


# Set your desired orders
cell_order <- c("WNT5A_BMP", "WNT5A_infl", "iCAF", "stromal_fibroblast", "mCAF")
region_order <- c("mucosa", "submucosa", "muscle", "metastasis")

# Prepare data
prop_df <- obj@meta.data %>%
  mutate(
    All_cell_type2 = factor(All_cell_type2, levels = cell_order),
    near_tmr_region = factor(near_tmr_region, levels = region_order)
  ) %>%
  group_by(near_tmr_region, All_cell_type2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(near_tmr_region) %>%
  mutate(prop = n / sum(n))

# Plot
p <- ggplot(prop_df,
            aes(x = near_tmr_region,
                y = prop,
                fill = All_cell_type2)) +
  geom_col(color = "black", width = 0.75) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "WNT5A_BMP" = "#8dd3c7",
    "WNT5A_infl" = "#ffffb3",
    "iCAF" = "#bebada",
    "stromal_fibroblast" = "#fb8072",
    "mCAF" = "#80b1d3"
  )) +
  labs(
    x = "Near TMR Region",
    y = "Cell proportion",
    fill = "Cell type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

pdf(file.path(out_dir, "barplot_fibroblast_subtype_by_regions.pdf"), width = 5, height = 4)
p
dev.off()

sample_prop <- obj@meta.data %>%
  filter(near_tmr_region %in% region_order) %>%      # remove mucosa
  mutate(
    All_cell_type2 = factor(All_cell_type2, levels = cell_order),
    near_tmr_region = factor(near_tmr_region, levels = region_order)
  ) %>%
  group_by(corrected_sample_id, near_tmr_region, All_cell_type2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(corrected_sample_id, near_tmr_region) %>%
  mutate(prop = n / sum(n))

# ---------------------------------------------------------
#  INPUT: sample_prop (already computed earlier)
#  sample_prop contains: corrected_sample_id, near_tmr_region,
#                        All_cell_type2, prop
# ---------------------------------------------------------

fib_types <- c("WNT5A_BMP", "WNT5A_infl", "iCAF", "stromal_fibroblast", "mCAF")

plots <- lapply(fib_types, function(ct) {

  df_ct <- sample_prop %>% filter(All_cell_type2 == ct)

  # ---- Dynamic y-axis based on data ----
  data_max <- max(df_ct$prop, na.rm = TRUE)

  # panel-specific fallback buffer
  buffer <- max(0.05, data_max * 0.30)

  ymax <- data_max + buffer

  # ---- SAFE p-value heights (fractions of ymax) ----
  y_pos1 <- ymax * 0.95   # top comparison
  y_pos2 <- ymax * 0.85   # middle comparison
  y_pos3 <- ymax * 0.75   # bottom comparison

  ggplot(df_ct,
         aes(x = near_tmr_region, y = prop,
             fill = near_tmr_region, color = near_tmr_region)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_point(size = 2, alpha = 0.8) +
    
    # ----- Three Comparisons (Wilcoxon) -----
    stat_compare_means(
      comparisons = list(c("submucosa", "muscle")),
      method = "wilcox.test",
      label = "p.format",
      label.y = y_pos1
    ) +
    stat_compare_means(
      comparisons = list(c("muscle", "metastasis")),
      method = "wilcox.test",
      label = "p.format",
      label.y = y_pos2
    ) +
    stat_compare_means(
      comparisons = list(c("submucosa", "metastasis")),
      method = "wilcox.test",
      label = "p.format",
      label.y = y_pos3
    ) +

    # ---- Scales ----
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, ymax)
    ) +
    scale_fill_manual(values = c(
      "submucosa"   = "#c7e9b4",
      "muscle"      = "#7fcdbb",
      "metastasis"  = "#41b6c4"
    )) +
    scale_color_manual(values = c(
      "submucosa"   = "#c7e9b4",
      "muscle"      = "#7fcdbb",
      "metastasis"  = "#41b6c4"
    )) +

    # ---- Labels & Theme ----
    labs(
      title = ct,
      x = "Near TMR Region",
      y = "Cell proportion"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(),
      axis.ticks = element_line(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
})

final_plot <- wrap_plots(plots, nrow = 1)


# ---------------------------------------------------------
#  Combine all 5 fibroblast panels into ONE ROW
# ---------------------------------------------------------

pdf(file.path(out_dir, "boxplot_fibroblast_subtype_by_regions.pdf"), width = 20, height = 4)
final_plot
dev.off()
