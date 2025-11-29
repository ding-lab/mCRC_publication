library(tidyverse)
library(ggplot2)
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
output_dir = getwd()

xenium_ranknet = read.csv('Xenium_cellchat_ranknet_liver_vs_lung_tumor_targeted.csv')
snRNA_ranknet = read.csv('snRNA_cellchat_ranknet_liver_vs_lung_tumor_targeted.csv')
common_pathways = intersect(xenium_ranknet$name, snRNA_ranknet$name)
xenium_ranknet <- xenium_ranknet %>% filter(name %in% common_pathways)
snRNA_ranknet <- snRNA_ranknet %>% filter(name %in% common_pathways)

xenium_ranknet_wide <- xenium_ranknet %>% select(name, group, contribution) %>% 
                       pivot_wider(id_cols = name, names_from = group, values_from = contribution, values_fn = mean)
xenium_ranknet_wide$difference = xenium_ranknet_wide$Liver - xenium_ranknet_wide$Lung
snRNA_ranknet_wide <- snRNA_ranknet %>% select(name, group, contribution) %>% 
                       pivot_wider(id_cols = name, names_from = group, values_from = contribution, values_fn = mean)
snRNA_ranknet_wide$difference = snRNA_ranknet_wide$Liver - snRNA_ranknet_wide$Lung

merged_ranknet <- merge(
  xenium_ranknet_wide[, c("name", "difference")],
  snRNA_ranknet_wide[, c("name", "difference")],
  by = "name",
  suffixes = c("_xenium", "_snRNA")
)

# Keep only pathways with the same sign in both
same_trend <- merged_ranknet %>%
  filter(sign(difference_xenium) == sign(difference_snRNA))

# If you only want the pathway names
same_trend_pathways <- same_trend$name

xenium_ranknet_wide <- xenium_ranknet_wide %>% filter(name %in% same_trend_pathways)
snRNA_ranknet_wide <- snRNA_ranknet_wide %>% filter(name %in% same_trend_pathways)

snRNA_ranknet_wide$color <- ifelse(snRNA_ranknet_wide$difference < 0, "steelblue1", "brown")

ggplot(snRNA_ranknet_wide, aes(x = reorder(name, difference), y = difference, fill = color)) +
  geom_col() +
  coord_flip() +
  scale_fill_identity() +
  labs(
    title = "Liver vs. Lung Signaling Difference",
    x = "Pathway",
    y = "Liver - Lung Contribution"
  ) + 
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Keep axes
    axis.ticks = element_line(color = "black")  # Optional: emphasize axis ticks
  ) -> p1

pdf('cellchat_consistent_snRNA_Xenium_barplot_liver_vs_lung.pdf', width = 3, height = 4)
p1
dev.off()