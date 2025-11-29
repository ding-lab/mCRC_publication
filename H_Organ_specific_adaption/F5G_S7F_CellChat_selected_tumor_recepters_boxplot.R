suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(ggpubr))
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/mCRC_colors.R")
output_dir = getwd()

# Load data
tumor_rds_name = '57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean5.rds'
tumor_rds_path = file.path('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/', 
                           tumor_rds_name)
tumor_obj = readRDS(tumor_rds_path)

metadata_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata = read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL
tumor_obj <- AddMetaData(tumor_obj, metadata)

# EGFR and ERBB2
proinvasive_cell_obj <- tumor_obj %>% 
                        subset(Site_of_Origin %in% c('colon', 'rectum', 'liver', 'lung')) %>% 
                        subset(Patient_ID != 'HT413C1') %>% 
                        subset(cell_type_xenium == 'Non_Canonical_CRC_1')

proinvasive_cell_obj@meta.data <- proinvasive_cell_obj@meta.data %>% 
                           mutate(Organ = case_when(
                                  Site_of_Origin %in% c('rectum', 'colon') ~ 'colorectum',
                                  TRUE ~ Site_of_Origin
                           ))

genes <- c('EGFR', 'ERBB2', 'MET')

# Extract expression matrix
expr_mat <- GetAssayData(proinvasive_cell_obj, slot = "data")[genes, ]

# Build metadata df
meta <- proinvasive_cell_obj@meta.data %>%
  select(orig.ident, Organ) %>%
  mutate(Cell = rownames(.))

# Convert expression to long form
expr_long <- expr_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Cell", values_to = "Expression")

# Merge with metadata
expr_merge <- expr_long %>%
  left_join(meta, by = "Cell")

# Compute per-sample median expression
sample_medians <- expr_merge %>%
  group_by(Organ, orig.ident, gene) %>%
  summarize(median_expr = median(Expression), .groups = "drop")

sample_means <- expr_merge %>%
  group_by(Organ, orig.ident, gene) %>%
  summarize(mean_expr = mean(Expression), .groups = "drop")

head(sample_medians)
head(sample_means)

organ_levels <- c("colorectum", "liver", "lung")

sample_means <- sample_means %>%
  mutate(
    Organ = factor(Organ, levels = organ_levels),
    gene = as.character(gene)
  )

  pvals <- sample_means %>%
  group_by(gene) %>%
  do({
    df <- .
    comps <- combn(organ_levels, 2, simplify = FALSE)

    map_df(seq_along(comps), function(i) {
      cp <- comps[[i]]
      df2 <- df %>% filter(Organ %in% cp)

      tibble(
        gene = df$gene[1],
        group1 = cp[1],
        group2 = cp[2],
        p = wilcox.test(mean_expr ~ Organ, data = df2)$p.value,
        x1 = match(cp[1], organ_levels),
        x2 = match(cp[2], organ_levels),
        x_mid = mean(c(match(cp[1], organ_levels), match(cp[2], organ_levels))),
        y = max(df$mean_expr) + 0.15 * i
      )
    })
  }) %>%
  ungroup()

organ_colors <- c(
  "colorectum" = "#C2B280",
  "liver"      = "brown",
  "lung"       = "steelblue1"
)

p <- ggplot(sample_means, aes(x = Organ, y = mean_expr, fill = Organ)) +
  geom_boxplot(width = 0.7, alpha = 0.7, color = "black") +
  geom_jitter(aes(color = Organ), width = 0.2, alpha = 0.9) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +

  # bracket line
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x1, y = y, yend = y - 0.02 * y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pvals,
    aes(x = x2, xend = x2, y = y, yend = y - 0.02 * y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = pvals,
    aes(x = x_mid, y = y + 0.02 * y,
        label = sprintf("p=%.3g", p)),
    size = 3,
    inherit.aes = FALSE
  ) +

  labs(x = "Organ", y = "Mean expression (per sample)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.border = element_blank(),      # remove full box
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.line.y = element_line(color = "black", linewidth = 0.8)
  ) +
  scale_fill_manual(values = organ_colors) +
  scale_color_manual(values = organ_colors)


pdf(file.path(output_dir, 'EGFR_ERBB2_MET_gene_expression_by_organs.pdf'), width = 10, height = 4)
print(p)
dev.off()

# ROBO1 and IGF1R
stem_tumor_obj <- tumor_obj %>% 
                  subset(Site_of_Origin %in% c('colon', 'rectum', 'liver', 'lung')) %>% 
                  subset(Patient_ID != 'HT413C1') %>% 
                  subset(cell_type_xenium == 'Canonical_CRC_Stem')

stem_tumor_obj@meta.data <- stem_tumor_obj@meta.data %>% 
                            mutate(Organ = case_when(
                                   Site_of_Origin %in% c('rectum', 'colon') ~ 'colorectum',
                                   TRUE ~ Site_of_Origin
                            ))

genes <- c('ROBO1', 'IGF1R')

# Extract expression matrix
expr_mat <- GetAssayData(stem_tumor_obj, slot = "data")[genes, ]

# Build metadata df
meta <- stem_tumor_obj@meta.data %>%
  select(orig.ident, Organ) %>%
  mutate(Cell = rownames(.))

# Convert expression to long form
expr_long <- expr_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "Cell", values_to = "Expression")

# Merge with metadata
expr_merge <- expr_long %>%
  left_join(meta, by = "Cell")

# Compute per-sample median expression
sample_medians <- expr_merge %>%
  group_by(Organ, orig.ident, gene) %>%
  summarize(median_expr = median(Expression), .groups = "drop")

sample_means <- expr_merge %>%
  group_by(Organ, orig.ident, gene) %>%
  summarize(mean_expr = mean(Expression), .groups = "drop")

head(sample_medians)
head(sample_means)

organ_levels <- c("colorectum", "liver", "lung")

sample_means <- sample_means %>%
  mutate(
    Organ = factor(Organ, levels = organ_levels),
    gene = as.character(gene)
  )

pvals <- sample_means %>%
  group_by(gene) %>%
  do({
    df <- .
    comps <- combn(organ_levels, 2, simplify = FALSE)

    map_df(seq_along(comps), function(i) {
      cp <- comps[[i]]
      df2 <- df %>% filter(Organ %in% cp)

      tibble(
        gene = df$gene[1],
        group1 = cp[1],
        group2 = cp[2],
        p = wilcox.test(mean_expr ~ Organ, data = df2)$p.value,
        x1 = match(cp[1], organ_levels),
        x2 = match(cp[2], organ_levels),
        x_mid = mean(c(match(cp[1], organ_levels), match(cp[2], organ_levels))),
        y = max(df$mean_expr) + 0.15 * i
      )
    })
  }) %>%
  ungroup()

p2 <- ggplot(sample_means, aes(x = Organ, y = mean_expr, fill = Organ)) +
  geom_boxplot(width = 0.7, alpha = 0.7, color = "black") +
  geom_jitter(aes(color = Organ), width = 0.2, alpha = 0.9) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +

  # bracket line
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pvals,
    aes(x = x1, xend = x1, y = y, yend = y - 0.02 * y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = pvals,
    aes(x = x2, xend = x2, y = y, yend = y - 0.02 * y),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = pvals,
    aes(x = x_mid, y = y + 0.02 * y,
        label = sprintf("p=%.3g", p)),
    size = 3,
    inherit.aes = FALSE
  ) +

  labs(x = "Organ", y = "Mean expression (per sample)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    panel.border = element_blank(),      # remove full box
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.line.y = element_line(color = "black", linewidth = 0.8)
  ) +
  scale_fill_manual(values = organ_colors) +
  scale_color_manual(values = organ_colors)

pdf(file.path(output_dir, 'IGFR1_ROBO1_gene_expression_by_organs.pdf'), width = 7, height = 4)
print(p2)
dev.off()