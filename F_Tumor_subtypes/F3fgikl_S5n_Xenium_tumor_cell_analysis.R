# ---- Libraries ----
library(tidyverse)
library(Seurat)
library(ggpubr)
library(qs)
library(purrr)
library(MAST)
library(ggrepel)
library(RColorBrewer)
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")

# ---- Setup ----
output_dir <- getwd()

# ---- Color Definitions ----
nb_col <- c(
  'NB2_tumor_core' = 'hotpink4', 
  'NB3_tumor_body' = 'violetred', 
  'NB4_tumor_stroma_interface' = 'orange1'
)

pm_col <- c(
  Primary = 'navy', 
  Metastasis = 'violetred'
)

# ---- Load Metadata ----
metadata_path <- '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
meta_df <- read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"
meta_df <- meta_df %>% column_to_rownames("barcode")

# ---- Filter and Annotate Tumor Cells ----
meta_tmr_df <- meta_df %>% 
               dplyr::filter(neighborhoods %in% c('NB2_tumor_core', 'NB3_tumor_body', 'NB4_tumor_stroma_interface')) %>% 
               dplyr::filter(All_cell_type2 %in% c('Intestine-like', 'Non-canonical', 'Proliferative-like', 'Stem-like')) %>% 
               mutate(tumor_subtype = case_when(
                   All_cell_type2 == 'Non-canonical' & neighborhoods == 'NB4_tumor_stroma_interface' ~ 'Border_non-canonical',
                   All_cell_type2 == 'Non-canonical' & neighborhoods == 'NB2_tumor_core' ~ 'Core_non-canonical',
                   TRUE ~ All_cell_type2
               ))

# ---- Figure 3G: Non-canonical cell proportion by neighborhood ----
nb_tumor_prop <- meta_tmr_df %>%
                 dplyr::count(neighborhoods, All_cell_type2) %>%
                 group_by(neighborhoods) %>%
                 mutate(proportion = n / sum(n)) %>%
                 ungroup() %>% 
                 filter(All_cell_type2 == 'Non-canonical')

ggplot(nb_tumor_prop, aes(x = neighborhoods, y = proportion, fill = All_cell_type2)) +
  geom_bar(stat = "identity") +
  labs(x = "NB", y = "Proportion", fill = "Cell type") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = c("Non-canonical" = "#be0032")) -> BarP

pdf(file=file.path(output_dir, 'F3g_mCRC_N26_TMR_NB_tumor_barplot2.pdf'), width=5, height=6) 
print(BarP)
dev.off()

# ---- Figure 3I: Neighborhood proportions by organ ----
nb_counts <- meta_tmr_df %>%
             filter(Organ != 'Breast') %>% 
             dplyr::count(sample_id, neighborhoods, Organ)

sample_info <- meta_tmr_df %>%
               dplyr::select(sample_id, Tissue_ID) %>%
               distinct()

nb_counts <- left_join(nb_counts, sample_info, by = "sample_id")

sample_total <- nb_counts %>%
  group_by(sample_id) %>%
  summarise(total = sum(n), .groups = "drop")

nb_prop_df <- nb_counts %>%
  left_join(sample_total, by = "sample_id") %>%
  mutate(proportion = n / total)

ggplot(nb_prop_df, aes(x = Organ, y = proportion, fill = neighborhoods, color = neighborhoods)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 1.5, aes(color = neighborhoods)) +
  geom_line(data = nb_prop_df,
            aes(x = Organ, y = proportion, group = Tissue_ID),
            linewidth = 0.4, alpha = 0.5, color = "gray40") +
  facet_wrap(~ neighborhoods, scales = "free_y") +
  scale_fill_manual(values = nb_col) +
  scale_color_manual(values = nb_col) +
  theme_minimal(base_size = 13) +
  labs(y = "Proportion", x = "Organ") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")) +
  stat_compare_means(
    comparisons = list(c("Colon", "Liver"), c("Colon", "Lung"), c("Liver", "Lung")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.01
  ) -> BoxP

pdf(file=file.path(output_dir, 'F3i_mCRC_N26_TMR_NB_tumor_boxplot_COLbyNB.pdf'), width=12, height=6) 
print(BoxP)
dev.off()



# ---- Figure 3F: Non-canonical cells in Primary vs Metastasis ----
cell_counts2 <- meta_tmr_df %>%
               dplyr::count(sample_id, All_cell_type2, Organ) %>% 
               mutate(PM_status = case_when(
                   Organ == 'Colon' ~ 'Primary',
                   TRUE ~ 'Metastasis'
               )) 

sample_info <- meta_tmr_df %>%
               dplyr::select(sample_id, Tissue_ID) %>%
               distinct()

cell_counts2 <- left_join(cell_counts2, sample_info, by = "sample_id")

cell_total2 <- cell_counts2 %>%
               group_by(sample_id) %>%
               summarise(total = sum(n), .groups = "drop")

cell_prop_df2 <- cell_counts2 %>%
                 left_join(cell_total2, by = "sample_id") %>%
                 mutate(proportion = n / total) %>% 
                 filter(All_cell_type2 == 'Non-canonical') 

cell_prop_df2$PM_status <- factor(cell_prop_df2$PM_status, levels = c('Primary', 'Metastasis'))    

ggplot(cell_prop_df2, 
       aes(x = PM_status, y = proportion, fill = PM_status, color = PM_status)) +
       geom_boxplot(outlier.shape = NA, alpha = 0.75) +
       geom_point(size = 1.5, aes(color = PM_status)) +
       geom_line(data = cell_prop_df2,
                 aes(x = PM_status, y = proportion, group = Tissue_ID),
                 linewidth = 0.4, alpha = 0.5, color = "gray40") +
       scale_fill_manual(values = pm_col) +
       scale_color_manual(values = pm_col) +
       theme_minimal(base_size = 13) +
       labs(y = "Proportion")  +
       theme(axis.text.x = element_text(angle = 45, hjust = 1),
             panel.grid = element_blank(),
             axis.line = element_line(color = "black"),
             axis.ticks = element_line(color = "black")) + 
       stat_compare_means(
           method = "t.test",
           label = "p.format",
           tip.length = 0.01
       ) -> BoxP3

pdf(file=file.path(output_dir, 'F3f_mCRC_N26_TMR_Primary_Met_tumor_boxplot.pdf'), width=5, height=4) 
print(BoxP3)
dev.off()

# ---- Figure 3K: Non-canonical cells in NB4 (tumor-stroma interface) ----
meta_nb4 <- meta_tmr_df %>%
  filter(neighborhoods == "NB4_tumor_stroma_interface")

cell_counts_nb4 <- meta_nb4 %>%
  dplyr::count(sample_id, All_cell_type2, Organ) %>%
  mutate(PM_status = case_when(
    Organ == "Colon" ~ "Primary",
    TRUE ~ "Metastasis"
  ))

sample_info <- meta_nb4 %>%
  dplyr::select(sample_id, Tissue_ID) %>%
  distinct()

cell_counts_nb4 <- left_join(cell_counts_nb4, sample_info, by = "sample_id")


cell_total_nb4 <- cell_counts_nb4 %>%
  group_by(sample_id) %>%
  summarise(total = sum(n), .groups = "drop")

cell_prop_nb4 <- cell_counts_nb4 %>%
  left_join(cell_total_nb4, by = "sample_id") %>%
  mutate(proportion = n / total) %>%
  filter(All_cell_type2 == "Non-canonical")

cell_prop_nb4$PM_status <- factor(cell_prop_nb4$PM_status, levels = c("Primary", "Metastasis"))

cell_prop_nb4 <- cell_prop_nb4 %>% dplyr::filter(total > 100)


ViolinP_nb4 <- ggplot(cell_prop_nb4,
       aes(x = PM_status, y = proportion, fill = PM_status, color = PM_status)) +
  geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
  geom_point(size = 1.8, alpha = 0.9) +
  geom_line(aes(group = Tissue_ID), color = "gray50", alpha = 0.6, linewidth = 0.4) +
  scale_fill_manual(values = pm_col) +
  scale_color_manual(values = pm_col) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Non-canonical cells in NB4_tumor_stroma_interface",
    y = "Proportion",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    tip.length = 0.01
  )

pdf(file=file.path(output_dir, 'F3k_mCRC_N26_NB4_Primary_Met_tumor_Proinvasive_boxplot.pdf'), width=6, height=4) 
print(ViolinP_nb4)
dev.off()

# ---- Figure 3L: Dot plot of border vs core markers ----
# Load tumor Seurat object
tumor_obj <- qread('/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tumor_SketchInt.qs')
tumor_obj <- AddMetaData(tumor_obj, meta_tmr_df)

# Note: DEG analysis is commented out but can be run if needed
# Idents(tumor_obj) <- tumor_obj$tumor_subtype
# sample_deg <- FindAllMarkers(tumor_obj, 
#                              test.use = 'MAST',
#                              thresh.use = 0.25,  
#                              min.pct = 0.2,  
#                              only.pos = TRUE, 
#                              return.thresh = 0.05) 
# write.csv(sample_deg, file.path(output_dir, 'Xenium_tumor_core_border_subtype_deg.csv'))

noncan_tumor <- tumor_obj %>% 
               subset(tumor_subtype %in% c('Border_non-canonical', 'Stem-like', 'Core_non-canonical')) 

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-1.5, -0.75, 0, 0.75, 1.5)

core_border_markers <- c(
  'UBD', 'RBP1', 'KLK10', 'CXCL2', 'LIF', 'LAMA3',
  'ANXA1', 'SLC2A1', 'EMP1', 'PLAUR', 'VEGFA',
  'NDRG1', 'PLOD2', 'P4HA1', 'HK2', 'KDM3A'
)

noncan_tumor$tumor_subtype <- factor(
  noncan_tumor$tumor_subtype,
  levels = c('Border_non-canonical', 'Stem-like', 'Core_non-canonical')
)

p <- DotPlot(noncan_tumor, features = core_border_markers, group.by = 'tumor_subtype') +
     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
     scale_color_gradientn(
       colors = rdylbu_colors, 
       limits = c(-1.5, 1.5),   
       breaks = color_breaks
     ) +         
     scale_size_area(limits = c(0, 100), oob = scales::squish)

pdf(file=file.path(output_dir, 'F3l_Dotplot_border_vs_core_proinvasive_cells.pdf'), width=10, height=4) 
print(p)
dev.off()


# ---- Figure S5M: Volcano plot of border vs core DEGs ----
all_subtype_deg <- read.csv(file.path(output_dir, 'Xenium_tumor_core_border_subtype_deg.csv'), row.names = 1)
border_core_noncan_deg <- read.csv(file.path(output_dir, 'Xenium_tumor_core_border_non_canonical_deg.csv'), row.names = 1)
border_core_noncan_deg$gene <- rownames(border_core_noncan_deg)

noncangenes <- all_subtype_deg %>% 
               filter(cluster %in% c('Border_non-canonical', 'Non-canonical', 'Core_non-canonical'))
noncangenes <- unique(noncangenes$gene)

border_core_noncan_deg2 <- border_core_noncan_deg %>% 
                          filter(gene %in% noncangenes)

# Handle zero p-values and define significance thresholds
border_core_noncan_deg2 <- border_core_noncan_deg2 %>%
  mutate(
    p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj)
  )

# Define thresholds
fc_thresh <- 0.25      # log2 fold-change cutoff
padj_thresh <- 0.05    # adjusted p-value cutoff

# Add significance category
border_core_noncan_deg2 <- border_core_noncan_deg2 %>%
  mutate(
    sig = case_when(
      p_val_adj < padj_thresh & avg_log2FC >  fc_thresh ~ "Up-regulated",
      p_val_adj < padj_thresh & avg_log2FC < -fc_thresh ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Pick top labeled genes for annotation
df_top <- border_core_noncan_deg2 %>%
  arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 5) %>%
  bind_rows(
    border_core_noncan_deg2 %>%
      arrange(avg_log2FC) %>% 
      slice_head(n = 5)
  )

# Define axis limits
x_lim <- max(abs(border_core_noncan_deg2$avg_log2FC), na.rm = TRUE)
y_lim <- max(-log10(border_core_noncan_deg2$p_val_adj), na.rm = TRUE)

# Volcano plot
p <- ggplot(border_core_noncan_deg2, 
            aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(data = subset(border_core_noncan_deg2, sig == "Not significant"),
             color = "gray80", size = 1.5, alpha = 0.5) +
  geom_point(data = subset(border_core_noncan_deg2, sig != "Not significant"),
             size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-fc_thresh, 0, fc_thresh), 
             linetype = c("dotted", "dashed", "dotted"), 
             color = c("gray60", "gray50", "gray60"), linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_thresh), 
             linetype = "dotted", color = "gray50", linewidth = 0.4) +
  geom_text_repel(
    data = df_top,
    aes(label = gene),
    size = 3.3,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "gray60",
    segment.size = 0.3,
    seed = 42
  ) +
  scale_color_manual(
    values = c(
      "Up-regulated"   = "#D73027",
      "Down-regulated" = "#4575B4",
      "Not significant" = "gray70"
    )
  ) +
  labs(
    title = "Volcano Plot — Spatial Expression Gradient",
    subtitle = paste0("Thresholds: |log2FC| > ", fc_thresh, ", adj.p < ", padj_thresh),
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("Adj. p-value")),
    color = NULL
  ) +
  coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim * 1.05)) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pdf(file=file.path(output_dir, 'S5n_Volcano_border_vs_core_proinvasive_cells.pdf'), width=8, height=4) 
print(p)
dev.off()
