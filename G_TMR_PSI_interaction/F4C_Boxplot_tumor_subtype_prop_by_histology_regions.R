library(qs)
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(Matrix)
library(circlize)
library(viridis)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

paired_CRC_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tumor_paired_CoLi_TMR_association.qs'
paired_CRC <- qread(paired_CRC_path)
paired_CRC

paired_CRC@meta.data <- paired_CRC@meta.data %>% 
                        mutate(tmr_region = case_when(
                                   corrected_sample_id == 'CM819C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM819C1_Co' & tn_label %in% c(7, 8, 12, 14)  ~ 'muscle',
                                   corrected_sample_id == 'CM819C1_Co' & tn_label %in% c(47, 48, 49)  ~ 'mucosa',
                                   corrected_sample_id == 'CM819C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM579C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM579C1_Co' & tn_label %in% c(24) ~ 'muscle',
                                   corrected_sample_id == 'CM579C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM798C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM798C1_Co' & tn_label %in% c(1:9) ~ 'muscle',
                                   corrected_sample_id == 'CM798C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM397C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM397C1_Co' & tn_label %in% c(9, 10, 12, 13) ~ 'muscle',
                                   corrected_sample_id == 'CM397C1_Co' ~ 'submucosa'
                                   ))

tumor_subtype_prop = 
    paired_CRC@meta.data %>%
    filter(All_cell_type2 %in% c('Intestine-like', 'Proliferative-like', 'Stem-like', 'Non-canonical')) %>% 
    filter(tmr_region != 'mucosa') %>% 
    count(Case_ID, tmr_region, All_cell_type2) %>%
    group_by(Case_ID, tmr_region) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup() 

tumor_subtype_prop$tmr_region <- factor(tumor_subtype_prop$tmr_region, 
                                        level = c("submucosa", "muscle", "metastasis"))

tumor_subtype_prop$All_cell_type2 <- factor(tumor_subtype_prop$All_cell_type2, 
                                        level = c('Intestine-like', 'Proliferative-like', 'Stem-like', 'Non-canonical'))

tmr_col <- c(submucosa = 'palegreen4', muscle = 'firebrick1', metastasis = 'violetred1')

BoxP1 <- ggplot(tumor_subtype_prop, 
                aes(x = tmr_region, y = proportion, fill = tmr_region, color = tmr_region)) +
         geom_boxplot(outlier.shape = NA, alpha = 0.6) +
         geom_point(size = 1.5, aes(color = tmr_region)) +
         geom_line(data = tumor_subtype_prop,
                   aes(x = tmr_region, y = proportion, group = Case_ID),
                   linewidth = 0.4, alpha = 0.5, color = "gray40") +
         facet_wrap(~ All_cell_type2, scales = "free_y") +
         scale_fill_manual(values = tmr_col) +
         scale_color_manual(values = tmr_col) +
         theme_minimal(base_size = 13) +
         labs(y = "Proportion", x = "TMR region") +
  stat_compare_means(
    comparisons = list(c("submucosa", "muscle"), c("muscle", "metastasis"), c("submucosa", "metastasis")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.01
  ) +
  theme(
    panel.grid = element_blank(),       # Remove all grid lines
    axis.line = element_line(),         # Add x and y axis lines
    axis.ticks = element_line(),        # Add tick marks
    axis.text.x = element_text(angle = 45, hjust = 1) # Optional: angled x-axis labels
  )

BoxP2 <- ggplot(tumor_subtype_prop %>% filter(All_cell_type2 %in% c('Stem-like', 'Non-canonical')), 
                aes(x = tmr_region, y = proportion, fill = tmr_region, color = tmr_region)) +
         geom_boxplot(outlier.shape = NA, alpha = 0.6) +
         geom_point(size = 1.5, aes(color = tmr_region)) +
         geom_line(data = tumor_subtype_prop %>% filter(All_cell_type2 %in% c('Stem-like', 'Non-canonical')),
                   aes(x = tmr_region, y = proportion, group = Case_ID),
                   linewidth = 0.4, alpha = 0.5, color = "gray40") +
         facet_wrap(~ All_cell_type2, scales = "free_y") +
         scale_fill_manual(values = tmr_col) +
         scale_color_manual(values = tmr_col) +
         theme_minimal(base_size = 13) +
         labs(y = "Proportion", x = "TMR region") +
  stat_compare_means(
    comparisons = list(c("submucosa", "muscle"), c("muscle", "metastasis"), c("submucosa", "metastasis")),
    method = "wilcox.test",
    label = "p.format",
    tip.length = 0.01
  ) +
  theme(
    panel.grid = element_blank(),       # Remove all grid lines
    axis.line = element_line(),         # Add x and y axis lines
    axis.ticks = element_line(),        # Add tick marks
    axis.text.x = element_text(angle = 45, hjust = 1) # Optional: angled x-axis labels
  )


out_dir = getwd()
pdf(file=file.path(out_dir, 'mCRC_N26_TMR_tumor_boxplot_by_tmr_regions.pdf'), width=8, height=8) 
BoxP1
dev.off()