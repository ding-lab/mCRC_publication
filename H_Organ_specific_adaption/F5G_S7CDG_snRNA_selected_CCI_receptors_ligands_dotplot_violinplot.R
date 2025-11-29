library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(RColorBrewer)
library(scales)
library(rstatix)
library(scales)
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
output_dir = getwd()

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
mCRC_sn_all = readRDS('/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean5.rds')
mCRC_sn_all

metadata_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata = read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL
mCRC_sn_all <- AddMetaData(mCRC_sn_all, metadata)

mCRC_sn_simplified_cell <- mCRC_sn_all %>% 
                           subset(broad_cell_type %in% c('T cell', 'Tumor', 'Macrophage', 
                                                         'Endothelial cell', 'Hepatocyte', 'Fibroblast', 'Plasma', 
                                                         'B cell', 'NK cell', 'Alveolar cell'
                            ))  %>%
                           subset(!(cell_type_xenium %in% c('APCDD1_CRC', 'mesenchymal_unassigned', 'Myeloid_doublets',
                                                            'TNK_doublets', 'mast', 'smooth_muscle'))) %>%
                           subset(Site_of_Origin %in% c('rectum', 'colon', 'liver', 'lung')) %>%
                           subset(Site_of_Origin %in% c('liver', 'lung'))

mCRC_sn_simplified_cell@meta.data <- mCRC_sn_simplified_cell@meta.data %>% 
                                     mutate(Organ = case_when(Site_of_Origin %in% c('rectum', 'colon') ~ 'colorectum',
                                                              TRUE ~ Site_of_Origin))

mCRC_sn_simplified_cell$cell_type_xenium <- factor(mCRC_sn_simplified_cell$cell_type_xenium,
level = c('Canonical_CRC_Intestine', 'Canonical_CRC_Intestine_Proliferation', 'Canonical_CRC_Stem_Proliferation', 'Canonical_CRC_Stem', 'Non_Canonical_CRC_1',
          'Hepatocyte', 'Cholangiocyte', 'Alveolar type1', 'Alveolar type2',
          'WNT5A_BMP', 'WNT5A_infl', 'mCAF', 'iCAF', 'stromal_fibroblast', 'stellate_cell', 'LSEC', 'Pericyte',
          'TAM_hypoxic-angiogenic', 'TAM', 'APC-like', 'B cell', 'Plasma', 'CD56_low_NK', 'CD56_high_NK',
          'CD8_NK_like_T', 'CD8_EM_T', 'CD8_TRM_T',  'CD8_Effector_T', 'CD4_Resting_T', 'CD4_naive/CM_T', 'CD4_Th17_T', 'CD4_Treg', 
          'Vascular endothelial', 'Lymphatic endothelial'))


mCRC_sn_simplified_cell@meta.data <- mCRC_sn_simplified_cell@meta.data %>% 
                                     mutate(broad_cell_type2 = 
                                            case_when(broad_cell_type %in% c('T cell', 'NK cell') ~ 'T/NK',
                                                      broad_cell_type %in% c('B cell', 'Plasma') ~ 'B/Plasma',
                                                      TRUE ~ broad_cell_type),
                                            broad_cell_type_organ = paste0(Organ, '_', broad_cell_type2),
                                            xenium_cell_type_organ = paste0(Organ, '_', cell_type_xenium)
                                           )


mCRC_sn_simplified_cell$broad_cell_type2 <- factor(mCRC_sn_simplified_cell$broad_cell_type2,
                                                   level = c(
                                                       'Tumor', 'Fibroblast', 'Macrophage', 'Hepatocyte', 'Alveolar cell',
                                                       'Endothelial cell', 'B/Plasma', 'T/NK')
                                                   )


p1 = DotPlot(mCRC_sn_simplified_cell, features=c('SLIT2', 'HGF', 'HBEGF', 'ROBO1', 'MET', 'EGFR', 'ERBB2'), 
            group.by='broad_cell_type2') +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
            coord_fixed() +
            coord_flip() + 
            scale_color_gradientn(colors = rdylbu_colors, 
                                  limits = c(-3, 3),   
                                  breaks = color_breaks) +         
            scale_size_area(limits = c(0, 100), oob = scales::squish) 


mCRC_sn_Fib = mCRC_sn_simplified_cell %>% 
                 subset(broad_cell_type %in% c('Fibroblast')) %>% 
                 subset(cell_type_xenium != 'Pericyte')

Idents(mCRC_sn_Fib) <- 'cell_type_xenium'

p2 = VlnPlot(
  mCRC_sn_Fib, 
  features = c('SLIT2'), 
  group.by = 'cell_type_xenium',
  split.by = 'Site_of_Origin',  
  idents = 'stromal_fibroblast',
  pt.size = 0.1   # no single-cell dots
) + 
  stat_summary(
    fun = median, 
    geom = "point", 
    size = 1.5, 
    color = "black"
  ) + 
  scale_fill_manual(
    values = c("liver" = "brown", "lung" = "steelblue1")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p3 = VlnPlot(
  mCRC_sn_Fib, 
  features = c('HGF'), 
  group.by = 'cell_type_xenium',
  split.by = 'Site_of_Origin',  
  idents = c('stromal_fibroblast', 'stellate_cell'),
  pt.size = 0.1   # no single-cell dots
) + 
  stat_summary(
    fun = median, 
    geom = "point", 
    size = 1.5, 
    color = "black"
  ) + 
  scale_fill_manual(
    values = c("liver" = "brown", "lung" = "steelblue1")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p4 = VlnPlot(
  mCRC_sn_Fib, 
  features = c('SLIT2'), 
  group.by = 'cell_type_xenium',
  pt.size = 0   # no single-cell dots
) + 
  stat_summary(
    fun = median, 
    geom = "point", 
    size = 1.5, 
    color = "black"
  ) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p5 = VlnPlot(
  mCRC_sn_Fib, 
  features = c('HGF'), 
  group.by = 'cell_type_xenium',
  pt.size = 0   # no single-cell dots
) + 
  stat_summary(
    fun = median, 
    geom = "point", 
    size = 1.5, 
    color = "black"
  ) + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(output_dir, "sn_cellchat_Ligand_expression_Fib_celltypes_violin_plot.pdf"), width = 12, height = 4)
patchwork::wrap_plots(p1, p4, p5, nrow = 1, byrow = TRUE)
dev.off()

pdf(file.path(output_dir, "sn_cellchat_Ligand_expression_Fib_SLIT2_HGF_violin_plot.pdf"), width = 5, height = 5)
patchwork::wrap_plots(p2,p3, nrow = 1, byrow = TRUE)
dev.off()