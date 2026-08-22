library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(forcats)
library(scales)
library(rstatix)
library(Seurat)
library(qs)
library(RColorBrewer)

output_dir = getwd()

mye_rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/avisani/Xenium/mCRC_LM_5K_objects/'
Xenium_mye_obj = qread(file.path(mye_rds_path, 'mCRC_Xenium_N26_stroma_harmonized_myeloid_harmonized_joined.qs'))

metadata_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
meta_df <- read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"
meta_df <- meta_df %>% column_to_rownames("barcode")
Xenium_mye_obj <- AddMetaData(Xenium_mye_obj, meta_df)

Xenium_mye_obj <- Xenium_mye_obj %>% 
                  subset(All_cell_type1 %in% c('classical_monocyte', 'RTM_KC', 'M2_TAM', 'TAM',
                                               'PMN-like_monocyte', 'RTM_alveolar', 'IFN_TAM',
                                               'SPP1_TAM'
                                              ))%>% 
                  subset(Organ %in% c('Liver', 'Lung'))

Xenium_mye_obj$All_cell_type1 <- factor(Xenium_mye_obj$All_cell_type1,
                                     level = c(
                                         'PMN-like_monocyte', 'SPP1_TAM', 'IFN_TAM',
                                         'RTM_KC', 'M2_TAM',  'classical_monocyte', 'RTM_alveolar',  'TAM'
                                           ))

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

p = DotPlot(Xenium_mye_obj, features=c('HBEGF') , 
            group.by='All_cell_type1') +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      coord_flip() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 30), oob = scales::squish)

pdf(file.path(output_dir, "Xenium_cellchat_macrophage_HBEGF_expressions_dotplot.pdf"), width = 6, height = 3)
p
dev.off()
