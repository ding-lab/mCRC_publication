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

metadata_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
meta_df <- read_csv(metadata_path)
colnames(meta_df)[1] <- "barcode"
meta_df <- meta_df %>% column_to_rownames("barcode")

all_rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/'
Xenium_all_obj = qread(file.path(all_rds_path, 'mCRC_Xenium_N26_all_cells_DB_removed_normalized.qs'))
Xenium_all_obj <- AddMetaData(Xenium_all_obj, meta_df)

Xenium_all_obj$Broad_cell_type3 <- factor(Xenium_all_obj$Broad_cell_type3,
                                          level = c('Tumor', 'Fibroblast', 'Myeloid cells', 'Hepatocyte',
                                                    'Aveolar_epithelium', 'Endothelial_cell', 'B/Plasma',
                                                    'T_NK_cell'))

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

p1 = DotPlot(Xenium_all_obj, features=c('SLIT2', 'HGF', 'HBEGF', 'ROBO1', 'MET', 'EGFR', 'ERBB2'), 
            group.by='Broad_cell_type3') +
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      		coord_fixed() +
      		coord_flip() + 
      		scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      		scale_size_area(limits = c(0, 100), oob = scales::squish) 

output_dir = getwd()
pdf(file.path(output_dir, "Xenium_cellchat_LR_expression_dotplot.pdf"), width = 5, height = 5)
p1
dev.off()