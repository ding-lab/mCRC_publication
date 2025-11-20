library(Seurat)
library(tidyverse)
library(data.table)
library(RColorBrewer)

out_dir = getwd()

ME_rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_MesenEndo_clean3_reINT.rds'
MEobj = readRDS(ME_rds_path)
MEobj

cell_type_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean6_metadata_cell_type_all_20251112.csv'
cell_type = read.table(cell_type_file, header = TRUE, sep=',', row.name=1)
MEobj <- AddMetaData(MEobj, cell_type, col.name = c('cell_type_all2'))

mecell.color <- c('WNT5A_BMP' = 'sienna3', 
                  'WNT5A_infl' = 'orange1', 
                  'mCAF' = 'palegreen3', 
                  'iCAF' = 'darkseagreen4', 
                  'stromal_fibroblast' = 'olivedrab1',
                  'stellate_cell' = 'pink3', 
                  'Pericyte' = 'peru', 
                  'smooth_muscle' = 'lightseagreen',
                  'Vascular endothelial' = 'thistle', 
                  'Lymphatic endothelial' = 'pink1', 
                  'LSEC' = 'tomato')

p1 = DimPlot(MEobj, reduction = 'mesen_endo_umap.scvi', group.by = 'cell_type_all2', label = FALSE, raster = NULL, raster.dpi = c(200, 200)) + 
     scale_color_manual(values=mecell.color, name='cell type') 


pdf(file.path(out_dir, 'snRNAseq_ME_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

markers = c('WNT2', 'BMP4', 'NRG1', 'WNT5A',
            'VEGFA', 'POSTN',
            'PDGFRA', 'IL6', 'C7', 'PLA2G2A',
            'PDGFRB', 'RGS5', 'MYH11',  'RELN',
            'PECAM1', 'VWF', 'PROX1', 'CLEC4G')

MEobj$cell_type_all2 <- factor(MEobj$cell_type_all2, 
                               level =rev(c('WNT5A_BMP',
                                            'WNT5A_infl',
                                            'mCAF',
                                            'stromal_fibroblast', 
                                            'iCAF',
                                            'Pericyte',
                                            'smooth_muscle',
                                            'stellate_cell',
                                            'Vascular endothelial',
                                            'Lymphatic endothelial',
                                            'LSEC')))

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)


p2 = DotPlot(MEobj, features=markers , group.by='cell_type_all2') +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)


pdf(file.path(out_dir, 'snRNAseq_MesenEndo_cells_dotplot.pdf'), width = 16, height = 4)
p2
dev.off()