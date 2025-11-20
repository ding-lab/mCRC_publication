library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)

out_dir = getwd()
cell_type_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean6_metadata_cell_type_all_20251112.csv'
cell_type = read.table(cell_type_file, header = TRUE, sep=',', row.name=1)

MY_rds_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/myeloid/AV_myeloid_subset4.rds_reINT.rds'
MYobj = readRDS(MY_rds_path)
MYobj <- AddMetaData(MYobj, cell_type, col.name = c('cell_type_all4'))

MYobj <- MYobj %>% 
         subset(cell_type_all4 != 'Doublet') 

mycell.color <- c('Classical-TIM' = 'violet', 
                  'Kupffer RTM-like' = 'violetred4',
                  'Alveolar RTM-like' = 'violetred1', 
                  'IFN-TAM' = 'violetred', 
                  'M2_TAM' = 'indianred1',
                  'SPP1_TAM' = 'salmon', 
                  'LA-TAM' = 'thistle2',    
                  'PMN-like_monocyte' = 'orchid3', 
                  'TAM' = 'orchid4', 
                  'cDC1' = 'tan', 
                  'cDC2' = 'khaki3',
                  'pDC' = 'tan1', 
                  'LAMP3+DC' = 'tan3', 
                  'mast' = 'tan4',
                  'Spatial_doublet' = 'gray'
                 )

p1 = DimPlot(MYobj, reduction = 'myeloid4_umap.scvi', group.by = 'cell_type_all4', label = TRUE) + 
     scale_color_manual(values=mycell.color , name='cell type') 

pdf(file.path(out_dir, 'snRNAseq_MY_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

markers = c('ITGAX', 'CD68', 'MSR1', 'CD163', 'CSF1R',
            'STAB1', 'SPP1', 'VEGFA', #'IFIT1', 
            'VCAN', 'FCN1', 'MARCO', 'TIMD4', #'NAMPT',
            'CLEC9A', 'CD1C', 'LAMP3',
            'CLEC4C', 'KIT')

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

MYobj$cell_type_all4 <- factor(MYobj$cell_type_all4, 
                               level =rev(c('TAM',
                                            'M2_TAM', 
                                            'SPP1_TAM',
                                            'Classical-TIM',
                                            'Alveolar RTM-like',
                                            'Kupffer RTM-like',
                                            'cDC1',
                                            'cDC2',
                                            'LAMP3+DC',
                                            'pDC',
                                            'mast')))

p2 = DotPlot(MYobj, features=markers , group.by='cell_type_all4') +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)


pdf(file.path(out_dir, 'snRNAseq_Myeloid_cells_dotplot.pdf'), width = 16, height = 4)
p2
dev.off()