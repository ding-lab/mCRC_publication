library(Seurat)
library(tidyverse)
library(data.table)
library(RColorBrewer)

out_dir = getwd()
T_NK_rds_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/T/57_Integrated_normalized_mCRC_snRNA_noDB_v7_T_NK_clean4.rds'
Tobj = readRDS(T_NK_rds_path)
Tobj

cell_type_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean6_metadata_cell_type_all_20251112.csv'
cell_type = read.table(cell_type_file, header = TRUE, sep=',', row.name=1)
Tobj <- AddMetaData(Tobj, cell_type, col.name = c('cell_type_all2'))

tcell.color <- c('CD4_naive/CM_T' = 'seagreen1', 
                'CD4_Resting_T' = 'yellowgreen',  
                'CD4_Treg' = 'seagreen4', 
                'CD4_Th17_T' = 'seagreen3', 
                'CD8_TRM_T' = 'steelblue1', 
                'CD8_EM_T'= 'cyan4', 
                'CD8_Effector_T' = 'steelblue3', 
                'CD8_NK_like_T' = 'steelblue4', 
                'CD56_low_NK' = 'turquoise', 
                'CD56_high_NK' = 'turquoise1')

options(repr.plot.width = 10, repr.plot.height = 8)
p1 = DimPlot(Tobj, reduction = 'T_NK_umap.scvi', group.by = 'cell_type_all2', label = FALSE, raster = TRUE, raster.dpi = c(200, 200)) + 
     scale_color_manual(values=tcell.color, name='cell type') 

pdf(file.path(out_dir, 'snRNAseq_TNK_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()


markers = c('CD4', 'IL7R', 'SELL', 'TCF7', 'FOXP3', 'CTLA4', 'TIGIT', 
            'IL17A', 'CCR6', 'CD8A', 'ITAGAE', 'ITGA1', 'IFNG', 'GZMK',
            'GNLY', 'NCAM1', 'GZMB')

Tobj$cell_type_all2 <- factor(Tobj$cell_type_all2,
                              level = rev(c('CD4_Resting_T', 'CD4_naive/CM_T', 'CD4_Treg', 'CD4_Th17_T', 
                                            'CD8_TRM_T', 'CD8_Effector_T', 'CD8_EM_T', 'CD8_NK_like_T',
                                            'CD56_high_NK', 'CD56_low_NK'
                                           )))

DefaultAssay(Tobj) <- 'SCT'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))

p2 <- DotPlot(Tobj, 
              group.by = 'cell_type_all2', 
              feature = markers) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors) +       
      scale_size_area(limits = c(0, 100), oob = scales::squish)


pdf(file.path(out_dir, 'snRNAseq_TNK_cells_dotplot.pdf'), width = 16, height = 4)
p2
dev.off()