library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)

out_dir = getwd()
Tobj = readRDS('/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_5K/subset_objects/mCRC_Xenium_N26_stroma_harmonized_T_NK_DB_removed.rds')

Tobj <- Tobj %>% 
        subset(TNK_cell_type1 %in% c('NK', 'CD4_Treg', 'CD4_Th17_T', 'CD8_Effector_T', 
                                     'CD4_naive/CM_T', 'CD4_Resting_T', 'Proliferative_T', 'CD4_Th1_T', 'CD8_TRM_T'))

Tobj

tcell.color <- c('CD4_naive/CM_T' = 'seagreen1', 
                 'CD4_Resting_T' = 'yellowgreen',  
                 'CD4_Treg' = 'seagreen4', 
                 'CD4_Th17_T' = 'seagreen3', 
                 'CD8_TRM_T' = 'steelblue1', 
                 'CD8_EM_T'= 'cyan4', 
                 'CD8_Effector_T' = 'steelblue3', 
                 'CD4_Th1_T' = 'steelblue4', 
                 'Proliferative_T' = 'turquoise', 
                 'NK' = 'turquoise1')

p3 <- DimPlot(Tobj, group.by = "TNK_cell_type1", reduction = "umap.harmony", label = FALSE, raster=TRUE, raster.dpi = c(250, 250)) +
      scale_color_manual(values=tcell.color, name='cell type') 
p3


pdf(file.path(out_dir, 'Xenium_TNK_cells_umap.pdf'), width = 8, height = 6)
p3
dev.off()

Tobj$TNK_cell_type1 <- factor(Tobj$TNK_cell_type1,
                              level = rev(c('CD4_Resting_T', 'CD4_naive/CM_T', 'CD4_Treg', 'CD4_Th1_T', 'CD4_Th17_T', 
                                            'CD8_TRM_T', 'CD8_Effector_T', 'CD8_EM_T', 'NK', 'Proliferative_T'
                                           )))

markers = c('CD4', 'TCF7', 'SELL', 'FOXP3', 'CTLA4', 'TIGIT', 
            'TNF', 'IL17A', 'CCR6', 'CD8A', 'ITGA1', 'IFNG', 
            'GZMK', 'NCAM1', 'GZMB', 'MKI67')

DefaultAssay(Tobj) <- 'Xenium'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))

p4 <- DotPlot(Tobj, 
              group.by = 'TNK_cell_type1', 
              feature = markers) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors) +       
      scale_size_area(limits = c(0, 100), oob = scales::squish)

pdf(file.path(out_dir, 'Xenium_TNK_cells_dotplot.pdf'), width = 16, height = 4)
p4
dev.off()
