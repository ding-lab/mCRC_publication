library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)

out_dir = getwd()
MEobj = qread('/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_ME_cells_noDB_normalized_noFOVs_sketch.qs')

mecell.color <- c('WNT5A_BMP' = 'sienna3', 
                  'WNT5A_infl' = 'orange1', 
                  'mCAF' = 'palegreen3', 
                  'iCAF' = 'darkseagreen4', 
                  'stromal_fibroblast' = 'olivedrab1',
                  'pericyte' = 'peru', 
                  'smooth_muscle' = 'lightseagreen',
                  'LSEC' = 'tomato',
                  'EC' = 'thistle', 
                  'Angiogenic_EC' = 'lightpink2',
                  'Tip_EC' = 'lightcoral',
                  'Activated_EC' = 'indianred1',
                  'Arterial_EC' = 'indianred3',
                  'Lymphatic_EC' = 'pink1',
                  'Alveolar_EC' = 'plum3'
                 )

p1 <- DimPlot(MEobj, group.by = "All_cell_type1", reduction = "umap", label = TRUE) +
      scale_color_manual(values=mecell.color, name='cell type') 


pdf(file.path(out_dir, 'Xenium_ME_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

Mobj <- MEobj %>% 
        subset(All_cell_type1 %in% c('WNT5A_BMP','WNT5A_infl','mCAF','stromal_fibroblast', 'iCAF','pericyte','smooth_muscle'))

M_markers = c('WNT2', 'BMP4', 'NRG1', 'WNT5A', 
              'VEGFA', 'POSTN',
              'PDGFRA', 'IL6', 'CCL2', 'C7', 
              'PDGFRB', 'RGS5', 'MYH11')

DefaultAssay(Mobj) <- 'Xenium'
p2 = DotPlot(Mobj, features=M_markers , group.by='All_cell_type1') +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = c(-2, -1, 0, 1, 2)) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

Eobj <- MEobj %>% 
        subset(All_cell_type1 %in% c('EC', 'Activated_EC', 'Arterial_EC', 'Alveolar_EC', 'Angiogenic_EC',
                                     'Tip_EC', 'Lymphatic_EC', 'LSEC'))

Eobj$All_cell_type1 <- factor(Eobj$All_cell_type1, 
                               level =rev(c('EC',
                                            'Angiogenic_EC',
                                            'Tip_EC',
                                            'Activated_EC',
                                            'Arterial_EC',
                                            'Alveolar_EC',
                                            'Lymphatic_EC',
                                            'LSEC')))
DefaultAssay(Eobj) <- 'Xenium'
p3 = DotPlot(Eobj, features=E_markers , group.by='All_cell_type1') +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = c(-2, -1, 0, 1, 2)) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

pdf(file.path(out_dir, 'Xenium_MesenEndo_cells_dotplot.pdf'), width = 17, height = 4)
patchwork::wrap_plots(p1,p2, nrow = 1, byrow = TRUE)
dev.off()