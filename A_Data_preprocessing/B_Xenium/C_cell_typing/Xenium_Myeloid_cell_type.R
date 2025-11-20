library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)

out_dir = getwd()
MYobj = qread('/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_MY_cells_noDB_normalized_noFOVs_sketch.qs')

mycell.color <- c('classical_monocyte' = 'violet', 
                  'RTM_KC' = 'violetred4',
                  'RTM_alveolar' = 'violetred1', 
                  'IFN_TAM' = 'violetred',     
                  'M2_TAM' = 'indianred1', 
                  'SPP1_TAM' = 'salmon', 
                  'PMN-like_monocyte' = 'orchid3', 
                  'TAM' = 'orchid4', 
                  'cDC1' = 'tan', 
                  'cDC2' = 'khaki3', 
                  'pDC' = 'tan1', 
                  'LAMP3_DC' = 'tan3')

p1 <- DimPlot(MYobj, group.by = "All_cell_type1", reduction = "umap", label = TRUE) +
      scale_color_manual(values=mycell.color, name='cell type') 


pdf(file.path(out_dir, 'Xenium_MY_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

MYobj@meta.data <- MYobj@meta.data %>%              
                   mutate(All_cell_type1 =
                          case_when(
                              All_cell_type1 == 'classical_monocyte' ~  'Classical-TIM',
                              All_cell_type1 == 'RTM_KC' ~ 'Kupffer RTM-like',
                              All_cell_type1 == 'RTM_alveolar' ~ 'Alveolar RTM-like',
                              All_cell_type1 == 'LAMP3_DC' ~ 'LAMP3+DC',
                              TRUE ~ All_cell_type1))

MYobj$All_cell_type1 <- factor(MYobj$All_cell_type1, 
                               level =rev(c('TAM',
                                            'M2_TAM', 
                                            'SPP1_TAM',
                                            'IFN_TAM',
                                            'Classical-TIM',
                                            'PMN-like_monocyte',
                                            'Alveolar RTM-like',
                                            'Kupffer RTM-like',
                                            'cDC1',
                                            'cDC2',
                                            'LAMP3+DC',
                                            'pDC')))

markers = c('ITGAX', 'CD68', 'MSR1', 'CD163', 'CSF1R',
            'STAB1', 'SPP1', 'VEGFA', 'IFIT1', 
            'VCAN', 'FCN1', 'NAMPT', 'MARCO', 'TIMD4', 
            'CLEC9A', 'CD1C', 'LAMP3',
            'CLEC4C')

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

p2 = DotPlot(MYobj, features=markers , group.by='All_cell_type1') +
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + 
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

pdf(file.path(out_dir, 'Xenium_Myeloid_cells_dotplot.pdf'), width = 16, height = 4)
p2
dev.off()
