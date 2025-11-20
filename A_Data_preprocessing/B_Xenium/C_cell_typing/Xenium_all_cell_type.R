library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)
library(harmony)

out_dir = getwd()
rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_all_cells_noDB_normalized_noFOVs_sketch2_harmonized.qs'

obj = qread(rds_path)
obj

obj$All_cell_type3 <- factor(obj$All_cell_type3,
                             level = rev(c('Tumor', 'Intestinal_epithelium', 'Mesenchymal', 
                                           'Myelocyte', 'T/NK cell', 'B_cell', 'Plasma_cell',  'Endothelial', 'LSEC', 
                                           'Hepatocyte', 'Cholangiocyte', 'Breast_duct', 'Aveolar_epithelium', 'Bronchial_gland',
                                           'Enteric_neuron')))

cell.color <- c('Tumor' = 'violetred', 
                'Intestinal_epithelium' = 'gold2', 
                'Mesenchymal' = 'olivedrab4',
                'T/NK cell' = 'steelblue3', 
                'B_cell' = 'lightpink2', 
                'Plasma_cell' = 'lightpink4',    
                'Myelocyte' = 'violet',
                'Endothelial' = 'red2', 
                'LSEC' = 'red4', 
                'Hepatocyte' = 'sienna4', 
                'Cholangiocyte' = 'chocolate3',
                'Bronchial_gland' = 'thistle', 
                'Aveolar_epithelium' = 'thistle4', 
                'Breast_duct' = 'tan1', 
                'Enteric_neuron' = 'burlywood3')

p1 <- DimPlot(obj, group.by = "All_cell_type3", reduction = "umap", label = FALSE) + 
      scale_color_manual(values=cell.color, name='cell type')
p1

pdf(file.path(out_dir, 'Xenium_simple_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

all_markers <- c(
  "EPCAM",      # epithelial / tumor / intestinal / cholangio
  "VCAN",     # mesenchymal / fibroblast
  "PTPRC",      # pan-immune  
  "CD68",       # myeloid
  "CD3E",       # T/NK
  "MS4A1",      # B cell
  "CD38",     # plasma cell
  "PECAM1",     # vascular endothelial
  "CD34",     # LSEC
  "CYP3A4",        # hepatocyte
  "KRT7",      # cholangiocyte
  "AGER",       # alveolar type 1
  "SCGB1A1",    # bronhial gland
  "L1CAM"      # adipocyte
)


DefaultAssay(obj) <- 'Xenium'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))

p2 <- DotPlot(obj, 
              group.by = 'All_cell_type3', 
              feature = all_markers) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) +
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors) +       
      scale_size_area(limits = c(0, 100), oob = scales::squish)

pdf(file.path(out_dir, 'Xenium_simple_cells_dotplot.pdf'), width = 16, height = 6)
p2
dev.off()