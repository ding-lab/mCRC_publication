library(Seurat)
library(tidyverse)
library(data.table)
library(RColorBrewer)

rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean6.rds'
cell_type_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean6_metadata_cell_type_all_20251112.csv'
out_dir = getwd()

obj = readRDS(rds_path)
cell_type = read.table(cell_type_file, header = TRUE, sep=',', row.name=1)
obj <- AddMetaData(obj, cell_type, col.name = c('broad_cell_type2'))

obj$broad_cell_type2 <- factor(obj$broad_cell_type2,
                               level = c('Tumor', 
                                        'Intestinal_epithelium',
                                        'Mesenchymal',
                                        'Myelocyte',
                                        'T_NK_cell',
                                        'B cell',
                                        'Plasma',
                                        'Vascular endothelial',
                                        'Lymphatic endothelial',
                                        'LSEC',
                                        'Hepatocyte',
                                        'Cholangiocyte',
                                        'Alveolar type1',
                                        'Alveolar type2',
                                        'Club',
                                        'Ciliated',
                                        'Adrenal cortex', 
                                        'Oligodendrocyte',
                                        'Astrocyte',
                                        'Adipocyte'
                                       ))

cell.color <- c('Tumor' = 'violetred', 
                 'Intestinal_epithelium' = 'gold2',
                 'Mesenchymal' = 'olivedrab4',
                 'Myelocyte' = 'violet',
                 'T_NK_cell' = 'steelblue3',
                 'B cell' = 'lightpink2', 
                 'Plasma' = 'lightpink4',         
                 'Vascular endothelial' = 'red2', 
                 'Lymphatic endothelial' = 'red3', 
                 'LSEC' = 'red4', 
                 'Hepatocyte' = 'sienna4', 
                 'Cholangiocyte' = 'chocolate3',
                 'Club' = 'thistle2', 
                 'Ciliated' = 'thistle3', 
                 'Alveolar type1' = 'indianred', 
                 'Alveolar type2' = 'thistle4', 
                 'Adrenal cortex' = 'hotpink4', 
                 'Oligodendrocyte' = 'burlywood3', 
                 'Astrocyte' = 'wheat4', 
                 'Adipocyte' = 'peachpuff')


all_markers <- c(
  "EPCAM",      # epithelial / tumor / intestinal / cholangio
  "COL1A1",     # mesenchymal / fibroblast
  "PTPRC",      # pan-immune  
  "LST1",       # myeloid
  "CD3D",       # T/NK
  "MS4A1",      # B cell
  "JCHAIN",     # plasma cell
  "PECAM1",     # vascular endothelial
  "PROX1",      # lymphatic endothelial
  "CLEC4G",     # LSEC
  "ALB",        # hepatocyte
  "KRT19",      # cholangiocyte
  "AGER",       # alveolar type 1
  "SFTPC",      # alveolar type 2
  "SCGB1A1",    # club cell
  "FOXJ1",      # ciliated epithelium
  "CYP11B1",    # adrenal cortex
  "PLP1",  
  "GFAP",       # astrocyte
  "ADIPOQ"      # adipocyte
)

p1 = DimPlot(obj, reduction = 'mCRCv6_umap.scvi', group.by = 'broad_cell_type2', label = FALSE) + 
     scale_color_manual(values=cell.color, name='cell type')



pdf(file.path(out_dir, 'snRNAseq_simple_cells_umap.pdf'), width = 8, height = 6)
p1
dev.off()

obj$broad_cell_type2 <- factor(obj$broad_cell_type2,
                               level = rev(c('Tumor', 
                                        'Intestinal_epithelium',
                                        'Mesenchymal',
                                        'Myelocyte',
                                        'T_NK_cell',
                                        'B cell',
                                        'Plasma',
                                        'Vascular endothelial',
                                        'Lymphatic endothelial',
                                        'LSEC',
                                        'Hepatocyte',
                                        'Cholangiocyte',
                                        'Alveolar type1',
                                        'Alveolar type2',
                                        'Club',
                                        'Ciliated',
                                        'Adrenal cortex', 
                                        'Oligodendrocyte',
                                        'Astrocyte',
                                        'Adipocyte'
                                       )))


DefaultAssay(obj) <- 'SCT'
rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))

p2 <- DotPlot(obj, 
              group.by = 'broad_cell_type2', 
              feature = all_markers) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors) +       
      scale_size_area(limits = c(0, 100), oob = scales::squish)


pdf(file.path(out_dir, 'snRNAseq_simple_cells_dotplot.pdf'), width = 16, height = 6)
p2
dev.off()