library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(MAST)
library(readr)
library(RColorBrewer)
library(scales)
library(rstatix)
library(googlesheets4)
source("../../../scripts/jupyter_support_functions.R")

# setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
# epithelial_reint = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean3_reINT.rds')
# epithelial_reint

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial')
epithelial_reint = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds')
epithelial_reint

colnames(epithelial_reint@meta.data)

unique(epithelial_reint$Site_of_Origin)

head(epithelial_reint@meta.data, 3)

epithelial_reint

wxs_epithelial <- subset(epithelial_reint, 
                         subset = (orig.ident == 'CM268C1-S1' | orig.ident == 'CM268C1-T1' | orig.ident == 'CM354C1-T1' |
                                   orig.ident == 'CM354C2-T1' | orig.ident == 'CM268C1-T1' | orig.ident == 'CM392C1-S1' |
                                   orig.ident == 'CM392C2-Th1' | orig.ident == 'CM426C1-Th1' | orig.ident == 'CM426C2-Tp1' |
                                   orig.ident == 'CM478C1-T1Y2' | orig.ident == 'HT213C1-Te1' | orig.ident == 'HT225C1-Th1' |
                                   orig.ident == 'HT230C1-Th1' | orig.ident == 'HT253C1-Th1' | orig.ident == 'HT260C1-Th1' |
                                   orig.ident == 'HT266C1-Tb1'))

p <- DimPlot(wxs_epithelial, group.by = "KRAS_mut", reduction = "epithelial_umap.scvi", order = T, label=F, pt.size=1, label.size=6, raster=FALSE) +
    labs(title = "KRAS_alleles", color="Allele_Detected") + 
    scale_color_manual(labels = c("NA", "KRAS reference", "KRAS variant"), values = c("grey", "blue", "red"))
pdf('KRAS_alleles_umap.pdf', height = 6, width = 6)
p
dev.off()

set_size(6,6)
p

p <- DimPlot(wxs_epithelial, group.by = "APC_mut", reduction = "epithelial_umap.scvi", order = T, label=F, pt.size=1, label.size=6, raster=FALSE) +
    labs(title = "APC_alleles", color="Allele_Detected") + 
    scale_color_manual(labels = c("NA", "APC reference", "APC variant"), values = c("grey", "blue", "red"))
pdf('APC_alleles_umap.pdf', height = 6, width = 6)
p
dev.off()

set_size(6,6)
p

p <- DimPlot(wxs_epithelial, group.by = "TP53_mut", reduction = "epithelial_umap.scvi", order = T, label=F, pt.size=1, label.size=6, raster=FALSE) +
    labs(title = "TP53_alleles", color="Allele_Detected") + 
    scale_color_manual(labels = c("NA", "TP53 variant"), values = c("grey", "red"))
pdf('TP53_alleles_umap.pdf', height = 6, width = 6)
p
dev.off()

set_size(6,6)
p

colnames(wxs_epithelial@meta.data)

set_size(6,6)
p1 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub6', label = TRUE)
p1

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('MUC2', 'TFF1', 'SYTL2'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('KLK1', 'RETNLB', 'CLCA1'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('RGMB', 'SMOC2', 'LGR5', 'ASCL2', 'OLFM4'), reduction = 'epithelial_umap.scvi')

set_size(8,12)
FeaturePlot(epithelial_reint, features = c('SH2D6', 'TRPM5', 'BMX', 'LRMP', 'HCK'), reduction = 'epithelial_umap.scvi')

set_size(8,12)
FeaturePlot(epithelial_reint, features = c('SCGN', 'FEV', 'CHGA', 'CHGB', 'PYY', 'GCG'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('TICRR', 'CDC25C', 'MKI67', 'TOP2A'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('LYZ', 'DEFA5'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('SLC26A3', 'KRT20'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('SLC2A1', 'VEGFA'), reduction = 'epithelial_umap.scvi')

organ_col <- c(colon = 'sienna3', rectum = 'firebrick3', liver = 'brown', 
               lung = 'steelblue1', brain = 'bisque1', adrenal = 'cyan4', 
               breast = 'violet', spleen = 'lightcoral')
tissue_col <- c(primary = 'gold3', metastasis = 'maroon3', normal ='lightblue1')

set_size(12,12)
p2 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)
p3 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'Site_of_Origin', 
             label=FALSE) + scale_color_manual(values=organ_col, name='organ')
p4 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'Tissue_Type', 
             label = FALSE) + scale_color_manual(values=tissue_col, name='tissue')
p1 + p2  + p3 + p4

pdf("Dimplot_mCRC_epithelial_organ.pdf", width=6, height=6)
print(p3)
dev.off()

pdf("Dimplot_mCRC_epithelial_tissue.pdf", width=6, height=6)
print(p4)
dev.off()

set_size(10,6)
p5 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'orig.ident', label = TRUE)
p5

pdf("Dimplot_mCRC_epithelial_samples.pdf", width=12, height=6)
print(p5)
dev.off()

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'MSI', label = TRUE)

library(scCustomize)

Idents(epithelial_reint) = 'MSI'
MSI_cells <- WhichCells(object = epithelial_reint, idents = 'MSI_H')
cells <- list(MSI_cells = MSI_cells)
Cell_Highlight_Plot(seurat_object = epithelial_reint,
                    reduction = 'epithelial_umap.scvi',
                    cells_highlight = cells, raster = FALSE)

Idents(epithelial_reint) = 'orig.ident'
CM329C1_S1_cells <- WhichCells(object = epithelial_reint, idents = 'CM329C1-S1')
cells <- list(CM329C1_S1_cells = CM329C1_S1_cells)
Cell_Highlight_Plot(seurat_object = epithelial_reint,
                    reduction = 'epithelial_umap.scvi',
                    cells_highlight = cells, raster = FALSE)

Idents(epithelial_reint) = 'orig.ident'
CM655C1_S1_cells <- WhichCells(object = epithelial_reint, idents = 'CM655C1-S1')
cells <- list(CM655C1_S1_cells = CM655C1_S1_cells)
Cell_Highlight_Plot(seurat_object = epithelial_reint,
                    reduction = 'epithelial_umap.scvi',
                    cells_highlight = cells, raster = FALSE)

# MSI_H = CM329C1-S1, CM655C1-S1

table(epithelial_reint$epithelial_clusters_0.1)

epithelial_reint <- epithelial_reint %>% 
                    subset(epithelial_clusters_0.1 == 0 |
                           epithelial_clusters_0.1 == 1 |
                           epithelial_clusters_0.1 == 10 |
                           epithelial_clusters_0.1 == 2 |
                           epithelial_clusters_0.1 == 3 |
                           epithelial_clusters_0.1 == 4 
                          ) 

epithelial_reint <- FindClusters(epithelial_reint, 
                                 resolution = 0.1, 
                                 cluster.name = 'epi_clusters_0.1', 
                                 graph = 'RNA_snn')

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_clusters_0.1', label = TRUE)

table(epithelial_reint$epi_clusters_0.1)

epithelial_reint <- epithelial_reint %>% 
                    subset(epi_clusters_0.1 == 0 |
                           epi_clusters_0.1 == 1 |
                           epi_clusters_0.1 == 2 |
                           epi_clusters_0.1 == 3 |
                           epi_clusters_0.1 == 4 
                          ) 

Idents(epithelial_reint) = "epi_clusters_0.1"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 0,
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epithelial_cluster_sub0",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epithelial_cluster_sub0)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub0', label = TRUE)

DotPlot(epithelial_reint, feature = c('LGR5', 'VEGFA'), group.by = 'epithelial_cluster_sub0')

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                    mutate(epithelial_cluster_sub0 = case_when(epithelial_cluster_sub0 == '0_3' ~ '0_1',
                                                               TRUE ~ epithelial_cluster_sub0
                                                              ))

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub0', label = TRUE)

Idents(epithelial_reint) = "epithelial_cluster_sub0"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 2,
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epithelial_cluster_sub2",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epithelial_cluster_sub2)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub2', label = TRUE)

Idents(epithelial_reint) = "epithelial_cluster_sub2"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = '2_2',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epithelial_cluster_sub2",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epithelial_cluster_sub2)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub2', label = TRUE)

DotPlot(epithelial_reint, feature = c('SH2D6', 'TRPM5', 'BMX', 'LRMP', 'HCK', 'CHGA', 'CHGB'), 
        group.by = 'epithelial_cluster_sub2')

Highlight_Cluster_UMAP(seurat_object = epithelial_reint, 
                       metadata_column = 'epithelial_cluster_sub3', 
                       cell_type ='2_3', 
                       umap_reduction = 'epithelial_umap.scvi', 
                       highlight_color = "red")

Idents(epithelial_reint) = "epithelial_cluster_sub2"
epithelial_2_3_markers = FindMarkers(epithelial_reint,
                                     ident.1 = '2_3',
                                     min.pct=0.1, 
                                     logfc.threshold=0.5, 
                                     only.pos = TRUE, 
                                     return.thresh = 0.01,
                                     test.use = 'MAST')
write.csv(epithelial_2_3_markers,  "mCRC_epithelial_clean3_reINT_cluster2_3_markers_mast_deg.csv", row.names = TRUE)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                    mutate(epithelial_cluster_sub2 = case_when(epithelial_cluster_sub2 == '2_1' ~ '2_0',
                                                               epithelial_cluster_sub2 == '2_3' ~ '2_0',
                                                               epithelial_cluster_sub2 == '2_2_2' ~ '2_2_0',
                                                               TRUE ~ epithelial_cluster_sub2
                                                              ))

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub2', label = TRUE)

Idents(epithelial_reint) = "epithelial_cluster_sub2"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 3,
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epithelial_cluster_sub3",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epithelial_cluster_sub3)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub3', label = TRUE)

Idents(epithelial_reint) = "epithelial_cluster_sub3"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = '3_1',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epithelial_cluster_sub3",
                                  resolution = 0.2,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epithelial_cluster_sub3)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub3', label = TRUE)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                    mutate(epithelial_cluster_sub3 = case_when(epithelial_cluster_sub3 == '3_1_1' ~ '3_1_0',
                                                               epithelial_cluster_sub3 == '3_1_3' ~ '3_1_0',
                                                               TRUE ~ epithelial_cluster_sub3
                                                              ))

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub3', label = TRUE)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub3', label = TRUE)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type = case_when(
                                  epithelial_cluster_sub3 == '0_0' ~ 'Stem-like tumor',
                                  epithelial_cluster_sub3 == '0_1' ~ 'Angiogenic tumor',
                                  epithelial_cluster_sub3 == '0_2' ~ 'Proliferative stem-like tumor',
                                  epithelial_cluster_sub3 == '1' ~ 'Proliferative tumor',
                                  epithelial_cluster_sub3 == '2_0' ~ 'APCDD1+ tumor',
                                  epithelial_cluster_sub3 == '2_2_0' ~ 'Enteroendocrine-like cells',
                                  epithelial_cluster_sub3 == '2_2_1' ~ 'Tuft cells',
                                  epithelial_cluster_sub3 == '3_0' ~ 'Enterocytes',
                                  epithelial_cluster_sub3 == '3_1_0' ~ 'Transit-amplifying cells',
                                  epithelial_cluster_sub3 == '3_1_2' ~ 'Stem cells',
                                  epithelial_cluster_sub3 == '4' ~ 'Goblet cells'
                              ))

epithelial_reint@meta.data %>% select(epi_cell_type) %>% 
                               write.csv('mCRC_epithelial_metadata.csv', row.names = TRUE)

saveRDS(epithelial_reint, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4.rds')

epithelial_reint$epi_cell_type <- factor(epithelial_reint$epi_cell_type, 
                                         levels = c('APCDD1+ tumor', 'Angiogenic tumor', 'Proliferative tumor',
                                                    'Proliferative stem-like tumor', 'Stem-like tumor',
                                                    'Stem cells', 'Enterocytes', 'Transit-amplifying cells', 
                                                    'Goblet cells', 'Tuft cells', 'Enteroendocrine-like cells'
                                                   ))

epi.cell.color <- c("APCDD1+ tumor" = '#008856',
                    "Angiogenic tumor" = '#be0032',
                    "Proliferative tumor" = '#f3c300',
                    "Proliferative stem-like tumor" = '#8db600',
                    "Stem-like tumor" = '#e68fac',
                    "Stem cells" = '#e25822',
                    "Enterocytes" = '#882d17',
                    "Transit-amplifying cells" = '#a1caf1',
                    "Goblet cells" = '#c2b280',
                    "Tuft cells" = '#f99379',
                    "Enteroendocrine-like cells" = '#2b3d26'
                   )

set_size(7,6)
p <- DimPlot(epithelial_reint, 
             reduction = 'epithelial_umap.scvi',  
             group.by = 'epi_cell_type', 
             label=FALSE) + scale_color_manual(values=epi.cell.color, name='cell type')
p

pdf("Dimplot_mCRC_epithelial_cell_type.pdf", width=7, height=6)
print(p)
dev.off()

Idents(epithelial_reint) = "epi_cell_type"
epithelial_markers = FindAllMarkers(epithelial_reint, 
                                    min.pct=0.1, 
                                    logfc.threshold=0.5, 
                                    only.pos = TRUE, 
                                    return.thresh = 0.01,
                                    test.use = 'MAST'
                                    )

write.csv(epithelial_markers, "mCRC_epithelial_mast_deg.csv", row.names = TRUE)

epithelial_reint$epi_cell_type <- factor(epithelial_reint$epi_cell_type, 
                                         levels = rev(c('APCDD1+ tumor', 'Angiogenic tumor', 'Proliferative tumor',
                                                    'Proliferative stem-like tumor', 'Stem-like tumor',
                                                    'Stem cells', 'Enterocytes', 'Transit-amplifying cells', 
                                                    'Goblet cells', 'Tuft cells', 'Enteroendocrine-like cells'
                                                   )))

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

pdf("Dotplot_mCRC_epithelial_cell_type_RNA.pdf", width=8, height=8)
print(p6)
dev.off()

DefaultAssay(epithelial_reint) <- 'SCT'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-2.5, 2.5),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

pdf("Dotplot_mCRC_epithelial_cell_type_SCT.pdf", width=8, height=8)
print(p6)
dev.off()

set_size(6,4)
DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type', 
              feature = c('CCND1', 'EEF2', 'EPHA2', 'ERBB2', 'FRK', 'GART', 'HDAC2', 'HDAC3',
                          'HPRT1', 'IDH1', 'IDH2', 'JUN', 'KRAS', 'MET', 'PARP1', 'PSMB1', 'RPL3', 
                          'TOP1', 'TOP1MT', 'TUBB', 'TXNRD1', 'VEGFA')) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors) +       
      scale_size_area(limits = c(0, 100), oob = scales::squish)

Idents(epithelial_reint) = "epi_cell_type"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 'Angiogenic tumor',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epi_cell_type2",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epi_cell_type2)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type2', label = TRUE)

table(epithelial_reint$epi_cell_type2)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                    mutate(epi_cell_type2 = case_when(epi_cell_type2 == 'Angiogenic tumor_4' ~ 'Angiogenic tumor_1',
                                                      epi_cell_type2 == 'Angiogenic tumor_5' ~ 'Angiogenic tumor_1',
                                                      epi_cell_type2 == 'Angiogenic tumor_2' ~ 'Angiogenic tumor_0',
                                                      epi_cell_type2 == 'Angiogenic tumor_3' ~ 'Angiogenic tumor_0',
                                                      epi_cell_type2 == 'Angiogenic tumor_6' ~ 'Angiogenic tumor_0',
                                                      TRUE ~ epi_cell_type2
                                                              ))

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type2', label = TRUE)

DotPlot(epithelial_reint, feature = c('SLC2A1'), 
        group.by = 'epi_cell_type2')

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/snRNA/Epithelial')

colnames(epithelial_reint@meta.data)

epi_metadata = epithelial_reint@meta.data %>% select('orig.ident', 'Patient_ID', 'Site_of_Origin', 
                                                     'Tissue_Type', 'Primary_Side', 'MSI', 'Tx_in_6mo', 'epi_cell_type')

tumor_metadata = epi_metadata %>% filter(epi_cell_type == 'APCDD1+ tumor' |
                                         epi_cell_type == 'Angiogenic tumor' |
                                         epi_cell_type == 'Proliferative tumor' |
                                         epi_cell_type == 'Proliferative stem-like tumor' |
                                         epi_cell_type == 'Stem-like tumor'
                                        ) %>% 
                                  filter(Tissue_Type != 'normal')

tumor_tbl <- tumor_metadata %>%
               group_by(orig.ident) %>% 
               summarise(cell_count = n(),
                         APCDD1_tumor_count = sum(epi_cell_type == 'APCDD1+ tumor', na.rm = TRUE),
                         APCDD1_tumor_proportion = round(100*APCDD1_tumor_count/n(),2),
                         Angiogenic_tumor_count = sum(epi_cell_type == 'Angiogenic tumor', na.rm = TRUE),
                         Angiogenic_tumor_proportion = round(100*Angiogenic_tumor_count/n(),2),
                         Proliferative_tumor_count = sum(epi_cell_type == 'Proliferative tumor', na.rm = TRUE),
                         Proliferative_tumor_proportion = round(100*Proliferative_tumor_count/n(),2),
                         Proliferative_stem_tumor_count = sum(epi_cell_type == 'Proliferative stem-like tumor', na.rm = TRUE),
                         Proliferative_stem_tumor_proportion = round(100*Proliferative_stem_tumor_count/n(),2),
                         Stem_like_tumor_count = sum(epi_cell_type == 'Stem-like tumor', na.rm = TRUE),
                         Stem_like_tumor_proportion = round(100*Stem_like_tumor_count/n(),2),
                         tissue_type = dplyr::first(Tissue_Type),
                         Organ = dplyr::first(Site_of_Origin),
                         Patient_ID = dplyr::first(Patient_ID),
                         Primary_Side = dplyr::first(Primary_Side),
                         MSI = dplyr::first(MSI),
                         Tx_in_6mo = dplyr::first(Tx_in_6mo),
                         .groups = 'drop')

write.table(tumor_tbl, 'mCRC_tumor_proportion_tbl.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
tumor_tbl = tumor_tbl %>% filter(cell_count > 100)
write.table(tumor_tbl, 'mCRC_tumor_proportion_tbl_filter100.tsv', sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

plot_clinical_clusters <- function(data, columns, var, ncol = 1, nrow = 1, palette = c("#00AFBB", "#FC4E07"), show_pvalues = TRUE) {
  plot_list <- list()  
  
  for (col in columns) {
    formula <- paste(col, "~", var)
    
    if (show_pvalues) {
      stat <- data %>% wilcox_test(as.formula(formula)) %>% add_xy_position(x = var)
    }
    
    plot <- ggboxplot(data, var, col,
                      color = var, palette = palette,
                      add = "jitter")
    
    if (show_pvalues) {
      plot <- plot + stat_pvalue_manual(stat, label = "p", tip.length = 0.01)
    }
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  # Combine all plots into one figure
  combined_plot <- ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)
  
  return(combined_plot)
}

columns_to_plot <- c("APCDD1_tumor_proportion", "Angiogenic_tumor_proportion", "Proliferative_tumor_proportion",
                     "Proliferative_stem_tumor_proportion", "Stem_like_tumor_proportion"
                    )
combined_plot <- plot_clinical_clusters(data = tumor_tbl, 
                                        columns = columns_to_plot, 
                                        var = "tissue_type", 
                                        palette = c('gold3', 'maroon3'),
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_tumor_proportion_primary_vs_metastasis.pdf", width=12, height=4)
combined_plot
dev.off()

unique(tumor_tbl$Organ)


combined_plot <- plot_clinical_clusters(data = tumor_tbl, 
                                        columns = columns_to_plot, 
                                        var = "Organ", 
                                        palette = c('sienna3', 'brown', 'steelblue1',
                                                    'firebrick3', 'bisque1', 'cyan4',
                                                    'violet', 'lightcoral'
                                                   ),
                                        show_pvalues = FALSE,
                                        ncol = 2, nrow = 3)

set_size(8, 12)
print(combined_plot)

pdf("Boxplot_mCRC_tumor_proportion_organs.pdf", width=10, height=12)
combined_plot
dev.off()

unique(tumor_tbl$Primary_Side)

combined_plot <- plot_clinical_clusters(data = tumor_tbl, 
                                        columns = columns_to_plot, 
                                        var = "Primary_Side", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_tumor_proportion_Primary_Side.pdf", width=12, height=4)
combined_plot
dev.off()

liver_tbl = tumor_tbl %>% filter(Organ == 'liver')

primary_tbl = tumor_tbl %>% filter(tissue_type == 'primary')

combined_plot <- plot_clinical_clusters(data = liver_tbl, 
                                        columns = columns_to_plot, 
                                        var = "Primary_Side", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_liver_tumor_proportion_Primary_Side.pdf", width=12, height=4)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "Primary_Side", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_Primary_Side.pdf", width=12, height=4)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = tumor_tbl, 
                                        columns = columns_to_plot, 
                                        var = "MSI", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_tumor_proportion_MSI.pdf", width=12, height=4)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = liver_tbl, 
                                        columns = columns_to_plot, 
                                        var = "MSI", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = FALSE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_liver_tumor_proportion_MSI.pdf", width=12, height=4)
combined_plot
dev.off()

table(primary_tbl$MSI, primary_tbl$Primary_Side)

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "MSI", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_MSI.pdf", width=12, height=4)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = tumor_tbl, 
                                        columns = columns_to_plot, 
                                        var = "Tx_in_6mo", 
                                        palette = c('sienna3', 'brown', 'steelblue1'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_tumor_proportion_Tx_in_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

iri_6mo <- c("CM1563C1-T1Y1", "CM1563C1-S1Y1", "CM354C1-T1", "CM354C2-T1", "CM426C1-Th1", "CM426C2-Tp1", "CM492C1-S1", "CM492C2-T1", "CM492C2-T2", "CM556C1-T1", "CM556C1-T2", 
	"CM556C2-T1", "CM655C1-S1", "CM663C1-T1", "CM873C1-S1", "HT225C1-Th1", "HT253C1-Th1", "HT260C1-Th1", "HT291C1-M1", "HT525C1-Th1", "HT539C1-Th1")
ox_6mo <- c("CM1563C1-T1Y1", "CM1563C1-S1Y1", "CM354C1-T1", "CM354C2-T1", "CM478C1-T1Y2", "CM492C1-S1", "CM492C2-T1", "CM492C2-T2", "CM556C1-T1", "CM556C1-T2", "CM556C2-T1", 
	"CM618C2-T1", "CM655C1-S1", "CM663C1-T1", "CM724C1-S1", "CM743C1-S1", "HT230C1-Th1", "HT291C1-M1", "HT307C1-Th1", "HT472C1-Th1", "HT472C1-S1", "HT525C1-Th1", "HT539C1-Th1")
bev_6mo <- c("CM1563C1-T1Y1", "CM1563C1-S1Y1", "CM354C1-T1", "CM354C2-T1", "CM426C1-Th1", "CM426C2-Tp1", "CM492C1-S1", 
	"CM492C2-T1", "CM492C2-T2", "CM556C1-T1", "CM556C1-T2", "CM556C2-T1", "CM655C1-S1", "CM663C1-T1", "HT225C1-Th1",
	"HT260C1-Th1", "HT291C1-M1", "HT307C1-Th1", "HT539C1-Th1")
egfr_6mo <- c("CM873C1-S1")
icb_6mo <- c("HT307C1-Th1", "HT525C1-Th1")
rt_6mo <- c("CM743C1-S1", "HT525C1-Th1")

tumor_tbl$IRI_6mo <- ifelse(tumor_tbl$orig.ident %in% iri_6mo, 'Yes', 'No')
tumor_tbl$OX_6mo <- ifelse(tumor_tbl$orig.ident %in% ox_6mo, 'Yes', 'No')
tumor_tbl$BEV_6mo <- ifelse(tumor_tbl$orig.ident %in% bev_6mo, 'Yes', 'No')
tumor_tbl$EGFR_6mo <- ifelse(tumor_tbl$orig.ident %in% egfr_6mo, 'Yes', 'No')
tumor_tbl$ICB_6mo <- ifelse(tumor_tbl$orig.ident %in% icb_6mo, 'Yes', 'No')
tumor_tbl$RT_6mo <- ifelse(tumor_tbl$orig.ident %in% rt_6mo, 'Yes', 'No')

colnames(tumor_tbl)

liver_tbl <- tumor_tbl %>% filter(Organ == 'liver')

primary_tbl <- tumor_tbl %>% filter(tissue_type == 'primary')

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "IRI_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_IRI_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "OX_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_OX_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "BEV_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_BEV_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "EGFR_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_EGFR_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = liver_tbl, 
                                        columns = columns_to_plot, 
                                        var = "ICB_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_liver_tumor_proportion_ICB_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

combined_plot <- plot_clinical_clusters(data = primary_tbl, 
                                        columns = columns_to_plot, 
                                        var = "RT_6mo", 
                                        palette = c('sienna3', 'brown'
                                                   ),
                                        show_pvalues = TRUE,
                                        ncol = 5, nrow = 1)

set_size(12, 4)
print(combined_plot)

pdf("Boxplot_mCRC_primary_tumor_proportion_RT_6mo.pdf", width=10, height=12)
combined_plot
dev.off()

gs4_deauth()
clinical_info.df <- read_sheet('https://docs.google.com/spreadsheets/d/1kKb-DC3jobabcp29VNkkrN024XcHBDavkiFEBxKsBTg/edit#gid=981909369', 
                    sheet = "Clinical")
clinical_info_selected <- clinical_info.df %>% dplyr::select('Patient ID', Dx_to_death_or_last_fu_days, 'Vital status')
head(clinical_info_selected, 1)

tumor_tbl_survival <- tumor_tbl %>% left_join(clinical_info_selected, by = c("Patient_ID" = "Patient ID"))
tumor_tbl_survival$Vital_status <- as.numeric(tumor_tbl_survival$'Vital status' == "dead")
head(tumor_tbl_survival, 10)

colnames(tumor_tbl_survival)

prepare_survival_data <- function(data, tissue_type_col = "tissue_type", tissue_type_value = "metastasis", 
                                    time_col_days = "Dx_to_death_or_last_fu_days", 
                                    proportions_cols = list("Angiogenic_tumor_proportion", 
                                                            "APCDD1_tumor_proportion", 
                                                            "Proliferative_tumor_proportion", 
                                                            "Proliferative_stem_tumor_proportion", 
                                                            "Stem_like_tumor_proportion")) {
  
  # Filter the data for the specified tissue type
  filtered_data <- data %>% filter(!!sym(tissue_type_col) == tissue_type_value)
  
  # Calculate the median values for each specified tumor proportion column
  medians <- lapply(proportions_cols, function(col) {
    median(filtered_data[[col]], na.rm = TRUE)
  })
  names(medians) <- proportions_cols
  
  # Create grouping variables based on the median values
  for (col in proportions_cols) {
    group_col_name <- paste0(gsub("_tumor_proportion", "", col), "_tumor_group")
    print(paste(col, 'median:', medians[[col]]))  
    filtered_data[[group_col_name]] <- ifelse(filtered_data[[col]] <= medians[[col]], "Low", "High")
  }
  
  # Convert time from days to months
  time_col_months <- gsub("_days", "_months", time_col_days)
  filtered_data[[time_col_months]] <- filtered_data[[time_col_days]] / 30
  
  # Filter out rows with NA in the time column
  result_data <- filtered_data %>% filter(!is.na(filtered_data[[time_col_days]]))
  
  return(result_data)
}

# Example usage:
# metastasis_sample_tbl <- prepare_metastasis_data(tumor_tbl_survival)


metastasis_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                               tissue_type_col = "tissue_type", 
                                               tissue_type_value = "metastasis",
                                               time_col_days = "Dx_to_death_or_last_fu_days"
                                              )

primary_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                               tissue_type_col = "tissue_type", 
                                               tissue_type_value = "primary",
                                               time_col_days = "Dx_to_death_or_last_fu_days"
                                              )

liver_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                          tissue_type_col = "Organ", 
                                          tissue_type_value = "liver",
                                          time_col_days = "Dx_to_death_or_last_fu_days"
                                         )

liver_sample_tbl <- liver_sample_tbl %>% mutate(Angiogenic_low_APCDD1 = if_else(Angiogenic_tumor_group == 'High' & 
                                                                                APCDD1_tumor_group == 'Low', 
                                                                                'High_angiogenic_Low_APCDD1',
                                                                                'Others'))

liver_sample_tbl <- liver_sample_tbl %>% 
                    mutate(Angiogenic_APCDD1_cat = 
                           case_when(Angiogenic_tumor_group == 'High' & APCDD1_tumor_group == 'Low' ~'High_angiogenic_Low_APCDD1',
                                     Angiogenic_tumor_group == 'High' & APCDD1_tumor_group == 'High' ~'High_angiogenic_High_APCDD1',
                                     Angiogenic_tumor_group == 'Low' & APCDD1_tumor_group == 'Low' ~'Low_angiogenic_Low_APCDD1',
                                     Angiogenic_tumor_group == 'Low' & APCDD1_tumor_group == 'High' ~'Low_angiogenic_High_APCDD1'))

liver_sample_tbl

lung_sample_tbl <- prepare_survival_data(tumor_tbl_survival,
                                          tissue_type_col = "Organ", 
                                          tissue_type_value = "lung",
                                          time_col_days = "Dx_to_death_or_last_fu_days"
                                         )

library(survival)
library(survminer)

# Prepare your survival object
surv_obj <- Surv(time = metastasis_sample_tbl$Dx_to_death_or_last_fu_months, event = metastasis_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Angiogenic_tumor_group, data = metastasis_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_tumor_group, data = metastasis_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = metastasis_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = metastasis_sample_tbl$Dx_to_death_or_last_fu_months, event = metastasis_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ APCDD1_tumor_group, data = metastasis_sample_tbl)

surv_diff <- survdiff(surv_obj ~ APCDD1_tumor_group, data = metastasis_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = metastasis_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = metastasis_sample_tbl$Dx_to_death_or_last_fu_months, event = metastasis_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Proliferative_tumor_group, data = metastasis_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_tumor_group, data = metastasis_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = metastasis_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = metastasis_sample_tbl$Dx_to_death_or_last_fu_months, event = metastasis_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Proliferative_stem_tumor_group, data = metastasis_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_stem_tumor_group, data = metastasis_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = metastasis_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = metastasis_sample_tbl$Dx_to_death_or_last_fu_months, event = metastasis_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Stem_like_tumor_group, data = metastasis_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Stem_like_tumor_group, data = metastasis_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = metastasis_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Angiogenic_tumor_group, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_tumor_group, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

pdf("SurvPlot_liver_meta_angiogenic_tumor.pdf", width=8, height=6)
surv_plot
dev.off()

# plot_obj <- surv_plot$plot
# ggsave(plot = plot_obj, "SurvPlot_liver_meta_angiogenic_tumor.pdf", width=6, height=6)

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ APCDD1_tumor_group, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ APCDD1_tumor_group, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

# plot_obj <- surv_plot$plot
# ggsave(plot = plot_obj, "SurvPlot_liver_meta_APCDD1_tumor.pdf", width=6, height=6)

pdf("SurvPlot_liver_meta_APCDD1_tumor.pdf", width=8, height=6)
surv_plot
dev.off()

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Proliferative_tumor_group, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_tumor_group, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

pdf("SurvPlot_liver_meta_Proliferative_tumor.pdf", width=8, height=6)
surv_plot
dev.off()

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Proliferative_stem_tumor_group, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_stem_tumor_group, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

pdf("SurvPlot_liver_meta_Proliferative_stem_tumor.pdf", width=8, height=6)
surv_plot
dev.off()

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Stem_like_tumor_group, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Stem_like_tumor_group, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

pdf("SurvPlot_liver_meta_Stem_like_tumor_tumor.pdf", width=8, height=6)
surv_plot
dev.off()

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Angiogenic_low_APCDD1, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_low_APCDD1, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        ylim = c(0.25, 1),
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,4)
surv_plot

pdf("SurvPlot_liver_meta_Angiogenic_low_APCDD1_tumor.pdf", width=7, height=4)
surv_plot
dev.off()

surv_obj <- Surv(time = liver_sample_tbl$Dx_to_death_or_last_fu_months, event = liver_sample_tbl$'Vital_status')

fit <- survfit(surv_obj ~ Angiogenic_APCDD1_cat, data = liver_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_APCDD1_cat, data = liver_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = liver_sample_tbl, 
                        palette = c("#7570B3", "#E7298A", "firebrick2", "navyblue"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(24, 36),
                        break.x.by = 3,    
                        ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,4)
surv_plot

pdf("SurvPlot_liver_meta_Angiogenic_combined_APCDD1_tumor.pdf", width=8, height=6)
surv_plot
dev.off()



surv_obj <- Surv(time = lung_sample_tbl$Dx_to_death_or_last_fu_months, event = lung_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Angiogenic_tumor_group, data = lung_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_tumor_group, data = lung_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = lung_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        #xlim = c(24, 36),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

surv_obj <- Surv(time = lung_sample_tbl$Dx_to_death_or_last_fu_months, event = lung_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ APCDD1_tumor_group, data = lung_sample_tbl)

surv_diff <- survdiff(surv_obj ~ APCDD1_tumor_group, data = lung_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = lung_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        #xlim = c(24, 36),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(33, 0.9)  
                       )

set_size(6,6)
surv_plot

# Prepare your survival object
surv_obj <- Surv(time = primary_sample_tbl$Dx_to_death_or_last_fu_months, event = primary_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Angiogenic_tumor_group, data = primary_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Angiogenic_tumor_group, data = primary_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = primary_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(72, 84),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(75, 0.9)  
                       )

set_size(6,6)
surv_plot

# Prepare your survival object
surv_obj <- Surv(time = primary_sample_tbl$Dx_to_death_or_last_fu_months, event = primary_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ APCDD1_tumor_group, data = primary_sample_tbl)

surv_diff <- survdiff(surv_obj ~ APCDD1_tumor_group, data = primary_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = primary_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(72, 84),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(75, 0.9)  
                       )

set_size(6,6)
surv_plot

# Prepare your survival object
surv_obj <- Surv(time = primary_sample_tbl$Dx_to_death_or_last_fu_months, event = primary_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Proliferative_tumor_group, data = primary_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_tumor_group, data = primary_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = primary_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(72, 84),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(75, 0.9)  
                       )

set_size(6,6)
surv_plot

# Prepare your survival object
surv_obj <- Surv(time = primary_sample_tbl$Dx_to_death_or_last_fu_months, event = primary_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Proliferative_stem_tumor_group, data = primary_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Proliferative_stem_tumor_group, data = primary_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = primary_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(72, 84),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(75, 0.9)  
                       )

set_size(6,6)
surv_plot

# Prepare your survival object
surv_obj <- Surv(time = primary_sample_tbl$Dx_to_death_or_last_fu_months, event = primary_sample_tbl$'Vital_status')

# Fit a survival curve
fit <- survfit(surv_obj ~ Stem_like_tumor_group, data = primary_sample_tbl)

surv_diff <- survdiff(surv_obj ~ Stem_like_tumor_group, data = primary_sample_tbl)
p_value <- round(1 - pchisq(surv_diff$chisq, df = 1), 2)

surv_plot <- ggsurvplot(fit, 
                        data = primary_sample_tbl, 
                        palette = c("#7570B3", "#E7298A"),
                        xlab = "Months",
                        ylab = "Survival Probability",
                        xlim = c(72, 84),
                        break.x.by = 3,    
                        #ggtheme = theme_minimal(),
                        censor = FALSE,
                        risk.table = TRUE,
                        pval = paste0("p value=", p_value),
                        pval.coord = c(75, 0.9)  
                       )

set_size(6,6)
surv_plot


