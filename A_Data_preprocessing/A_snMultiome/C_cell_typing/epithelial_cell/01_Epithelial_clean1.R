library(Seurat)
library(dplyr)
library(ggplot2)
library(MAST)
source("../../../scripts/jupyter_support_functions.R")

stem_cell_genes = c('LGR5', 'ASCL2', 'OLFM4', 'PPP1R1B', 'RNF43', 'SATB2', 'SELENBP1',
                    'SLC12A2', 'SLC26A2', 'SMOC2', 'TCEA3', 'TSPAN8', 'UGT2B17',
                    'ATP1A1', 'CCDC14', 'CD44', 'CDHR1', 'CFTR', 'EPHB2', 'HNF4A',
                    'KLF4')

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/tumor')
tumor = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_reINT.rds')
tumor

DefaultAssay(tumor) = 'SCT'
tumor <- FindClusters(tumor, resolution = 0.1, cluster.name = 'tumor_clusters_0.1', graph = 'RNA_snn')

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_0.1', label = TRUE)

table(tumor$tumor_clusters_0.1)

tumor <- tumor %>% subset(tumor_clusters_0.1 == '0' |
                          tumor_clusters_0.1 == '1' |
                          tumor_clusters_0.1 == '2' |
                          tumor_clusters_0.1 == '3' |
                          tumor_clusters_0.1 == '4' 
                         )

set_size(7,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'cell_type_simple3', label = TRUE)

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_0.1', label = TRUE)

set_size(6,6)
FeaturePlot(tumor, features = 'SLC2A1', reduction = 'tumor_umap.scvi')

set_size(6,6)
Idents(tumor)<- 'tumor_clusters_0.1'
DotPlot(tumor, features = c("MKI67", "CENPF", "TOP2A", "BRCA1", "LAMA3", "LAMC2", "LAMB3",
                            "FAM3D", "NDRG1", "VEGFA", "HSPH1",'RASSF6', 'BCAS1', 'SLC2A1', 'CHGA'))+ 
theme(axis.text.x = element_text(angle = 90),  panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))

Idents(tumor) = "tumor_clusters_0.1"
tumor = FindSubCluster(object = tumor,
                       cluster = '4',
                       graph.name = 'RNA_snn',
                       subcluster.name = "tumor_clusters_sub4",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(tumor$tumor_clusters_sub4)

table(tumor$tumor_clusters_sub4)

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_sub4', label = TRUE)

Idents(tumor) = "tumor_clusters_sub4"
tumor = FindSubCluster(object = tumor,
                       cluster = '2',
                       graph.name = 'RNA_snn',
                       subcluster.name = "tumor_clusters_sub4_2",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(tumor$tumor_clusters_sub4_2)

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_sub4_2', label = TRUE)

Idents(tumor) = "tumor_clusters_sub4_2"
tumor = FindSubCluster(object = tumor,
                       cluster = '0',
                       graph.name = 'RNA_snn',
                       subcluster.name = "tumor_clusters_sub4_2_0",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(tumor$tumor_clusters_sub4_2_0)

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_sub4_2_0', label = TRUE)

Idents(tumor) = "tumor_clusters_sub4_2_0"
tumor = FindSubCluster(object = tumor,
                       cluster = '1',
                       graph.name = 'RNA_snn',
                       subcluster.name = "tumor_clusters_0.1_sub",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(tumor$tumor_clusters_0.1_sub)

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_clusters_0.1_sub', label = TRUE)

table(tumor$tumor_clusters_0.1_sub)

tumor@meta.data <- tumor@meta.data %>% 
                   mutate(tumor_cluster = case_when(
                       tumor_clusters_0.1_sub == '0_0' ~ 0,
                       tumor_clusters_0.1_sub == '0_2' ~ 0,
                       tumor_clusters_0.1_sub == '0_3' ~ 0,
                       tumor_clusters_0.1_sub == '0_4' ~ 0,
                       tumor_clusters_0.1_sub == '0_1' ~ 1,
                       tumor_clusters_0.1_sub == '1_0' ~ 2,
                       tumor_clusters_0.1_sub == '1_2' ~ 2,
                       tumor_clusters_0.1_sub == '1_3' ~ 2,
                       tumor_clusters_0.1_sub == '1_1' ~ 3,
                       tumor_clusters_0.1_sub == '2_0' ~ 4,
                       tumor_clusters_0.1_sub == '2_1' ~ 4,
                       tumor_clusters_0.1_sub == '2_2' ~ 5,
                       tumor_clusters_0.1_sub == '3' ~ 6,
                       tumor_clusters_0.1_sub == '4_0' ~ 7,
                       tumor_clusters_0.1_sub == '4_2' ~ 7,
                       tumor_clusters_0.1_sub == '4_1' ~ 8
                   ))

set_size(6,6)
DimPlot(tumor, reduction = 'tumor_umap.scvi', group.by = 'tumor_cluster', label = TRUE)

set_size(6,6)
Idents(tumor)<- 'tumor_cluster'
DotPlot(tumor, features = c("MKI67", "CENPF", "TOP2A", "BRCA1", "LAMA3", "LAMC2", "LAMB3",
                            "FAM3D", "NDRG1", "VEGFA", "HSPH1",'RASSF6', 'BCAS1', 'SLC2A1', 'CHGA'))+ 
theme(axis.text.x = element_text(angle = 90),  panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"))

Idents(tumor) = "tumor_cluster"
tumor_markers = FindAllMarkers(tumor, 
                               min.pct=0.1, 
                               logfc.threshold=0.5, 
                               only.pos = TRUE, 
                               return.thresh = 0.01,
                               test.use = 'MAST'
                              )

write.csv(tumor_markers, "mCRC_tumor_mast_deg.csv", row.names = TRUE)

tumor_markers %>% filter(gene == 'MKI67')


# cluster 7 with T cell genes, favored doublets with T cells
# cluster 8 with abundant of mito genes, favored dead cells

tumor2 <- tumor %>% subset(tumor_cluster != '7' &
                           tumor_cluster != '8'
                          ) 

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/tumor')
saveRDS(tumor2, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean1.rds') 

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
epithelial = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_reINT.rds')
epithelial

DefaultAssay(epithelial) = 'SCT'
epithelial <- FindClusters(epithelial, resolution = 0.1, cluster.name = 'epithelial_clusters_0.1', graph = 'RNA_snn')

set_size(6,6)
DimPlot(epithelial, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)

tumor.metadata = tumor@meta.data

epithelial<- AddMetaData(object = epithelial, metadata = tumor.metadata, col.name = 'tumor_cluster')

colnames(epithelial@meta.data)

epithelial@meta.data <- epithelial@meta.data %>% 
  mutate(epithelial_cluster = if_else(
    !is.na(tumor_cluster), 
    as.character(tumor_cluster), 
    cell_type_simple3
  ))

set_size(6,6)
DimPlot(epithelial, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster', label = TRUE)

epithelial2 = epithelial %>% subset(epithelial_cluster != '7' &
                                    epithelial_cluster != '8' &
                                    epithelial_cluster != 'Tumor'
                                   )

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster', label = TRUE)

Idents(epithelial2) = "epithelial_cluster"
epithelial2 = FindSubCluster(object = epithelial2,
                       cluster = 'Enterocyte',
                       graph.name = 'RNA_snn',
                       subcluster.name = "epithelial_cluster_esub",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(epithelial2$epithelial_cluster_gsub)

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_esub', label = TRUE)

table(epithelial2$epithelial_cluster_esub)

Idents(epithelial2) = "epithelial_cluster_esub"
epithelial2 = FindSubCluster(object = epithelial2,
                       cluster = 'Goblet',
                       graph.name = 'RNA_snn',
                       subcluster.name = "epithelial_cluster_egsub",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(epithelial2$epithelial_cluster_egsub)

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub', label = TRUE)

table(epithelial2$epithelial_cluster_egsub)

epithelial2 <- epithelial2 %>% subset(epithelial_cluster_egsub == '0' |
                                      epithelial_cluster_egsub == '1' |
                                      epithelial_cluster_egsub == '2' |
                                      epithelial_cluster_egsub == '3' |
                                      epithelial_cluster_egsub == '4' |
                                      epithelial_cluster_egsub == '5' |
                                      epithelial_cluster_egsub == '6' |
                                      epithelial_cluster_egsub == 'Enterocyte_0' |
                                      epithelial_cluster_egsub == 'Enterocyte_1' |
                                      epithelial_cluster_egsub == 'Goblet_0' |
                                      epithelial_cluster_egsub == 'Goblet_1' |
                                      epithelial_cluster_egsub == 'Goblet_2'
                                     )

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub', label = TRUE)

Idents(epithelial2) = "epithelial_cluster_egsub"
epithelial2 = FindSubCluster(object = epithelial2,
                       cluster = 'Enterocyte_1',
                       graph.name = 'RNA_snn',
                       subcluster.name = "epithelial_cluster_egsub2",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(epithelial2$epithelial_cluster_egsub2)

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub2', label = TRUE)

Idents(epithelial2) = "epithelial_cluster_egsub2"
epithelial2 = FindSubCluster(object = epithelial2,
                       cluster = 'Goblet_2',
                       graph.name = 'RNA_snn',
                       subcluster.name = "epithelial_cluster_egsub3",
                       resolution = 0.1,
                       algorithm = 1
                      )
unique(epithelial2$epithelial_cluster_egsub3)

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub3', label = TRUE)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
saveRDS(epithelial2, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1.rds') 

Idents(epithelial2) = "epithelial_cluster_egsub3"
epithelial2_markers = FindAllMarkers(epithelial2, 
                                     min.pct=0.1, 
                                     logfc.threshold=0.5, 
                                     only.pos = TRUE, 
                                     return.thresh = 0.01,
                                     test.use = 'MAST'
                                    )

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
write.csv(epithelial2_markers, "mCRC_epithelial_mast_deg.csv", row.names = TRUE)

epithelial_cluster = epithelial2@meta.data %>% select('epithelial_cluster_egsub3')
write.csv(epithelial_cluster,  "mCRC_epithelial_cluster.metadata.csv", row.names = TRUE)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
epithelial2 = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1.rds') 

set_size(6,6)
DimPlot(epithelial2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub3', label = TRUE)

unique(epithelial2$orig.ident)

set_size(6,6)
Highlight_Cluster_UMAP(epithelial2, 
                       metadata_column='orig.ident', 
                       cell_type='HT413C1-Th1K4A2Nd1_2Bma1_1', 
                       umap_reduction = 'epithelial_umap.scvi')

set_size(6,6)
Highlight_Cluster_UMAP(epithelial2, 
                       metadata_column='orig.ident', 
                       cell_type='HT413C1-Th1K2A2Nd1_2Bma1_1', 
                       umap_reduction = 'epithelial_umap.scvi')

set_size(6,6)
FeaturePlot(tumor, features = 'SLC2A1', reduction = 'tumor_umap.scvi')

#BEST4+ enterocytes (BEST4, OTOP2), 
#goblet cells (MUC2, TFF1, SYTL2), 
#immature goblet cells (KLK1, RETNLB, CLCA1), 
#stem cells (RGMB, SMOC2, LGR5, ASCL2), 
#tuft cells (SH2D6, TRPM5, BMX, LRMP, HCK), 
#enteroendocrine cells (SCGN, FEV, CHGA, PYY, GCG), 
#cycling transit-amplifying cells (TICRR, CDC25C),
#Paneth cells (LYZ, DEFA5)
# Nature volume 619, pages572–584 (2023)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
epithelial = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1_reINT.rds')
epithelial

set_size(6,6)
DimPlot(epithelial, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub3', label = TRUE)

DefaultAssay(epithelial) = 'SCT'
epithelial <- FindClusters(epithelial, resolution = 0.1, cluster.name = 'epithelial_clusters_0.1', graph = 'RNA_snn')

set_size(6,6)
DimPlot(epithelial, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)

set_size(8,8)
FeaturePlot(epithelial, features = c('RGMB', 'SMOC2', 'LGR5', 'ASCL2'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
DefaultAssay(epithelial) <- 'RNA'
FeaturePlot(epithelial, features = c('SLC2A1', 'VEGFA'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
DefaultAssay(epithelial) <- 'SCT'
FeaturePlot(epithelial, features = c('SLC2A1', 'VEGFA'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial, features = c('TICRR', 'CDC25C', 'MKI67', 'KRAS'), reduction = 'epithelial_umap.scvi')


