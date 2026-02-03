library(Seurat)
library(dplyr)
library(ggplot2)
library(MAST)
library(readr)
source("../../../scripts/jupyter_support_functions.R")

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
epithelial_clean = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1.rds')
epithelial_clean

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial')
epithelial_reint = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1_reINT.rds')
epithelial_reint

mutation_mapping_path <- "/diskmnt/Projects/MetNet_analysis_2/Colorectal/epeng/10Xmapping/metadata_mapped_mutations.tsv"
mutation_mapping <- as.data.frame(read_tsv(mutation_mapping_path, col_names = TRUE))
rownames(mutation_mapping) <- mutation_mapping[[1]]
mutation_mapping <- mutation_mapping %>% select(APC_ref, APC_var, APC_mut, KRAS_ref, KRAS_var, KRAS_mut, TP53_mut, TP53_var, TP53_mut)
head(mutation_mapping)

type(mutation_mapping$APC_mut)

wxs_epithelial <- subset(epithelial_reint, 
                         subset = (orig.ident == 'CM268C1-S1' | orig.ident == 'CM268C1-T1' | orig.ident == 'CM354C1-T1' |
                                   orig.ident == 'CM354C2-T1' | orig.ident == 'CM268C1-T1' | orig.ident == 'CM392C1-S1' |
                                   orig.ident == 'CM392C2-Th1' | orig.ident == 'CM426C1-Th1' | orig.ident == 'CM426C2-Tp1' |
                                   orig.ident == 'CM478C1-T1Y2' | orig.ident == 'HT213C1-Te1' | orig.ident == 'HT225C1-Th1' |
                                   orig.ident == 'HT230C1-Th1' | orig.ident == 'HT253C1-Th1' | orig.ident == 'HT260C1-Th1' |
                                   orig.ident == 'HT266C1-Tb1'))

mutation_mapping_toadd <- mutation_mapping %>% select(APC_mut, KRAS_mut, TP53_mut)
wxs_epithelial <- AddMetaData(wxs_epithelial, metadata = mutation_mapping_toadd)                                  

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

unique(wxs_epithelial$TP53_mut)

p <- DimPlot(wxs_epithelial, group.by = "TP53_mut", reduction = "epithelial_umap.scvi", order = T, label=F, pt.size=1, label.size=6, raster=FALSE) +
    labs(title = "TP53_alleles", color="Allele_Detected") + 
    scale_color_manual(labels = c("NA", "TP53 variant"), values = c("grey", "red"))
pdf('TP53_alleles_umap.pdf', height = 6, width = 6)
p
dev.off()

set_size(6,6)
p

set_size(6,6)
p1 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_egsub3', label = TRUE)
p1

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('OTOP2', 'BEST4'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('MUC2', 'TFF1', 'SYTL2'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('KLK1', 'RETNLB', 'CLCA1'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('RGMB', 'SMOC2', 'LGR5', 'ASCL2'), reduction = 'epithelial_umap.scvi')

set_size(8,12)
FeaturePlot(epithelial_reint, features = c('SH2D6', 'TRPM5', 'BMX', 'LRMP', 'HCK'), reduction = 'epithelial_umap.scvi')

set_size(8,12)
FeaturePlot(epithelial_reint, features = c('SCGN', 'FEV', 'CHGA', 'CHGB', 'PYY', 'GCG'), reduction = 'epithelial_umap.scvi')

set_size(8,8)
FeaturePlot(epithelial_reint, features = c('TICRR', 'CDC25C', 'MKI67', 'TOP2A'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('LYZ', 'DEFA5'), reduction = 'epithelial_umap.scvi')

set_size(8,4)
FeaturePlot(epithelial_reint, features = c('SLC2A1', 'VEGFA'), reduction = 'epithelial_umap.scvi')

epithelial_reint <- FindClusters(epithelial_reint, 
                                 resolution = 0.1, 
                                 cluster.name = 'epithelial_clusters_0.1', 
                                 graph = 'RNA_snn')

set_size(12,12)
p2 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)
p3 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'Site_of_Origin', label = TRUE)
p4 = DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'Tissue_Type', label = TRUE)
p1 + p2  + p3 + p4

set_size(10,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'orig.ident', label = TRUE)

set_size(6,6)
Highlight_Cluster_UMAP(seurat_object = epithelial_reint,
                       metadata_column = 'orig.ident',
                       cell_type = 'CM426C1-Th1',
                       umap_reduction = 'epithelial_umap.scvi',
                       highlight_color = 'forestgreen'
                      )

Idents(epithelial_reint) = "epithelial_cluster_egsub3"
cluster6_markers = FindMarkers(epithelial_reint,
                               ident.1 = '6',
                               ident.2 = '2',
                               min.pct=0.1, 
                               logfc.threshold=0.5, 
                               only.pos = TRUE, 
                               return.thresh = 0.01,
                               test.use = 'MAST')

write.csv(cluster6_markers,  "mCRC_epithelial_cluster6_2_mast_deg.csv", row.names = TRUE)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)

set_size(4,4)
DefaultAssay(epithelial_clean2) <- 'RNA' 
DotPlot(epithelial_clean2, feature = c('ITGAM', 'LYZ', 'PTPRC', 'SPP1'), group.by='epithelial_clusters_0.1')

# cluster 5: a doublet cluster of microglia

table(epithelial_reint$epithelial_clusters_0.1)

epithelial_clean2 <- epithelial_reint %>% subset(epithelial_clusters_0.1 == '0' |
                                                 epithelial_clusters_0.1 == '1' |
                                                 epithelial_clusters_0.1 == '2' |
                                                 epithelial_clusters_0.1 == '3' |
                                                 epithelial_clusters_0.1 == '4' |
                                                 epithelial_clusters_0.1 == '9' 
                                                 )

epithelial_clean2 <- FindClusters(epithelial_clean2, 
                                  resolution = 0.1, 
                                  cluster.name = 'epithelial_clusters_0.1', 
                                  graph = 'RNA_snn')

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)

table(epithelial_clean2$epithelial_clusters_0.1)

epithelial_clean2 <- epithelial_clean2 %>% subset(epithelial_clusters_0.1 == '0' |
                                                  epithelial_clusters_0.1 == '1' |
                                                  epithelial_clusters_0.1 == '2' |
                                                  epithelial_clusters_0.1 == '3' |
                                                  epithelial_clusters_0.1 == '4' |
                                                  epithelial_clusters_0.1 == '5' )

DefaultAssay(epithelial_clean2) <- 'SCT'

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_clusters_0.1', label = TRUE)

Idents(epithelial_clean2) = "epithelial_clusters_0.1"
epithelial_clean2 = FindSubCluster(object = epithelial_clean2,
                                   cluster = '3',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub1",
                                   resolution = 0.1,
                                   algorithm = 1
                                  )
unique(epithelial_clean2$epithelial_cluster_sub1)

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub1', label = TRUE)

Idents(epithelial_clean2) = "epithelial_cluster_sub1"
epithelial_clean2 = FindSubCluster(object = epithelial_clean2,
                                   cluster = '3_1',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub2",
                                   resolution = 0.3,
                                   algorithm = 1
                                  )
unique(epithelial_clean2$epithelial_cluster_sub2)

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub2', label = TRUE)

epithelial_clean2@meta.data <- epithelial_clean2@meta.data %>% 
                               mutate(epithelial_cluster_sub2 = case_when(epithelial_cluster_sub2 == '3_1_0' ~ '3_1',
                                                                          epithelial_cluster_sub2 == '3_1_2' ~ '3_1',
                                                                          epithelial_cluster_sub2 == '3_1_4' ~ '3_1',
                                                                          TRUE ~ epithelial_cluster_sub2))

DotPlot(epithelial_clean2, feature = 'LGR5', group.by = 'epithelial_cluster_sub2')

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub2', label = TRUE)

Idents(epithelial_clean2) = "epithelial_cluster_sub2"
epithelial_clean2 = FindSubCluster(object = epithelial_clean2,
                                   cluster = '0',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub3",
                                   resolution = 0.1,
                                   algorithm = 1
                                  )
unique(epithelial_clean2$epithelial_cluster_sub3)

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub3', label = TRUE)

DotPlot(epithelial_clean2, feature = 'LGR5', group.by = 'epithelial_cluster_sub3')

Idents(epithelial_clean2) = "epithelial_cluster_sub3"
epithelial_clean2 = FindSubCluster(object = epithelial_clean2,
                                   cluster = '1',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub4",
                                   resolution = 0.1,
                                   algorithm = 1
                                  )
unique(epithelial_clean2$epithelial_cluster_sub4)

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub4', label = TRUE)

table(epithelial_clean2$epithelial_cluster_sub4)

epithelial_clean2 <- epithelial_clean2 %>% subset(epithelial_cluster_sub4 != '1_3')

epithelial_clean2@meta.data <- epithelial_clean2@meta.data %>% 
                               mutate(epithelial_cluster_sub4 = case_when(epithelial_cluster_sub4 == '1_2' ~ '1_0',
                                                                          TRUE ~ epithelial_cluster_sub4))

set_size(6,6)
DimPlot(epithelial_clean2, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub4', label = TRUE)

Idents(epithelial_clean2) = "epithelial_cluster_sub4"
cluster_1_1_markers = FindMarkers(epithelial_clean2,
                                  ident.1 = '1_1',
                                  ident.2 = '1_0',
                                  min.pct=0.1, 
                                  logfc.threshold=0.5, 
                                  only.pos = TRUE, 
                                  return.thresh = 0.01,
                                  test.use = 'MAST')
write.csv(cluster_1_1_markers,  "mCRC_epithelial_cluster_1_1_mast_deg.csv", row.names = TRUE)

epithelial_clean3 <- epithelial_clean2 %>% 
                     subset(epithelial_cluster_sub4 != '1_1') %>% 
                     RunUMAP(reduction = "integrated.scvi", 
                             dims = 1:30, 
                             reduction.name = "epithelial_umap2.scvi")

set_size(6,6)
DimPlot(epithelial_clean3, reduction = 'epithelial_umap2.scvi', group.by = 'epithelial_cluster_sub4', label = TRUE)

Idents(epithelial_clean3) = "epithelial_cluster_sub4"
epithelial_clean3 = FindSubCluster(object = epithelial_clean3,
                                   cluster = '2',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub5",
                                   resolution = 0.2,
                                   algorithm = 1
                                  )
unique(epithelial_clean3$epithelial_cluster_sub5)

set_size(6,6)
DimPlot(epithelial_clean3, reduction = 'epithelial_umap2.scvi', group.by = 'epithelial_cluster_sub5', label = TRUE)

table(epithelial_clean3$epithelial_cluster_sub5)

epithelial_clean3 <- epithelial_clean3 %>% subset(epithelial_cluster_sub5 != '2_4')

Idents(epithelial_clean3) = "epithelial_cluster_sub5"
epithelial_clean3 = FindSubCluster(object = epithelial_clean3,
                                   cluster = '2_3',
                                   graph.name = 'RNA_snn',
                                   subcluster.name = "epithelial_cluster_sub6",
                                   resolution = 0.1,
                                   algorithm = 1
                                  )
unique(epithelial_clean3$epithelial_cluster_sub6)

set_size(6,6)
DimPlot(epithelial_clean3, reduction = 'epithelial_umap2.scvi', group.by = 'epithelial_cluster_sub6', label = TRUE)

DotPlot(epithelial_clean3, feature = c('TRPM5', 'LRMP'), group.by = 'epithelial_cluster_sub6')

epithelial_clean3@meta.data <- epithelial_clean3@meta.data %>% 
                               mutate(epithelial_cluster_sub6 = case_when(epithelial_cluster_sub6 == '2_2' ~ '2_0',
                                                                          epithelial_cluster_sub6 == '2_1' ~ '2_0',
                                                                          epithelial_cluster_sub6 == '2_3_2' ~ '2_3',
                                                                          epithelial_cluster_sub6 == '2_3_0' ~ '2_3',
                                                                          TRUE ~ epithelial_cluster_sub6))

set_size(6,6)
DimPlot(epithelial_clean3, reduction = 'epithelial_umap2.scvi', group.by = 'epithelial_cluster_sub6', label = TRUE)

colnames(epithelial_clean3@meta.data)

# epithelial_clean3$mCRCv5_clusters_2 <- NULL
# epithelial_clean3$tumor_cluster <- NULL
# epithelial_clean3$epithelial_cluster <- NULL
# epithelial_clean3$epithelial_cluster_gsub <- NULL
# epithelial_clean3$epithelial_cluster_esub <- NULL
# epithelial_clean3$epithelial_cluster_egsub <- NULL
# epithelial_clean3$epithelial_cluster_egsub2 <- NULL
# epithelial_clean3$epithelial_cluster_egsub3 <- NULL
# epithelial_clean3$epithelial_clusters_0.2 <- NULL
# epithelial_clean3$epithelial_clusters_0.4 <- NULL
# epithelial_clean3$epithelial_clusters_0.5 <- NULL
# epithelial_clean3$epithelial_clusters_0.6 <- NULL
# epithelial_clean3$epithelial_clusters_0.8 <- NULL
# epithelial_clean3$epithelial_clusters_1 <- NULL
# epithelial_clean3$epithelial_clusters_1.2 <- NULL
# epithelial_clean3$epithelial_clusters_1.4 <- NULL
# epithelial_clean3$epithelial_clusters_1.6 <- NULL
# epithelial_clean3$epithelial_clusters_1.8 <- NULL
# epithelial_clean3$epithelial_clusters_2 <- NULL
# epithelial_clean3$mCRCv6_clusters_2_sub65_16 <- NULL
# epithelial_clean3$mCRCv6_clusters_2_sub65_16_34 <- NULL
# epithelial_clean3$tumor_in_normal <- NULL
# epithelial_clean3$epithelial_cluster_sub1 <- NULL
# epithelial_clean3$epithelial_cluster_sub2 <- NULL
# epithelial_clean3$epithelial_cluster_sub3 <- NULL
# epithelial_clean3$epithelial_cluster_sub4 <- NULL
# epithelial_clean3$epithelial_cluster_sub5 <- NULL

saveRDS(epithelial_clean3, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean3.rds')

Idents(epithelial_clean3) = "epithelial_cluster_sub6"
epithelial_clean3_markers = FindAllMarkers(epithelial_clean3,
                                           min.pct=0.1, 
                                           logfc.threshold=0.5, 
                                           only.pos = TRUE, 
                                           return.thresh = 0.01,
                                           test.use = 'MAST')
write.csv(epithelial_clean3_markers,  "mCRC_epithelial_clean3_markers_mast_deg.csv", row.names = TRUE)


