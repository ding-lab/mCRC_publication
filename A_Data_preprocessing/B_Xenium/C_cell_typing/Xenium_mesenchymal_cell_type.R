library(Seurat)
library(SeuratObject)
library(ggplot2)
library(circlize)
library(tidyverse)
library(patchwork)
library(qs)
library(reticulate)
library(harmony)
library(RColorBrewer)
library(BPCells)
library(future)
options(future.globals.maxSize = 150000 * 1024^2)
plan(multisession, workers = 30)

#Load R code and conda environment for secuer clustering

xenium_mesenchymal <- qread("/path/to/non_integrated_xenium_mesenchymal_subset.qs")
xenium_mesenchymal <- JoinLayers(xenium_mesenchymal)
xenium_mesenchymal <- NormalizeData(xenium_mesenchymal)
xenium_mesenchymal[["Xenium"]] <- split(xenium_mesenchymal[["Xenium"]], f = xenium_mesenchymal$orig.ident)
xenium_mesenchymal <- FindVariableFeatures(xenium_mesenchymal)
xenium_mesenchymal <- ScaleData(xenium_mesenchymal)
xenium_mesenchymal <- RunPCA(xenium_mesenchymal)
xenium_mesenchymal <- IntegrateLayers(
  object = xenium_mesenchymal, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)
xenium_mesenchymal <- RunUMAP(xenium_mesenchymal, reduction = "pca", dims = 1:30, reduction.name = "umap.pca")
xenium_mesenchymal <- RunUMAP(xenium_mesenchymal, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
fea = Embeddings(object = xenium_mesenchymal, reduction = "harmony")
dim(fea)
res <- sr$secuer(fea = fea,
                 Knn = as.integer(20),
                 multiProcessState = TRUE,
                 eskResolution = 4,
                 num_multiProcesses = as.integer(15))
mesenchymal_markers <- c("WNT5A", "WNT5B", "NRG1", "BMP4", "PDGFRA", "COL7A1",  "LOXL2","NDRG1", "TNC",  "POSTN", "FAP", "ACTA2", "CTHRC1", "THY1", "MYH11", "PDGFRB",  "RGS5", "PECAM1",
                         "CCL2", "C3", "CXCL12", "IL6ST")  

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
p2 = DotPlot(xenium_mesenchymal, features= mesenchymal_markers, group.by='secuer_res4', cluster.idents = TRUE) +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

xenium_mesenchymal@meta.data <- xenium_mesenchymal@meta.data %>% 
                        mutate(mesenchymal_cell_type = case_when(secuer_res4 %in% c( '14', '13', '5' ) ~ "smooth_muscle",
                                                     secuer_res4 == '3' ~ 'vascular_smooth_muscle',
                                                     secuer_res4 %in% c( '1', '17', '10' )  ~ "stromal_fibroblast",
                                                     secuer_res4 %in% c( '20', '9' ) ~ 'CD3_Fib_doublet',
                                                     secuer_res4 == '8' ~ 'CD68_Fib_doublet',          
                                                     secuer_res4 %in% c( '0', '2', '4' )~ 'iCAF',
                                                     secuer_res4 == '19' ~ 'WNT5A_BMP',
                                                     secuer_res4 %in% c( '18', '15', '22' ) ~ 'WNT5A_infl',
                                                     secuer_res4 %in% c( '11', '16', '6' ) ~ 'mCAF',
                                                     secuer_res4 %in% c( '21', '12', '7' ) ~ 'pericyte'))  

MetaData <- xenium_mesenchymal[[]]
MetaData$cell_id <- row.names(MetaData)
MetaData$cell_id <- row.names(MetaData)
MetaData <- MetadData[,c("cell_id", "mesenchymal_cell_type")]
write.csv(Metadata, "path/to/output/mesenchymal_cell_type.csv")

