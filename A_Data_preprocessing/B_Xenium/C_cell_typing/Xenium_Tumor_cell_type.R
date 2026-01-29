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

# Sketch Integration
xenium_tumor <- qread("/path/to/non_integrated_xenium_tumor_subset.qs")
xenium_tumor <- JoinLayers(xenium_tumor)
xenium_tumor <- NormalizeData(xenium_tumor)
xenium_tumor[["Xenium"]] <- split(xenium_tumor[["Xenium"]], f = xenium_tumor$orig.ident)
xenium_tumor <- FindVariableFeatures(xenium_tumor, selection.method = "vst", nfeatures = 1000)
xenium_tumor <- SketchData(object = xenium_tumor, ncells = 5000, method = "LeverageScore", features = VariableFeatures(xenium_tumor), sketched.assay = "sketch", seed = 123L,
  cast = "dgCMatrix", verbose = TRUE)
DefaultAssay(xenium_tumor) <- "sketch"
xenium_tumor <- FindVariableFeatures(xenium_tumor)
xenium_tumor <- ScaleData(xenium_tumor)
xenium_tumor  <- RunPCA(xenium_tumor)
xenium_tumor <- IntegrateLayers(
  object = xenium_tumor, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony"
)

# Tumor cell typing in sketch integrated object
xenium_tumor <- FindNeighbors(xenium_tumor, reduction = "harmony", dims = 1:30)
xenium_tumor <- FindClusters(xenium_tumor, resolution = 0.2)
DefaultAssay(xenium_tumor) <- "sketch"
Idents(xenium_tumor) <- "sketch_snn_res.0.2"
xenium_tumor <- FindSubCluster(
  xenium_tumor,
  cluster=4,
  graph.name="sketch_snn",
  subcluster.name = "subcluster4_0.1",
  resolution = 0.1,
  algorithm = 1
)
xenium_tumor <- RunUMAP(xenium_tumor, reduction = "harmony", dims = 1:30, return.model = T)
tumor_marker_genes <- c("LGR5", "SMOC2", "ASCL2", "RGMB","EPCAM","EPHB2","SOX9", "CDX2", "CDX1", "SLC26A3", "NOX1", "MUC2" , "MKI67", "TOP2A", "MYC", "BIRC5", "APCDD1", "VEGFA" , "SLC2A1", "EMP1","TGFBI", "ANXA1", "CLDN4", "YAP1", "L1CAM", "ACTA2", "POSTN", "CD19", "CD3E", "CD68")

DefaultAssay(xenium_tumor) <- "sketch"
Idents(xenium_tumor) = "sketch_snn_res.0.3"
rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
p2 = DotPlot(xenium_tumor, features= tumor_marker_genes , group.by='sketch_snn_res.0.3') +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

xenium_tumor@meta.data <- xenium_tumor@meta.data %>% 
                        mutate(tumor_cell_sketch2 = case_when(subcluster4_0.1 == '0' ~ 'Proliferative_like',
                                                     subcluster4_0.1 == '1' ~ 'Stem_like',
                                                     subcluster4_0.1 == '2' ~ 'Non_canonical',
                                                     subcluster4_0.1 == '3' ~ 'Intestine_like',
                                                     subcluster4_0.1 == '4_0' ~ 'Mesenchymal_doublet',
                                                     subcluster4_0.1 == '4_1' ~ 'Non_canonical',
                                                     subcluster4_0.1 == '4_2' ~ 'Proliferative_like',
                                                     subcluster4_0.1 == '5' ~ 'APCDD1',
                                                     subcluster4_0.1 == '6' ~ 'Non_canonical',
                                                     subcluster4_0.1 == '7' ~ 'Non_canonical'))
xenium_tumor$tumor_cell_sketch2 <- factor(xenium_tumor$tumor_cell_sketch2, levels = rev(c("Stem_like", "Intestine_like", "Proliferative_like", "Non_canonical", "APCDD1", "Mesenchymal_doublet")))

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
p1 = DotPlot(xenium_tumor, features= tumor_marker_genes , group.by='tumor_cell_sketch') +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

# Project tumor cell identities from sketch object to full object                                                    
DefaultAssay(xenium_tumor) <- "Xenium"
xenium_tumor <- ProjectIntegration(object = xenium_tumor, sketched.assay = "sketch", assay = "Xenium", reduction = "harmony")
xenium_tumor <- ProjectData(object = xenium_tumor, sketched.assay = "sketch", assay = "Xenium", sketched.reduction = "harmony.full",
    full.reduction = "harmony.full", dims = 1:30, refdata = list(celltype_subcluster = "tumor_cell_sketch2"))
xenium_tumor <- RunUMAP(xenium_tumor, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
    reduction.key = "UMAP_full_")

MetaData <- xenium_tumor[[]]
MetaData$tumor_cell_type <- MetaData$celltype_subcluster
MetaData$cell_id <- row.names(MetaData)
MetaData <- MetadData[,c("cell_id", "tumor_cell_type")]
write.csv(Metadata, "path/to/output/tumor_cell_type.csv"
