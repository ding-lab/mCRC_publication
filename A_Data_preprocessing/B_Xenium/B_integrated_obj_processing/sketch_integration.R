library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)
library(harmony)
options(future.globals.maxSize = 20 * 1024^3)

rds_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_all_cells.qs'
obj = qread(rds_path)
obj@images <- list()

xenium_annotation_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/'
xenium_annotation_path = file.path(xenium_annotation_dir, 'mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv')
xenium_annotation = read.csv(xenium_annotation_path, header = TRUE)
rownames(xenium_annotation) <- xenium_annotation$barcode
obj <- AddMetaData(obj, xenium_annotation)

obj <- obj %>% subset(All_cell_type1 %in% c('NK', 'smooth_muscle', 'CD4_Treg', 'CD4_Th17_T', 'CD8_Effector_T', 
                                            'CD4_naive/CM_T', 'B_cell', 'stromal_fibroblast', 'CD4_Resting_T', 
                                            'Plasma_cell', 'LAMP3_DC', 'Proliferative_T', 'classical_monocyte', 
                                            'cDC2', 'cDC1', 'pDC', 'CD4_Th1-like', 'Intestinal_epithelium', 'Enteric_neuron', 
                                            'RTM_KC', 'Activated_EC', 'M2_TAM', 'Alveolar_EC', 'TAM', 'CD8_TRM_T', 'iCAF', 
                                            'mCAF', 'Proliferative-like', 'Stem-like', 'Non-canonical', 'Intestine-like', 
                                            'PMN-like_monocyte', 'RTM_alveolar', 'pericyte', 'WNT5A_infl', 
                                            'IFN_TAM', 'SPP1_TAM', 'EC', 'Angiogenic_EC', 'WNT5A_BMP', 
                                            'Tip_EC', 'Arterial_EC', 'Lymphatic_EC', 'Hepatocyte', 'LSEC', 'Aveolar_epithelium', 
                                            'Bronchial_gland', 'Cholangiocyte', 'APCDD1', 'Breast_duct'))

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures=2000)

obj <- SketchData(object = obj, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch", features = VariableFeatures(obj))

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj, verbose = F)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, verbose = F)

obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = TRUE
)

obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, return.model = T, verbose = F)

p1 <- DimPlot(obj, group.by = "All_cell_type1", reduction = "umap")

obj <- ProjectIntegration(object = obj, sketched.assay = "sketch", assay = "Xenium", reduction = "harmony")
obj <- ProjectData(object = obj, 
                   sketched.assay = "sketch", 
                   assay = "Xenium", 
                   sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", 
                   dims = 1:30, 
                   refdata = list(celltype.full = "All_cell_type1"))

qsave(obj, '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_all_cells_noDB_normalized_noFOVs_sketch.qs')