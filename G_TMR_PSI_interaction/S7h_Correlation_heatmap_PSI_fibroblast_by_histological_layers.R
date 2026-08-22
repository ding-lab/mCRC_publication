library(qs)
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(Matrix)
library(circlize)
library(viridis)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## Data preprocessing
# all_cell_obj_path = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_all_cells.qs"
# all_cell_obj <- qread(all_cell_obj_path)

# xenium_anno = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
# xenium_metadata = read.csv(xenium_anno, header = TRUE)
# rownames(xenium_metadata) <- xenium_metadata$barcode
# xenium_metadata <- xenium_metadata %>% select(-c(barcode))

# all_cell_obj <- AddMetaData(all_cell_obj, xenium_metadata)

# target_cells_metadata = all_cell_obj@meta.data %>% 
#                         dplyr::filter(corrected_sample_id %in% c('CM819C1_Co', 'CM819C1_Li', 'CM579C1_Co', 'CM579C1_Li',
#                                                                  'CM397C1_Co', 'CM397C1_Li', 'CM798C1_Co', 'CM798C1_Li')) %>% 
#                         dplyr::filter((tn_label != 0 & !is.na(tn_label) & neighborhoods == "NB4_tumor_stroma_interface") | 
#                                       (tn_label == 0 & tn_outward %in% c(1,2,3,4,5))) %>% 
#                         dplyr::filter(Broad_cell_type1== 'Fibroblast')

# target_cells = row.names(target_cells_metadata)
# obj <- subset(all_cell_obj, cells = target_cells)

# obj <- NormalizeData(obj)
# obj <- JoinLayers(obj)
# obj <- FindVariableFeatures(obj)
# obj <- ScaleData(obj)

# obj_path = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tmr_association_fibroblast_obj.qs"
# qsave(obj, obj_path)

obj_path = "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tmr_association_fibroblast_obj.qs"
obj <- qread(obj_path)

obj@meta.data <- obj@meta.data %>% 
                 mutate(near_tmr_region = case_when(
                     Organ == 'Liver' ~ 'metastasis',
                     Organ == 'Colon' & (near_mucosa_TMR == TRUE | tmr_region == 'mucosa') ~ 'mucosa',
                     Organ == 'Colon' & (near_submucosa_TMR == TRUE | tmr_region == 'submucosa')  ~ 'submucosa',
                     Organ == 'Colon' & (near_muscle_TMR == TRUE | tmr_region == 'muscle')  ~ 'muscle'))

obj@meta.data <- obj@meta.data %>% 
                        mutate(sample_region = paste0(Tissue_ID, '_', near_tmr_region))


DefaultAssay(obj) <- "Xenium"

var_genes <- VariableFeatures(obj)
top_genes <- head(var_genes, 1000) 

avg_list <- AverageExpression(
  obj,
  assays   = "Xenium",
  slot     = "data",           
  group.by = "sample_region"
)

mat <- avg_list$Xenium

mat <- as.matrix(mat)
mat_top <- mat[top_genes, , drop = FALSE]
cor_mat <- cor(mat_top, method = "pearson")

col_fun <- colorRamp2(c(0.5, 0.75, 1), c("navy", "white", "firebrick"))

HpAll = Heatmap(
  cor_mat,
  name                   = "Pearson r",
  col                    = col_fun,
  clustering_distance_rows    = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_distance_columns = "euclidean",
  clustering_method_columns   = "ward.D2"
)

out_dir = getwd()
pdf(file.path(out_dir, "Pooled_samples_tmr_fibroblast_correlation_heatmap_top1000_genes_sample_region.pdf"), width = 12, height = 12)
draw(HpAll)
dev.off()





