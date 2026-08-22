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

out_dir = getwd()
paired_CRC_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tumor_paired_scaled_CoLi_TMR_association.qs'
paired_CRC <- qread(paired_CRC_path)
paired_CRC

# paired_CRC <- NormalizeData(paired_CRC)
# paired_CRC <- FindVariableFeatures(paired_CRC)
# paired_CRC <- ScaleData(paired_CRC)
# paired_CRC

xenium_anno = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
xenium_metadata = read.csv(xenium_anno, header = TRUE)
rownames(xenium_metadata) <- xenium_metadata$barcode
xenium_metadata <- xenium_metadata %>% select(-c(barcode))

paired_CRC@meta.data <- paired_CRC@meta.data %>% 
                        mutate(tmr_region = case_when(
                                   corrected_sample_id == 'CM819C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM819C1_Co' & tn_label %in% c(7, 8, 12, 14)  ~ 'muscle',
                                   corrected_sample_id == 'CM819C1_Co' & tn_label %in% c(47, 48, 49)  ~ 'mucosa',
                                   corrected_sample_id == 'CM819C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM579C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM579C1_Co' & tn_label %in% c(24) ~ 'muscle',
                                   corrected_sample_id == 'CM579C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM798C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM798C1_Co' & tn_label %in% c(1:9) ~ 'muscle',
                                   corrected_sample_id == 'CM798C1_Co' ~ 'submucosa',
                                   corrected_sample_id == 'CM397C1_Li' ~ 'metastasis',
                                   corrected_sample_id == 'CM397C1_Co' & tn_label %in% c(9, 10, 12, 13) ~ 'muscle',
                                   corrected_sample_id == 'CM397C1_Co' ~ 'submucosa'
                                   ))

CRC_samples = c('CM819C1_Co', 'CM579C1_Co', 'CM798C1_Co', 'CM397C1_Co')
mCRC_samples = c('CM819C1_Li', 'CM579C1_Li', 'CM798C1_Li', 'CM397C1_Li')

paired_CRC_tmr_count <- xenium_metadata %>%
                        filter(Tissue_ID %in% c('CM819C1', 'CM579C1', 'CM798C1', 'CM394C1')) %>%  
                        filter(tn_label !=0) %>% 
                        filter(is.na(tn_label) == FALSE) %>% 
                        mutate(tmr_anno = case_when(
                            Organ == 'Liver' ~ paste0(Tissue_ID, '_met_', tn_label),
                            Organ == 'Colon' & tmr_region == 'mucosa'  ~ paste0(Tissue_ID, '_muc_', tn_label),
                            Organ == 'Colon' & tmr_region == 'submucosa'  ~ paste0(Tissue_ID, '_sub_', tn_label),
                            Organ == 'Colon' & tmr_region == 'muscle'  ~ paste0(Tissue_ID, '_mus_', tn_label)))%>% 
                        count(tmr_anno)

paired_CRC_TMR_less_1000 = unique((paired_CRC_tmr_count %>% filter(n < 1000))$tmr_anno)

paired_CRC2 <- paired_CRC %>% subset(tn_label != 0) %>%
               subset(!(tmr_anno %in% paired_CRC_TMR_less_1000))  
paired_CRC2

DefaultAssay(paired_CRC) <- "Xenium"

var_genes <- VariableFeatures(paired_CRC)
top_genes <- head(var_genes, 5000)  

avg_list <- AverageExpression(
  paired_CRC,
  assays   = "Xenium",
  slot     = "data",           
  group.by = "sample_region"
)

mat <- avg_list$Xenium

mat <- as.matrix(mat)
mat <- mat[top_genes, , drop = FALSE]

var_genes <- VariableFeatures(paired_CRC)
top_genes <- head(var_genes, 1000)  

# Assume mat_top is genes (rows) x regions (columns)
mat_top <- mat[top_genes, , drop = FALSE]

# Extract tissue IDs
col_tissue_id <- sub("-.*$", "", colnames(mat_top))

# Create an empty matrix to store scaled values
scaled_mat <- matrix(NA, nrow = nrow(mat_top), ncol = ncol(mat_top))
rownames(scaled_mat) <- rownames(mat_top)
colnames(scaled_mat) <- colnames(mat_top)

# Loop through each tissue and scale the columns for each gene
for (tid in unique(col_tissue_id)) {
  cols <- which(col_tissue_id == tid)
  # Extract sub-matrix for this tissue
  sub_mat <- mat_top[, cols, drop = FALSE]
  # For each gene, scale across columns (regions)
  scaled_sub_mat <- t(apply(sub_mat, 1, scale))
  # Remove dimnames to avoid mismatch warning
  dimnames(scaled_sub_mat) <- NULL
  scaled_mat[, cols] <- scaled_sub_mat
}

cor_mat <- cor(scaled_mat, method = "pearson")

col_fun <- colorRamp2(c(-0.3, 0, 0.3), c("navy", "white", "firebrick"))

HpAll = Heatmap(
  cor_mat,
  name                   = "Pearson r",
  col                    = col_fun,
  clustering_distance_rows    = "euclidean",
  clustering_method_rows      = "ward.D2",
  clustering_distance_columns = "euclidean",
  clustering_method_columns   = "ward.D2"
)


pdf(file.path(out_dir, "Pooled_samples_tmr_annotation_correlation_heatmap_all_genes_sample_region_scale.data.pdf"), width = 12, height = 12)
draw(HpAll)
dev.off()