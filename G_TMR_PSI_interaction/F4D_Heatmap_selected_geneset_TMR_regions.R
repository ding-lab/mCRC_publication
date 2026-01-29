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

tmr_col <- c(submucosa = 'palegreen4', muscle = 'firebrick1', metastasis = 'violetred1')

out_dir = getwd()
paired_CRC_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_tumor_paired_scaled_CoLi_TMR_association.qs'
paired_CRC <- qread(paired_CRC_path)
paired_CRC

xenium_anno = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
xenium_metadata = read.csv(xenium_anno, header = TRUE)
rownames(xenium_metadata) <- xenium_metadata$barcode
xenium_metadata <- xenium_metadata %>% select(-c(barcode))

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

tumor_auc_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_Xenium5K_N26_tumor_genesets_auc_score_SampleID_corrected.csv'
tumor_auc = read.csv(tumor_auc_file, header = TRUE)
#tumor_auc = tumor_auc %>% separate(X, into = c("SampleID", "cell_id"), sep = "__")

genesets = c('HRC', 'revCSC',
             'Hallmark_angiogenesis', 'Hallmark_EMT', 'Hallmark_TGFB',
             'Hallmark_PI3K', 'Hallmark_KRAS_up')

tumor_auc_subcols = tumor_auc[, c('SampleID', 'SampleID2', 'cell_id', genesets)]
tumor_auc_subcols = tumor_auc_subcols %>% 
						filter(SampleID2 %in% c('CM819C1_Co', 'CM819C1_Li', 'CM397C1_Co', 'CM397C1_Li', 'CM798C1_Co', 'CM798C1_Li', 'CM579C1_Co', 'CM579C1_Li'))
tumor_auc_subcols$barcode <- paste0(tumor_auc_subcols$SampleID, '__', tumor_auc_subcols$cell_id)
tumor_auc_subcols <- tumor_auc_subcols %>% select(-c(cell_id, SampleID, SampleID2))
rownames(tumor_auc_subcols) <- tumor_auc_subcols$barcode

paired_CRC_data = paired_CRC@meta.data
paired_CRC_data$barcode = rownames(paired_CRC_data)

paired_CRC_data = left_join(paired_CRC_data, tumor_auc_subcols, by = 'barcode') %>% 
                  mutate(sample_tmr_region = paste0(Tissue_ID, '_', tmr_region))

paired_CRC_median <- paired_CRC_data %>%
  filter(All_cell_type2 %in% c('Intestine-like', 'Proliferative-like', 'Stem-like', 'Non-canonical')) %>%
  filter(tmr_region != 'mucosa') %>%
  group_by(Tissue_ID, sample_tmr_region) %>%
  summarise(across(genesets,
    \(x) median(x, na.rm = TRUE)
  ), .groups = "drop")

scaled_CRC_median <- paired_CRC_median %>%
  group_by(Tissue_ID) %>%
  mutate(across(
    .cols = -all_of(c("sample_tmr_region")),
    .fns = ~ as.numeric(scale(.)),
    .names = "scaled_{.col}"
  )) %>%
  ungroup()

scaled_matrix <- scaled_CRC_median %>%
  select(sample_tmr_region, starts_with("scaled_")) %>%
  tibble::column_to_rownames("sample_tmr_region") %>%
  as.matrix()

# Transpose to get gene sets as rows, samples as columns
heat_mat <- t(scaled_matrix)

# Extract TMR region (submucosa, muscle, metastasis) from column names
tmr_region_vector <- sub(".*_", "", colnames(heat_mat))
names(tmr_region_vector) <- colnames(heat_mat)

# Order columns by desired TMR region order
region_order <- c("submucosa", "muscle", "metastasis")
col_order <- unlist(lapply(region_order, function(r) which(tmr_region_vector == r)))

# Apply ordering
heat_mat_ordered <- heat_mat[, col_order]
tmr_region_vector_ordered <- tmr_region_vector[col_order]

# Create top annotation
column_ha <- HeatmapAnnotation(
  Region = tmr_region_vector_ordered,
  col = list(
    Region = tmr_col),
  show_annotation_name = TRUE
)


piyg_colors <- rev(brewer.pal(11, "PiYG"))

# Create color mapping (e.g., for Z-scores)
col_fun <- colorRamp2(c(-2, 0, 2), c(piyg_colors[1], piyg_colors[6], piyg_colors[11]))

                           
# Draw heatmap with column split
Heatmap(
  heat_mat_ordered,
  name = "Z-score",
  col = col_fun,
  column_split = factor(tmr_region_vector_ordered, levels = region_order),
  top_annotation = column_ha,
  column_names_rot = 45,
  column_title = "Samples (Grouped by TMR Region)",
  row_title = "Gene Signatures",
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE
) -> heatP1


pdf(file=file.path(out_dir, 'mCRC_N26_TMR_tumor_heatmap_by_tmr_regions.pdf'), width=8, height=3) 
draw(heatP1)
dev.off()

# Extract Case ID and TMR region
case_id_vector <- sub("^(CM\\d+C\\d+)_.*", "\\1", colnames(heat_mat_ordered))
tmr_region_vector_ordered <- sub(".*_", "", colnames(heat_mat_ordered))

# Column annotation (by region)
column_ha <- HeatmapAnnotation(
  Region = tmr_region_vector_ordered,
  col = list(
    Region = c(
      submucosa = "#4DBBD5",
      muscle = "#00A087",
      metastasis = "#E64B35"
    )
  ),
  show_annotation_name = TRUE
)

# Draw heatmap split by Case ID
Heatmap(
  heat_mat_ordered,
  name = "Z-score",
  col = col_fun,
  column_split = factor(case_id_vector, levels = unique(case_id_vector)),
  top_annotation = column_ha,
  column_names_rot = 45,
  column_title = "Samples (Grouped by Case)",
  row_title = "Gene Signatures",
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = TRUE,
  show_row_names = TRUE
) -> heatP2

pdf(file=file.path(out_dir, 'mCRC_N26_TMR_tumor_heatmap2_by_tmr_regions.pdf'), width=8, height=4) 
draw(heatP2)
dev.off()