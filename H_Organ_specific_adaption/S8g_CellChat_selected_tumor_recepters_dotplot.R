suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(ggpubr))
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/mCRC_colors.R")
output_dir = getwd()

tumor_rds_name = '57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean5.rds'
tumor_rds_path = file.path('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/', 
                           tumor_rds_name)
tumor_obj = readRDS(tumor_rds_path)

metadata_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata = read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL
tumor_obj <- AddMetaData(tumor_obj, metadata)

tumor_cell_obj <- tumor_obj %>% 
                 subset(Site_of_Origin %in% c('colon', 'rectum', 'liver', 'lung')) %>% 
                 subset(Patient_ID != 'HT413C1') %>% 
                 subset(!(cell_type_xenium %in% c('APCDD1_CRC', 'Canonical_CRC_Stem_Proliferation')))

# Create Organ variable combining colon and rectum into colorectum
tumor_cell_obj@meta.data <- tumor_cell_obj@meta.data %>% 
                            mutate(Organ = case_when(Site_of_Origin %in% c('rectum', 'colon') ~ 'colorectum',
                                                     TRUE ~ Site_of_Origin))

tumor_cell_obj$cell_type_xenium <- factor(tumor_cell_obj$cell_type_xenium,
                                          level = c('Non_Canonical_CRC_1',
                                                    'Canonical_CRC_Stem',
                                                    'Canonical_CRC_Intestine_Proliferation',
                                                    'Canonical_CRC_Intestine'
                                                   ))

# Create combined grouping variable: tumor_type_organ
tumor_cell_obj@meta.data <- tumor_cell_obj@meta.data %>% 
                            mutate(tumor_type_organ = paste0(Organ, '_', cell_type_xenium))

# Define factor levels for the 12 groups (4 tumor types × 3 organs)
tumor_type_levels <- c('Non_Canonical_CRC_1',
                       'Canonical_CRC_Stem',
                       'Canonical_CRC_Intestine_Proliferation',
                       'Canonical_CRC_Intestine')
organ_levels <- c('colorectum', 'liver', 'lung')
tumor_type_organ_levels <- paste0(rep(organ_levels, each = length(tumor_type_levels)), '_', 
                                  rep(tumor_type_levels, times = length(organ_levels)))

tumor_cell_obj$tumor_type_organ <- factor(tumor_cell_obj$tumor_type_organ,
                                          levels = tumor_type_organ_levels)

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)

p = DotPlot(tumor_cell_obj, features=rev(c('MET', 'EGFR', 'ERBB2', 'ROBO1', 'IGF1R')) , group.by='tumor_type_organ') +
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      coord_flip() +
      #coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-2, 2),   
                            breaks = color_breaks) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)


pdf(file.path(output_dir, 'EGFR_ERBB2_MET_ROBO1_IGFR1_gene_expression_by_tumor_subtypes.pdf'), width = 7, height = 12)
print(p)
dev.off()