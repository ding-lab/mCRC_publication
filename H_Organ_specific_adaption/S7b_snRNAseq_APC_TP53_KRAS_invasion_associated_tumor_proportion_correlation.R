library(Seurat)
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)
library(cowplot)
library(clinfun)
library(broom)
library(car)     
library(emmeans)
library(scales)

output_dir = getwd()

rds_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/'
tumor_obj = readRDS(file.path(rds_path,'57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean5.rds'))
tumor_obj

metadata_path = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata = read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL

tumor_obj <- AddMetaData(tumor_obj, metadata)