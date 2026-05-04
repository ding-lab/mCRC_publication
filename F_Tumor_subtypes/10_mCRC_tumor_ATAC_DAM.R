# Differential Accessible Motifs (DAM) Analysis for mCRC Tumor Subtypes
# Analysis of transcription factor motif accessibility across tumor cell subtypes
# using chromVAR data from snATAC-seq

# Load required libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(TFBSTools)
library(JASPAR2020)
library(EnhancedVolcano)
library(GSEABase)
library(ChIPseeker)
library(ensembldb)
library(EnsDb.Hsapiens.v100)
library(BSgenome.Hsapiens.UCSC.hg38)
library(RColorBrewer)
library(ComplexHeatmap)
library(rstatix)
library(ggpubr)
library(patchwork)
library(scCustomize)
library(chromVAR)
library(SummarizedExperiment)
library(BiocParallel)
library(motifmatchr)
source('/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R')
set.seed(1234)

# Get JASPAR motif database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Create motif to transcription factor mapping
motif.tf <- data.frame(motif = names(pfm@listData),
           TF = purrr::map_chr(names(pfm@listData), ~pfm@listData[[.x]]@name))

# Helper functions for motif-TF conversion
getTF <- function(motif.code) {
    motif.tf %>% filter(motif==motif.code) %>% pull(TF)
}

getMotif <- function(tf) {
    motif.tf %>% filter(TF==tf) %>% pull(motif)
}

# Load tumor ATAC-seq object with chromVAR data
rds_path <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/ATAC/tumor/'
tumor_obj <- readRDS(file.path(rds_path, 'chromvar', '41_merged_normalized_mCRC_Combo_clean_tumor.chromvar.rds'))

# Set output directory
output_dir <- getwd()

# Load and add metadata
metadata_path <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/metadata/RNA/'
metadata <- read.csv(file.path(metadata_path, 'mCRC_57_samples_clean3_metadata_cell_type_all_20250505.csv'), header=TRUE)
rownames(metadata) <- metadata$X
metadata$X <- NULL
tumor_obj <- AddMetaData(tumor_obj, metadata)

# Calculate average chromVAR accessibility per cell type
DefaultAssay(tumor_obj) <- 'chromvar'
tumor.chromvar.cell.type <- AverageExpression(tumor_obj, 
                                              group.by = 'cell_type_xenium', 
                                              assays = 'chromvar', 
                                              return.seurat = TRUE, 
                                              layer = 'counts')
tumor.chromvar.cell.type$orig.ident <- rownames(tumor.chromvar.cell.type@meta.data)

# Perform differential accessible motifs (DAM) analysis
Idents(tumor_obj) <- 'cell_type_xenium'

tumor_dam <- FindAllMarkers(object = tumor_obj,
                          logfc.threshold = 0.1,
                          only.pos = FALSE,
                          mean.fxn = rowMeans,
                          fc.name = "avg_diff")

# Add motif and TF information to DAM results
tumor_dam$motif <- rownames(tumor_dam)
tumor_dam$TF <- motif.tf$TF[match(tumor_dam$gene, motif.tf$motif)]

# Save DAM results
write.csv(tumor_dam, file.path(rds_path, 'mCRC_tumor_dam.csv'), row.names = FALSE)

# Example: Check specific TF (PROX1)
tumor_dam %>% filter(TF == 'PROX1')

# Identify highly significant TFs
interesting.tfs <- tumor_dam %>%  
  filter(abs(avg_diff) > 1 & p_val_adj < 1e-50) %>% 
  pull(gene)
length(interesting.tfs)

# Prepare data for heatmap visualization
toplot <- FetchData(tumor.chromvar.cell.type, 
                   vars = interesting.tfs, 
                   slot = 'counts')
colnames(toplot) <- motif.tf$TF[match(colnames(toplot), motif.tf$motif)]

# Define color palette for heatmap
breaks <- seq(-2, 2, by=0.001)
color.palette <- circlize::colorRamp2(
  breaks = breaks, 
  colors = colorRampPalette(c('#2b2d42', '#8d99ae', '#edf2f4', '#ff4d6d', '#d80032'))(length(breaks))
)

# Create and save heatmap
p <- Heatmap(toplot, 
             name = "Motif accessibility", 
             col = color.palette, 
             show_column_names = TRUE,
             clustering_method_columns = 'ward.D2',
             clustering_method_rows = 'ward.D2')

pdf(file.path(output_dir, 'Heatmap_mCRC_tumor_DAMs.pdf'), width = 24, height = 4)
draw(p, padding = unit(c(5, 5, 5, 10), "mm"))
dev.off()

# Function for pairwise cluster comparison and volcano plot
compare_clusters_and_plot <- function(
  seurat_obj, 
  ident1, 
  ident2, 
  logfc_threshold = 0.1, 
  fc_name = "avg_diff", 
  mean_function = rowMeans, 
  fc_cutoff = 0.25, 
  xlim_values = c(-2, 2),
  p_Cutoff = 1e-05,  
  title_text = 'Cluster Comparison',
  motif_tf_data
) {
  # Perform differential expression analysis
  markers <- FindMarkers(
    object = seurat_obj,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = logfc_threshold,
    only.pos = FALSE,
    mean.fxn = mean_function,
    pCutoff = p_Cutoff,  
    fc.name = fc_name
  )
  
  # Add motif and TF information
  markers$motif <- rownames(markers)
  markers$TF <- motif_tf_data$TF[match(markers$motif, motif_tf_data$motif)]
  
  # Create the volcano plot
  volcano_plot <- EnhancedVolcano(
    markers,
    lab = markers$TF,
    x = fc_name,
    y = 'p_val_adj',
    xlab = 'average difference', 
    xlim = xlim_values,        
    FCcutoff = fc_cutoff,
    pCutoff = p_Cutoff,  
    title = title_text,
    drawConnectors = TRUE,
    widthConnectors = 0.75,  
    subtitle = '',  
    titleLabSize = 14,  
    legendLabels = c("NS", "average difference", "p-value", 
                     expression(p - value ~ and ~ 'average difference'))
  )
  
  return(volcano_plot)
}

# Example: Compare Non_Canonical_CRC_1 vs Canonical_CRC_Stem
volcano_plot <- compare_clusters_and_plot(
  seurat_obj = tumor_obj, 
  ident1 = 'Non_Canonical_CRC_1', 
  ident2 = 'Canonical_CRC_Stem', 
  motif_tf_data = motif.tf, 
  xlim_values = c(-5, 5),
  p_Cutoff = 5e-02,    
  title_text = 'Non_Canonical_CRC versus Canonical-CRC-Stem'
)

# Save volcano plot
pdf(file.path(output_dir, 'Volcano_Non_canonical_vs_Stem_DAMs.pdf'), width = 12, height = 8)
print(volcano_plot)
dev.off()

