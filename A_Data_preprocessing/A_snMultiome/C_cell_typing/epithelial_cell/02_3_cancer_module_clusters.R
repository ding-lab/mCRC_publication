library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(MAST)
library(readr)
library(RColorBrewer)
library(scales)
library(rstatix)
library(data.table)
library(googlesheets4)
source("/diskmnt/Projects/Users/Evan.p/scripts/Rscript/jupyter_support_functions.R")
source("/diskmnt/Projects/Users/Evan.p/scripts/Rscript/visiualization_support_functions.R")

geneset_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/genesets'
colon_epithelial_geneset = readRDS(file.path(geneset_dir, 'colon_epithelium_genesets.rds'))
module_geneset = readRDS(file.path(geneset_dir, 'gene_modules.rds'))
Hallmark_EMT = readRDS(file.path(geneset_dir, 'Hallmark_EMT_genesets.rds'))
Hallmark_hypoxia = readRDS(file.path(geneset_dir, 'Hallmark_hpoxia_genesets.rds'))
Kohn_EMT = readRDS(file.path(geneset_dir, 'Kohn_Emt_genesets.rds'))
SenMayo = readRDS(file.path(geneset_dir, 'SenMayo_genesets.rds'))
module_geneset <- c(module_geneset, Hallmark_hypoxia)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial')
epithelial_reint = readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds')
epithelial_reint

colnames(epithelial_reint@meta.data)
unique(epithelial_reint$epi_cell_type)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epithelial_cluster_sub6', label = TRUE)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type', label = TRUE)

Idents(epithelial_reint) = "epi_cell_type"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 'Stem-like tumor',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epi_cell_type_sub0",
                                  resolution = 0.3,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epi_cell_type_sub0)

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type_sub0', label = TRUE)

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type_sub0', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 'SLC6A6', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/subcluster_umaps')
# Loop through all unique clusters in the epi_cell_type_sub0 metadata column
for (cluster in unique(epithelial_reint$epi_cell_type_sub0)) {
  
  # Print the current cluster for debugging purposes
  cat("Processing cluster:", cluster, "\n")
  
  # Create the file name for saving the plot
  file_name <- paste0("UMAP_Highlight_", gsub(" ", "_", cluster), ".png")
  
  # Generate the UMAP plot for the current cluster
  p = Highlight_Cluster_UMAP(
    seurat_object = epithelial_reint,
    metadata_column = 'epi_cell_type_sub0',
    cell_type = cluster,  # Pass the cluster name dynamically
    umap_reduction = 'epithelial_umap.scvi',
    highlight_color = "#B56727", 
    background_color = "lightgray", 
    background_alpha = 0.3
  )
    
  ggsave(
  filename = paste0("Highlight_UMAP_", gsub(" ", "_", cluster), ".jpg"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300)
}


epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type2 = case_when(
                                  epi_cell_type_sub0 == 'Stem-like tumor_0' ~ 'Intestine-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_1' ~ 'Tumor-ISC-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_2' ~ 'Tumor-ISC-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_3' ~ 'Tumor-ISC-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_4' ~ 'Intestine-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_5' ~ 'Secretory-intestine',
                                  epi_cell_type_sub0 == 'Stem-like tumor_6' ~ 'Intestine-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_7' ~ 'Intestine-like',
                                  epi_cell_type_sub0 == 'Stem-like tumor_8' ~ 'Intestine-like',
                                  TRUE ~ epi_cell_type
                              ))

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type2', label = TRUE)

Idents(epithelial_reint) = "epi_cell_type2"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 'Secretory-intestine',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epi_cell_type2_0",
                                  resolution = 0.1,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epi_cell_type2_0)

set_size(8,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type2_0', label = TRUE)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/subcluster_umaps')
# Loop through all unique clusters in the epi_cell_type_sub0 metadata column
for (cluster in unique(epithelial_reint$epi_cell_type2_0)) {
  
  # Print the current cluster for debugging purposes
  cat("Processing cluster:", cluster, "\n")
  
  # Create the file name for saving the plot
  file_name <- paste0("UMAP_Highlight_", gsub(" ", "_", cluster), ".png")
  
  # Generate the UMAP plot for the current cluster
  p = Highlight_Cluster_UMAP(
    seurat_object = epithelial_reint,
    metadata_column = 'epi_cell_type2_0',
    cell_type = cluster,  # Pass the cluster name dynamically
    umap_reduction = 'epithelial_umap.scvi',
    highlight_color = "#B56727", 
    background_color = "lightgray", 
    background_alpha = 0.3
  )
    
  ggsave(
  filename = paste0("Highlight_UMAP_", gsub(" ", "_", cluster), ".jpg"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300)
}


epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type3 = case_when(
                                  epi_cell_type2_0 == 'Secretory-intestine_1' ~ 'Intestine-like',
                                  epi_cell_type2_0 == 'Secretory-intestine_0' ~ 'Secretory-intestine',
                                  TRUE ~ epi_cell_type2_0
                              ))

set_size(8,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type3', label = TRUE)

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type3', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 'SLC6A6', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

Idents(epithelial_reint) = "epi_cell_type3"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 'Angiogenic tumor',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epi_cell_type3_0",
                                  resolution = 0.3,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epi_cell_type3_0)

set_size(8,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type3_0', label = TRUE)

table(epithelial_reint$epi_cell_type3_0)

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type3_0', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 'SLC6A6', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/subcluster_umaps')
# Loop through all unique clusters in the epi_cell_type_sub0 metadata column
for (cluster in unique(epithelial_reint$epi_cell_type3_0)) {
  
  # Print the current cluster for debugging purposes
  cat("Processing cluster:", cluster, "\n")
  
  # Create the file name for saving the plot
  file_name <- paste0("UMAP_Highlight_", gsub(" ", "_", cluster), ".png")
  
  # Generate the UMAP plot for the current cluster
  p = Highlight_Cluster_UMAP(
    seurat_object = epithelial_reint,
    metadata_column = 'epi_cell_type3_0',
    cell_type = cluster,  # Pass the cluster name dynamically
    umap_reduction = 'epithelial_umap.scvi',
    highlight_color = "#B56727", 
    background_color = "lightgray", 
    background_alpha = 0.3
  )
    
  ggsave(
  filename = paste0("Highlight_UMAP_", gsub(" ", "_", cluster), ".jpg"),
  plot = p,
  width = 6,
  height = 6,
  dpi = 300)
}

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type4 = case_when(
                                  epi_cell_type3_0 == 'Angiogenic tumor_0' ~ 'Hypoxic-Repair',
                                  epi_cell_type3_0 == 'Angiogenic tumor_1' ~ 'Hypoxic',
                                  epi_cell_type3_0 == 'Angiogenic tumor_2' ~ 'Hypoxic-EMT',
                                  epi_cell_type3_0 == 'Angiogenic tumor_3' ~ 'Hypoxic',
                                  epi_cell_type3_0 == 'Angiogenic tumor_4' ~ 'Hypoxic-Endoderm',
                                  epi_cell_type3_0 == 'Angiogenic tumor_5' ~ 'Hypoxic',
                                  epi_cell_type3_0 == 'Angiogenic tumor_6' ~ 'Hypoxic',
                                  epi_cell_type3_0 == 'Angiogenic tumor_7' ~ 'Hypoxic-Repair',
                                  epi_cell_type3_0 == 'Angiogenic tumor_8' ~ 'Hypoxic-Repair',
                                  epi_cell_type3_0 == 'Angiogenic tumor_9' ~ 'Hypoxic',
                                  epi_cell_type3_0 == 'Angiogenic tumor_10' ~ 'Hypoxic',
                                  TRUE ~ epi_cell_type3
                              ))

epi.cell.color <- c("APCDD1+ tumor" = '#008856',
                    "Hypoxic" = '#be0032',
                    "Hypoxic-Endoderm" = '#dcd300',
                    "Hypoxic-Repair" = '#b3446c',
                    "Hypoxic-EMT" = '#e25822',
                    "Proliferative tumor" = '#f3c300',
                    "Proliferative stem-like tumor" = '#8db600',
                    "Tumor-ISC-like" = '#e68fac',
                    "Intestine-like" = '#f99379',
                    "Stem cells" = '#e25822',
                    "Enterocytes" = '#882d17',
                    "Transit-amplifying cells" = '#a1caf1',
                    "Goblet cells" = '#c2b280',
                    "Tuft cells" = '#f99379',
                    "Enteroendocrine-like cells" = '#2b3d26'
                   )

set_size(6,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', 
        group.by = 'epi_cell_type4', 
        label = TRUE) + scale_color_manual(values=epi.cell.color, name='cell type')

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type4', 
              feature = module_geneset$`Injury Repair`
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(30,8)
pdf("DotPlot_Injury_Repair.pdf", width=24, height=12)
print(p6)
dev.off()

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type4', 
              feature = c('APCDD1', 'BAMBI', 'VEGFA', 'SLC2A1', 'BRCA2', 'MKI67', 'TOP2A', 'LGR5', 'SMOC2', 'OLFM4', 'RGMB', 'SLC6A6', 
                          'CLCA4', 'SLC26A3', 'TICRR', 'CDC25C', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

Idents(epithelial_reint) = "epi_cell_type4"
epithelial_reint = FindSubCluster(object = epithelial_reint,
                                  cluster = 'APCDD1+ tumor',
                                  graph.name = 'RNA_snn',
                                  subcluster.name = "epi_cell_type4_0",
                                  resolution = 0.2,
                                  algorithm = 1
                                  )
unique(epithelial_reint$epi_cell_type4_0)

set_size(8,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', group.by = 'epi_cell_type4_0', label = TRUE)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type5 = case_when(
                                  epi_cell_type4_0 == 'APCDD1+ tumor_0' ~ 'APCDD1 tumor',
                                  epi_cell_type4_0 == 'APCDD1+ tumor_1' ~ 'APCDD1 tumor',
                                  epi_cell_type4_0 == 'APCDD1+ tumor_2' ~ 'APCDD1-Neuroendocrine',
                                  epi_cell_type4_0 == 'APCDD1+ tumor_3' ~ 'APCDD1 tumor',
                                  TRUE ~ epi_cell_type4
                              ))

epithelial_reint$epi_cell_type5 <- factor(epithelial_reint$epi_cell_type5, 
                                         levels = c('APCDD1 tumor', "APCDD1-Neuroendocrine", 
                                                    "Hypoxic", "Hypoxic-Endoderm", "Hypoxic-Repair", "Hypoxic-EMT",   
                                                    'Proliferative tumor', "Intestine-like", "Secretory-intestine",
                                                    "Tumor-ISC-like", 'Proliferative stem-like tumor',
                                                    "Stem cells",
                                                    'Enterocytes', 'Transit-amplifying cells', 
                                                    'Goblet cells', 'Tuft cells', 'Enteroendocrine-like cells'
                                                   ))

epi.cell.color2 <- c("APCDD1 tumor" = '#008856',
                    "APCDD1-Neuroendocrine" = '#604197',
                    "Hypoxic" = '#be0032',
                    "Hypoxic-Endoderm" = '#dcd300',
                    "Hypoxic-Repair" = '#b3446c',
                    "Hypoxic-EMT" = '#e25822',
                    "Proliferative tumor" = '#f3c300',
                    "Proliferative stem-like tumor" = '#8db600',
                    "Tumor-ISC-like" = '#e68fac',
                    "Intestine-like" = '#f99379',
                    "Secretory-intestine" = 'firebrick',
                    "Stem cells" = 'khaki',
                    "Enterocytes" = '#882d17',
                    "Transit-amplifying cells" = '#a1caf1',
                    "Goblet cells" = '#c2b280',
                    "Tuft cells" = '#f99379',
                    "Enteroendocrine-like cells" = '#2b3d26'
                   )

set_size(7,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', 
        group.by = 'epi_cell_type5', label.size = 5,
        label = TRUE) + scale_color_manual(values=epi.cell.color2, name='cell type')

pdf("Dimplot_mCRC_epithelial_cell_type5.pdf", width=7, height=6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', 
        group.by = 'epi_cell_type5', label.size = 5,
        label = FALSE) + scale_color_manual(values=epi.cell.color2, name='cell type')
dev.off()

set_size(9,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', 
        group.by = 'orig.ident', 
        label = TRUE) 

colnames(epithelial_reint@meta.data)

output_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/snRNA/01_Cohort/Metadata/'
celltype_update <- epithelial_reint@meta.data %>% select(epi_cell_type, epi_cell_type5)
fwrite(celltype_update, file.path(output_dir, "57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean5_epithelial_cell_type_20250126.csv"), row.names=TRUE)

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial')
saveRDS(epithelial_reint, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds')

setwd('/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial')
epithelial_reint <- readRDS('57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds')
epithelial_reint

colnames(epithelial_reint@meta.data)

unique(epithelial_reint$epi_cell_type6)

output_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/snRNA/01_Cohort/Metadata/'
celltype_update <- epithelial_reint@meta.data %>% select(epi_cell_type, epi_cell_type6)
fwrite(celltype_update, file.path(output_dir, "57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean5_epithelial_cell_type_20250322.csv"), row.names=TRUE)

Idents(epithelial_reint) = "epi_cell_type6"
epithelial_markers = FindAllMarkers(epithelial_reint, 
                                    min.pct=0.1, 
                                    logfc.threshold=0.25, 
                                    only.pos = TRUE, 
                                    return.thresh = 0.01,
                                    test.use = 'MAST'
                                    )

write.csv(epithelial_markers, "mCRC_epithelial_type6_MAST_deg.csv", row.names = TRUE)

tumor <- epithelial_reint %>% 
         subset(epi_cell_type6 %in% c('Canonical_CRC_Intestine_Proliferation', 
                                      'Canonical_CRC_Stem_Proliferation', 
                                      'Canonical_CRC_Intestine', 
                                      'Canonical_CRC_Stem',
                                      'Non_Canonical_CRC_2',
                                      'Non_Canonical_CRC_1'))
tumor

Idents(tumor) = "epi_cell_type6"
tumor_markers = FindAllMarkers(tumor, 
                               min.pct=0.1, 
                               logfc.threshold=0.25, 
                               only.pos = TRUE, 
                               return.thresh = 0.01,
                               test.use = 'MAST')

write.csv(tumor, "mCRC_tumor_MAST_deg.csv", row.names = TRUE)

epithelial_reint@meta.data <- epithelial_reint@meta.data %>% 
                              mutate(epi_cell_type6 = case_when(
                                  epi_cell_type5 == 'Tumor-ISC-like' ~ 'Canonical_CRC_Stem',
                                  epi_cell_type5 == 'Intestine-like' ~ 'Canonical_CRC_Intestine',
                                  epi_cell_type5 == 'Secretory-intestine' ~ 'Canonical_CRC_Intestine',
                                  epi_cell_type5 == 'Proliferative tumor' ~ 'Canonical_CRC_Intestine_Proliferation',
                                  epi_cell_type5 == 'Proliferative stem-like tumor' ~ 'Canonical_CRC_Stem_Proliferation',
                                  epi_cell_type5 == 'Hypoxic' ~ 'Non_Canonical_CRC_1',
                                  epi_cell_type5 == 'Hypoxic-Repair' ~ 'Non_Canonical_CRC_1',
                                  epi_cell_type5 == 'Hypoxic-EMT' ~ 'Non_Canonical_CRC_1',
                                  epi_cell_type5 == 'Hypoxic-Endoderm' ~ 'Non_Canonical_CRC_1',
                                  epi_cell_type5 == 'APCDD1 tumor' ~ 'Non_Canonical_CRC_2',
                                  epi_cell_type5 == 'APCDD1-Neuroendocrine' ~ 'Non_Canonical_CRC_2',
                                  TRUE~epi_cell_type5
             ))

epithelial_reint$epi_cell_type6 <- factor(epithelial_reint$epi_cell_type6, 
                                         levels = rev(c('Non_Canonical_CRC_1', 'Non_Canonical_CRC_2',
                                                        'Canonical_CRC_Stem', 'Canonical_CRC_Intestine', 
                                                        'Canonical_CRC_Stem_Proliferation', 
                                                        'Canonical_CRC_Intestine_Proliferation', 
                                                        'Stem cells', 'Transit-amplifying cells', 'Enterocytes',
                                                        'Goblet cells', 'Tuft cells', 'Enteroendocrine-like cells'
                                                   )))

DefaultAssay(epithelial_reint) <- 'RNA'

rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))
color_breaks <- c(-2, -1, 0, 1, 2)
#dotplot.color <- colorRampPalette(c('#eae2b7','#fcbf49','#f77f00','#d62828', '#003049'))(10) 

p6 <- DotPlot(epithelial_reint, 
              group.by = 'epi_cell_type6', 
              feature = c('EMP1', 'KRT20', 'VEGFA', 'APCDD1', 'PROX1', 'LGR5', 'RGMB', 'MKI67', 'TOP2A', 
                          'CLCA4', 'SLC26A3', 'MUC2', 'CLCA1', 'SH2D6', 'TRPM5', 'CHGA', 'CHGB')
             ) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + #, legend.byrow = FALSE, legend.direction = 'vertical', legend.position = 'bottom') +
      coord_fixed() +
      scale_color_gradientn(colors = rdylbu_colors, 
                            limits = c(-3, 3),    # Set limits to -2 and 2
                            breaks = color_breaks # Specify the breaks
      ) +         
      scale_size_area(limits = c(0, 100), oob = scales::squish)

set_size(6,4)
p6

output_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Manuscript/Figure_3'
pdf(file.path(output_dir, "Dotplot_mCRC_epithelial_cell_type_snRNA.pdf"), width=8, height=8)
print(p6)
dev.off()

epi.cell.color2 <- c("Non_Canonical_CRC_2" = '#008856',
                     "Non_Canonical_CRC_1" = '#be0032',
                     "Canonical_CRC_Intestine" = '#e25822',
                     "Canonical_CRC_Stem_Proliferation" = '#f3c300',
                     "Canonical_CRC_Intestine_Proliferation" = '#8db600',
                     "Canonical_CRC_Stem" = '#e68fac',
                     "Stem cells" = '#604197',
                     "Enterocytes" = '#882d17',
                     "Transit-amplifying cells" = '#a1caf1',
                     "Goblet cells" = '#c2b280',
                     "Tuft cells" = '#f99379',
                     "Enteroendocrine-like cells" = '#2b3d26'
                   )

epithelial_reint$epi_cell_type6 <- factor(epithelial_reint$epi_cell_type6, 
                                         levels = c('Non_Canonical_CRC_1', 'Non_Canonical_CRC_2',
                                                        'Canonical_CRC_Stem', 'Canonical_CRC_Intestine', 
                                                        'Canonical_CRC_Stem_Proliferation', 
                                                        'Canonical_CRC_Intestine_Proliferation', 
                                                        'Stem cells', 'Transit-amplifying cells', 'Enterocytes',
                                                        'Goblet cells', 'Tuft cells', 'Enteroendocrine-like cells'
                                                   ))

set_size(7,6)
DimPlot(epithelial_reint, reduction = 'epithelial_umap.scvi', 
        group.by = 'epi_cell_type6', label.size = 5,
        label = FALSE) + scale_color_manual(values=epi.cell.color2, name='cell type') -> p2

p2

output_dir = '/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/Manuscript/Figure_3'
pdf(file.path(output_dir, "Dimplot_mCRC_epithelial_cell_type_snRNA.pdf"), width=8, height=6)
print(p2)
dev.off()


