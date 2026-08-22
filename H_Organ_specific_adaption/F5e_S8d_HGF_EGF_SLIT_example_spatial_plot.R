suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(qs))
suppressMessages(library(future))
options(future.globals.maxSize = 30 * 1024^3)
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Plot/ImagePlotFaster_v1.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool_github/ImageBinPlot/R/ImageBinPlot.R')
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")
Xenium_5K_path = '/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_5K/'
output_dir = getwd()

xenium_annotation_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/'
xenium_annotation_path = file.path(xenium_annotation_dir, 'mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv')
xenium_annotation = read.csv(xenium_annotation_path, header = TRUE) %>% 
                    mutate(GTR = if_else(tns_label > 0, 'GTR', 'NAT')) %>% 
                    mutate(TMR = if_else(tn_label > 0, 'TMR', 'non-TMR')) %>% 
                    mutate(GTR_cat = case_when(
                        tn_label > 0 ~ 'TMR',
                        tns_label == 0 ~ 'NAT',
                        TRUE ~ 'TSR'
                    )) %>% 
                    mutate(GTR_NB = case_when(
                        GTR == 'NAT' ~ 'NAT',
                        TRUE ~ neighborhoods
                    )) %>% 
                    mutate(tns_mask = case_when(
                           tns_label > 0 ~ 'GTR',
                           TRUE ~ 'APT')) %>% 
                    mutate(GTR_cat2 = case_when(
                           tns_label == 0 ~ 'APT',
                           tn_label == 0 ~ 'PSIZ',
                           t_label > 0 ~ 'Tumor',
                           TRUE ~ 'Necrosis')) %>% 
                    mutate(TMR_tumor_cat = case_when(
                           tn_label > 0 & All_cell_type1 == 'Non-canonical' ~ 'Pro-invasive',
                           tn_label > 0 & All_cell_type1 %in% c('Proliferative-like', 'Stem-like', 'Intestine-like') ~ 'Canonical',
                           TRUE ~ 'Other'))


NB_colors = c('NB1_necrosis' = '#604e97', 'NB2_tumor_core' = 'hotpink4', 
              'NB3_tumor_body' = 'violetred', 'NB4_tumor_stroma_interface' = 'orange1',
              'NB5_desmoplastic_stroma' = 'darkgreen', 'NB6_perivascular_stroma' = 'yellowgreen',
              'NB7_myelocyte_enriched_stroma' = '#be0032', 'NB8_lymphocyte_enriched_stroma' = '#dcd300', 
              'NB9_plasma_cell_enriched_stroma' = '#0067a5', 'NB10_smooth_muscle_stroma' = 'lightseagreen',
              'NB11_colon_epithelium' = 'gold2', 'NB12_germinal_center' = 'yellow4', 'NB13_liver' = '#882d17', 
              'NB14_ductal_epithelium' = '#2b3d26',  'Unassined' = KellyPalette$grey)
GTR_NB_colors = NB_colors
GTR_NB_colors['NAT'] <- '#7f7e80'
TMR_tumor_cat_colors = c('Pro-invasive' = "#be0032", 'Canonical' = "darkseagreen", 'Other' = '#7f7e80')

gtr_col <- list(
    GTR_cat = c(
        'NAT' = "grey80",  
        'TSR' = "palegreen4",
        'TMR' = "firebrick2"),
    neighborhoods = NB_colors,
    GTR_NB = GTR_NB_colors,
    tns_mask = c('GTR' = 'steelblue1', 'APT' = KellyPalette$grey),
    GTR_cat2 = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey)
)

CM478C2_A3_path = file.path(Xenium_5K_path, 'CM478C2-A3-Tp1Fp1U1/CM478C2-A3-Tp1Fp1U1_processed.rds')
CM1799C1_A1_path = file.path(Xenium_5K_path, 'CM1799-A1-Th1Fp1U1/CM1799-A1-Th1Fp1U1_processed.rds')
CM1691C1_A3_path = file.path(Xenium_5K_path, 'CM1691-A3-Th1Fp1U1/CM1691-A3-Th1Fp1U1_processed.rds')

# CM478C2-A3

CM478C2_A3 = readRDS(CM478C2_A3_path)
CM478C2_A3_meta = xenium_annotation %>% filter(sample_id == 'CM478C2-A3-Tp1Fp1U1')
rownames(CM478C2_A3_meta) <- CM478C2_A3_meta$cell_id
CM478C2_A3 <- AddMetaData(CM478C2_A3, CM478C2_A3_meta)
DefaultAssay(CM478C2_A3) <- 'Xenium'
DefaultBoundary(CM478C2_A3[["fov"]]) <- 'centroids'

# SLIT2 count 14112 on Xenium Explorer
# ROBO1 count 51219 on Xenium Explorer
p1 = ImageDimPlot(CM478C2_A3, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('ROBO1', 'SLIT2'),
                  mols.col = c('yellow', 'cyan'),
                  mols.size = 0.3,
                  nmols = 1411,
                  coord.fixed = TRUE)


ggplot2::ggsave(
  file.path(output_dir, "M478C2-A3-Tp1Fp1U1_ROBO1_SLIT2_spatial_plot.png"),
  plot = p1, width = 6, height = 6, units = "in", dpi = 300
)


# HGF count 12532 on Xenium Explorer
# MET count 90458 on Xenium Explorer
p2 = ImageDimPlot(CM478C2_A3, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('MET', 'HGF'),
                  mols.col = c('green', 'orange'),
                  nmols = 1253,
                  mols.size = 0.3,
                  coord.fixed = TRUE)
ggplot2::ggsave(
  file.path(output_dir, "M478C2-A3-Tp1Fp1U1_MET_HGF_spatial_plot.png"),
  plot = p2, width = 6, height = 6, units = "in", dpi = 300
)

# HBEGF count 10589 on Xenium Explorer
# EGFR count 24468 on Xenium Explorer

p7 = ImageDimPlot(CM478C2_A3, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('HBEGF', 'EGFR'),
                  mols.col = c('gold', 'green'),
                  nmols = 1058,
                  mols.size = 0.3,
                  coord.fixed = TRUE)

ggplot2::ggsave(
  file.path(output_dir, "M478C2-A3-Tp1Fp1U1_HBEGF_EGFR_spatial_plot.png"),
  plot = p7, width = 6, height = 6, units = "in", dpi = 300
)


#CM1799C1-A1
CM1799C1_A1 = readRDS(CM1799C1_A1_path)
CM1799C1_A1_meta = xenium_annotation %>% filter(sample_id == 'CM1799-A1-Th1Fp1U1')
rownames(CM1799C1_A1_meta) <- CM1799C1_A1_meta$cell_id
CM1799C1_A1 <- AddMetaData(CM1799C1_A1, CM1799C1_A1_meta)

DefaultAssay(CM1799C1_A1) <- 'Xenium'
DefaultBoundary(CM1799C1_A1[["fov"]]) <- 'centroids'

# SLIT2 count 9604 on Xenium Explorer
# ROBO1 count 44517 on Xenium Explorer

p3 = ImageDimPlot(CM1799C1_A1, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('ROBO1', 'SLIT2'),
                  mols.col = c('yellow', 'cyan'),
                  nmols = 960,
                  mols.size = 0.3,
                  coord.fixed = TRUE)

ggplot2::ggsave(
  file.path(output_dir, "CM1799-A1-Th1Fp1U1_ROBO1_SLIT2_spatial_plot.png"),
  plot = p3, width = 6, height = 6, units = "in", dpi = 300
)


# HGF count 19604 on Xenium Explorer
# MET count 96722 on Xenium Explorer

p5 = ImageDimPlot(CM1799C1_A1, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('MET', 'HGF'),
                  mols.col = c('green', 'orange'),
                  nmols = 1960,
                  mols.size = 0.3,
                  coord.fixed = TRUE)


ggplot2::ggsave(
  file.path(output_dir, "CM1799-A1-Th1Fp1U1_MET_HGF_spatial_plot.png"),
  plot = p5, width = 6, height = 6, units = "in", dpi = 300
)

# HBEGF count 26266 on Xenium Explorer
# EGFR count 34972 on Xenium Explorer

p8 = ImageDimPlot(CM1799C1_A1, 
                  fov = "fov",  
                  group.by = 'GTR_cat2',
                  cols = c("Tumor" = "violetred", "PSIZ" = "palegreen4", "Necrosis" = "#604e97", 'APT' = KellyPalette$grey),
                  molecules = c('HBEGF', 'EGFR'),
                  mols.col = c('gold', 'green'),
                  nmols = 2626,
                  mols.size = 0.3,
                  coord.fixed = TRUE)

ggplot2::ggsave(
  file.path(output_dir, "CM1799-A1-Th1Fp1U1_HBEGF_EGFR_spatial_plot.png"),
  plot = p8, width = 6, height = 6, units = "in", dpi = 300
)