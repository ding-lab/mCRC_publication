library(Seurat)
library(tidyverse)
library(qs)
library(future)
options(future.globals.maxSize = 30 * 1024^3)
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Plot/ImagePlotFaster_v1.R')
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool_github/ImageBinPlot/R/ImageBinPlot.R')
source("/diskmnt/Users2/epeng/Projects/mCRC/scripts/jupyter_support_functions.R")

# Read input CSV file with sample_id and rds_file_path
# Modify this path to your CSV file
input_csv_path = './input_sample2.csv'  # Update this path to your CSV file
samples_df = read.csv(input_csv_path, header = TRUE)

# Verify required columns exist
if (!all(c('sample_id', 'rds_file_path') %in% colnames(samples_df))) {
    stop("CSV file must contain 'sample_id' and 'rds_file_path' columns")
}

# Load annotation metadata
xenium_annotation_dir = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/'
xenium_annotation_path = file.path(xenium_annotation_dir, 'mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv')
xenium_annotation = read.csv(xenium_annotation_path, header = TRUE) %>% 
                    mutate(GTR = if_else(tns_label > 0, 'GTR', 'APT')) %>% 
                    mutate(TMR = if_else(tn_label > 0, 'TMR', 'non-TMR')) %>% 
                    mutate(GTR_cat = case_when(
                        tn_label > 0 ~ 'TMR',
                        tns_label == 0 ~ 'APT',
                        TRUE ~ 'PSI'
                    ))

NB_colors = c('NB1_necrosis' = '#604e97', 'NB2_tumor_core' = 'hotpink4', 
              'NB3_tumor_body' = 'violetred', 'NB4_tumor_stroma_interface' = 'orange1',
              'NB5_desmoplastic_stroma' = 'darkgreen', 'NB6_perivascular_stroma' = 'yellowgreen',
              'NB7_myelocyte_enriched_stroma' = '#be0032', 'NB8_lymphocyte_enriched_stroma' = '#dcd300', 
              'NB9_plasma_cell_enriched_stroma' = '#0067a5', 'NB10_smooth_muscle_stroma' = 'lightseagreen',
              'NB11_colon_epithelium' = 'gold2', 'NB12_germinal_center' = 'yellow4', 'NB13_liver' = '#882d17', 
              'NB14_ductal_epithelium' = '#2b3d26',  'Unassined' = KellyPalette$grey)

gtr_col <- list(
    GTR_cat = c(
        'APT' = "grey80",  
        'PSI' = "palegreen4",
        'TMR' = "firebrick2"),
    neighborhoods = NB_colors
)

# Get script directory for output
script_dir = getwd()

# Process each sample
for (i in 1:nrow(samples_df)) {
    sample_id = samples_df$sample_id[i]
    rds_file_path = samples_df$rds_file_path[i]
    
    cat(sprintf("Processing sample: %s\n", sample_id))
    cat(sprintf("Reading RDS file: %s\n", rds_file_path))
    
    # Check if RDS file exists
    if (!file.exists(rds_file_path)) {
        warning(sprintf("RDS file not found: %s. Skipping sample %s.\n", rds_file_path, sample_id))
        next
    }
    
    # Read Seurat object
    if (endsWith(rds_file_path, '.qs')) {
        obj = qread(rds_file_path)
    } else {
        obj = readRDS(rds_file_path)
    }
    
    # Add annotation metadata
    sample_meta = xenium_annotation %>% filter(corrected_sample_id == !!sample_id)
    if (nrow(sample_meta) > 0) {
        rownames(sample_meta) <- sample_meta$cell_id
        obj <- AddMetaData(obj, sample_meta)
    } else {
        warning(sprintf("No annotation metadata found for sample %s. Skipping metadata addition.\n", sample_id))
    }
    
    # Generate spatial plot
    cat(sprintf("Generating spatial plot for sample: %s\n", sample_id))
    
    # Get FOV name - try to detect it or use default
    fov_names = GetFovNames(obj, ident_column = 'orig.ident', fov = 'fov')
    if (length(fov_names) > 0) {
        fov_use = fov_names[1]  # Use first FOV
    } else {
        fov_use = 'fov'  # Default
    }
    
    p1 = ImagePlotFaster(obj=obj, 
                    group.by='GTR_cat', 
                    fov=fov_use, 
                    groupby_palette = gtr_col,
                    pointsize = 3,     
                    flip_y = TRUE)

    p2 = ImagePlotFaster(obj=obj, 
                group.by='neighborhoods', 
                fov=fov_use, 
                groupby_palette = gtr_col,
                pointsize = 3,     
                flip_y = TRUE)
    
    # Save plot to script directory
    output_file = file.path(script_dir, sprintf('%s_GTR_cat.pdf', sample_id))
    cat(sprintf("Saving plot to: %s\n", output_file))
    pdf(file=output_file, width=8, height=8) 
    print(p1)
    dev.off()

    output_file = file.path(script_dir, sprintf('%s_neighborhoods.pdf', sample_id))
    cat(sprintf("Saving plot to: %s\n", output_file))
    pdf(file=output_file, width=8, height=8) 
    print(p2)
    dev.off()
    
    cat(sprintf("Completed sample: %s\n\n", sample_id))
}

cat("All samples processed!\n")