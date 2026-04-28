# Load necessary libraries
library(Seurat)
library(infercnv)
library(tidyverse)
library(purrr)

# Set options
options(future.globals.maxSize = 1e9)
options("Seurat.object.assay.version" = "v3")
options(scipen = 100)

# Source additional functions if needed
#source("/diskmnt/Projects/Users/Evan.p/scripts/Rscript/jupyter_support_functions.R")

# Define paths
rds_path <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/integrated/RNA/'
output_base_dir <- '/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/20260420_inferCNV_results/'
dir.create(output_base_dir)

# Set working directory
setwd(output_base_dir)

# Load the Seurat object
message("Loading object from: ", paste0(rds_path, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean5.rds'))
mCRC <- readRDS(file = paste0(rds_path, '57_Integrated_normalized_mCRC_snRNA_noDB_v7_clean5.rds'))
print(mCRC)

# Subset to exclude normal tissue samples
message("Making subsets with tumor samples")
mCRC <- subset(mCRC, Tissue_Type != 'normal')

cell_type_file = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/snRNA_objects/mCRC_57_samples_clean3_metadata_cell_type_all_20260420.csv'
cell_type = read.table(cell_type_file, header = TRUE, sep=',', row.name=1)
head(cell_type, 3)

mCRC <- AddMetaData(mCRC, cell_type, col.name = 'cell_type_all4')

# Get unique sample identifiers
samples <- unique(mCRC$orig.ident)
message("Samples with tumors: ", paste(samples, collapse = ", "))

# Resume settings:
# - set start_after_sample to the last fully attempted sample in the previous run
# - keep NULL to process all samples
start_after_sample <- "CM357C1-S1"
if (!is.null(start_after_sample) && start_after_sample %in% samples) {
  start_idx <- match(start_after_sample, samples)
  if (start_idx < length(samples)) {
    samples <- samples[(start_idx + 1):length(samples)]
  } else {
    samples <- character(0)
  }
}
message("Samples queued for this run: ", paste(samples, collapse = ", "))

# Load gene order file
message("Loading gene order file")
gene_order_file <- read.delim(
  file = "/diskmnt/Projects/MetNet_analysis_2/Colorectal/InferCNV/src/gencode_v32_gene_name_with_band.txt",
  header = FALSE,
  stringsAsFactors = FALSE,
  sep = "\t"
)
rownames(gene_order_file) <- gene_order_file$V1
gene_order_file$V1 <- NULL
head(gene_order_file)

# Define the function to run inferCNV for each sample
run_inferCNV_by_samples <- function(sample, main_object, gene_order_file, output_base_dir) {
  message("Processing sample: ", sample)

  sample_output_dir <- file.path(output_base_dir, sample)
  infercnv_obj_post_path <- file.path(sample_output_dir, paste0(sample, '_infercnv_obj_post_run.rds'))

  if (file.exists(infercnv_obj_post_path)) {
    message("Skipping ", sample, " because post-run object already exists: ", infercnv_obj_post_path)
    return(NULL)
  }
  
  message("making output dir: ", sample)
  # Create output directory specific to the sample
  dir.create(sample_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset the main object for the current sample
  sample_object <- subset(main_object, orig.ident == sample)
  
  # Update cell type annotations
  sample_object@meta.data <- sample_object@meta.data %>%
    mutate(cell_type_all4 = case_when(
      cell_type_all4 %in% c(
        'APCDD1_CRC', 'Canonical_CRC_Intestine_Proliferation', 'Canonical_CRC_Intestine', 'CMETS', 
        'Canonical_CRC_Stem', 'Canonical_CRC_Stem_Proliferation', 'Neuroendocrine-like tumor'
      ) ~ 'Tumor',
      TRUE ~ cell_type_all
    ))
  
  # Check the number of tumor cells
  tumor_cells <- sum(sample_object@meta.data$cell_type_all4 == 'Tumor')
  message("Number of tumor cells in ", sample, ": ", tumor_cells)
  if (tumor_cells < 10) {
    message("Skipping ", sample, " due to insufficient tumor cells (<10).")
    return(NULL)
  }

  # Identify cell types with only one cell and remove them
  cell_type_counts <- sample_object@meta.data %>%
    group_by(cell_type_all4) %>%
    tally()

  # Get cell types with more than one cell
  valid_cell_types <- cell_type_counts %>%
    filter(n > 10) %>%
    pull(cell_type_all4)
  message("Valid cell types: ", paste(valid_cell_types, collapse = ", "))  

  if (length(valid_cell_types) < 2) {
    message("Skipping ", sample, " because less than two cell types have more than one cell.")
    return(NULL)
  }

  # Subset the sample_object to include only valid cell types
  sample_object <- sample_object %>% subset(cell_type_all4 %in% valid_cell_types)
  
  # Prepare data for inferCNV
  sample_object$cell_names <- colnames(sample_object)
  cell_type_file <- data.frame(
    cell_ID = sample_object$cell_names,
    cell_type = sample_object$cell_type_all4
  )
  rownames(cell_type_file) <- cell_type_file$cell_ID
  cell_type_file$cell_ID <- NULL
  RNA_counts_matrix <- sample_object@assays$RNA@counts
  Idents(sample_object) <- "cell_type_all4"
  
  message("Creating inferCNV object: ", sample)
  # Create inferCNV object
  ref_groups <- setdiff(unique(Idents(sample_object)), "Tumor")
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = RNA_counts_matrix,
    annotations_file = cell_type_file,
    gene_order_file = gene_order_file,
    ref_group_names = ref_groups
  )
  
  # Save the inferCNV object
  infercnv_obj_path <- file.path(sample_output_dir, paste0(sample, '_infercnv_obj.rds'))
  saveRDS(infercnv_obj, infercnv_obj_path)
  
  message("Running inferCNV object: ", sample)
  output_dir <- file.path(sample_output_dir, "out")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Run inferCNV
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    analysis_mode = 'subclusters',
    cutoff = 0.1,
    out_dir = output_dir,
    cluster_by_groups = TRUE,
    leiden_resolution = 0.0001,
    cluster_references = TRUE,
    plot_steps = FALSE,
    denoise = TRUE,
    HMM = TRUE,
    mask_nonDE_genes = TRUE,
    resume_mode = TRUE,
    num_threads = 20
  )
  
  # Save the post-run inferCNV object
  saveRDS(infercnv_obj, infercnv_obj_post_path)
}

# Run the function for each sample
walk(samples, run_inferCNV_by_samples, main_object = mCRC, gene_order_file = gene_order_file, output_base_dir = output_base_dir)


