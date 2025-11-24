#!/usr/bin/env Rscript

# Visium Test - Tumor Boundary Analysis
# This script analyzes correlations between CRC and FIB genesets in tumor boundary regions

# Load required libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)

# ============================================================================
# Data Loading
# ============================================================================

# Load gene sets
geneset_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/genesets"
CRC_genesets <- readRDS(file.path(geneset_path, "WASHU_snRNA_CRC_epithelial_genesets.rds"))
FIB_genesets <- readRDS(file.path(geneset_path, "WASHU_snRNA_CRC_fibroblast_genesets.rds"))

# Check gene set names
names(CRC_genesets)
names(FIB_genesets)

# Load Visium tracking information
Visium_tracking_path <- "/diskmnt/Users2/epeng/Projects/mCRC/scripts/Visium/CRC_Visium_tracking.csv"
Visium_tracking <- read_csv(Visium_tracking_path, show_col_types = FALSE)
head(Visium_tracking)

# Check column names
colnames(Visium_tracking)

# Check annotation file structure (first sample)
head(read_csv(Visium_tracking[[1, "MorphOutput_FINAL"]], show_col_types = FALSE))

# ============================================================================
# Tumor Boundary Analysis
# ============================================================================

# Read pseudobulk boundary AUCell scores
pseudobulk_auc_path <- "/diskmnt/Projects/MetNet_analysis/Colorectal/data/Visium_annotations/AUCell_score/Pseudobulk/mCRC_Visium_pseudobulk_boundary_genesets_auc_score.csv"
pseudobulk_auc <- read_csv(pseudobulk_auc_path, show_col_types = FALSE)
head(pseudobulk_auc)
dim(pseudobulk_auc)
colnames(pseudobulk_auc)

# Process pseudobulk AUC data
# Extract sample name and boundary type from row names
pseudobulk_auc <- pseudobulk_auc %>%
  mutate(
    sample_boundary = ...1,  # Assuming first column is row names
    sample_name = str_extract(sample_boundary, "^[^_]+"),
    boundary_type = str_extract(sample_boundary, "(?<=_)(Tumor|TME)(?=_)")
  ) %>%
  select(-...1)

# Get case_id and Organ from tracking file
sample_case_map <- Visium_tracking %>%
  select(LibraryName, case_id, Organ) %>%
  distinct()

# Add case_id to pseudobulk data
pseudobulk_auc <- pseudobulk_auc %>%
  left_join(sample_case_map, by = c("sample_name" = "LibraryName"))

head(pseudobulk_auc)

# Identify CRC (tumor) and FIB geneset columns
CRC_gene_set_names <- names(CRC_genesets)
CRC_selected <- c('Canonical_CRC_Intestine', 
                  'Canonical_CRC_Intestine_Proliferation',
                  'Canonical_CRC_Stem',
                  'Canonical_CRC_Stem_Proliferation',
                  'Non_Canonical_CRC_1',
                  'Non_Canonical_CRC_2')
# Clean names to match column names
CRC_selected_clean <- gsub(" ", "_", CRC_selected)
CRC_selected_clean <- gsub("-", "_", CRC_selected_clean)

FIB_gene_set_names <- names(FIB_genesets)
FIB_clean <- gsub(" ", "_", FIB_gene_set_names)
FIB_clean <- gsub("-", "_", FIB_clean)

cat("CRC genesets:", paste(CRC_selected_clean, collapse=", "), "\n")
cat("FIB genesets:", paste(FIB_clean, collapse=", "), "\n")

# Extract tumor boundary scores (CRC genesets) and TME boundary scores (FIB genesets)
# Keep individual geneset scores (don't calculate mean)
tumor_boundary_data <- pseudobulk_auc %>%
  filter(boundary_type == "Tumor") %>%
  select(sample_name, case_id, all_of(CRC_selected_clean))

tme_boundary_data <- pseudobulk_auc %>%
  filter(boundary_type == "TME") %>%
  select(sample_name, case_id, all_of(FIB_clean))

head(tumor_boundary_data)
head(tme_boundary_data)

# Average each individual geneset score by case_id
tumor_by_case <- tumor_boundary_data %>%
  group_by(case_id) %>%
  summarise(
    across(all_of(CRC_selected_clean), ~ mean(.x, na.rm = TRUE)),
    n_samples = n(),
    .groups = 'drop'
  )

fib_by_case <- tme_boundary_data %>%
  group_by(case_id) %>%
  summarise(
    across(all_of(FIB_clean), ~ mean(.x, na.rm = TRUE)),
    n_samples = n(),
    .groups = 'drop'
  )

# Get Organ information for each case_id
case_organ_map <- Visium_tracking %>%
  select(case_id, Organ) %>%
  distinct()

# Combine data for plotting
plot_data <- tumor_by_case %>%
  inner_join(fib_by_case, by = "case_id") %>%
  left_join(case_organ_map, by = "case_id") %>%
  filter(!is.na(case_id))

head(plot_data)
cat("Number of cases:", nrow(plot_data), "\n")
cat("CRC genesets:", length(CRC_selected_clean), "\n")
cat("FIB genesets:", length(FIB_clean), "\n")
cat("Total combinations:", length(CRC_selected_clean) * length(FIB_clean), "\n")

# Create scatter plots for each combination of CRC × FIB genesets
# Create all combinations
combinations <- expand.grid(
  CRC_geneset = CRC_selected_clean,
  FIB_geneset = FIB_clean,
  stringsAsFactors = FALSE
)

# Function to create a scatter plot for one combination
create_scatter <- function(crc_gs, fib_gs, data) {
  # Calculate correlation and p-value using Spearman
  cor_test <- cor.test(data[[crc_gs]], data[[fib_gs]])
  cor_val <- cor_test$estimate
  p_val <- cor_test$p.value
  
  # Format p-value
  if (p_val < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- paste0("p = ", round(p_val, 3))
  }
  
  # Define color palette for Organ
  organ_colors <- c("Colon" = "#C2B280", "Liver" = "brown", "Lung" = "steelblue1")
  
  p <- ggplot(data, aes_string(x = fib_gs, y = crc_gs, color = "Organ")) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(name = "Organ", values = organ_colors) +
    labs(
      x = paste0("FIB: ", fib_gs),
      y = paste0("Tumor: ", crc_gs),
      title = paste0("Tumor: ", crc_gs, " vs FIB: ", fib_gs)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
  
  # Add correlation and p-value to plot
  p <- p + annotate("text", x = Inf, y = Inf, 
                    label = paste0("r = ", round(cor_val, 3), "\n", p_label),
                    hjust = 1.1, vjust = 1.5, size = 3)
  
  return(list(plot = p, cor = cor_val, pval = p_val))
}

# Create plots for all combinations
plots_results <- lapply(1:nrow(combinations), function(i) {
  create_scatter(combinations$CRC_geneset[i], 
                 combinations$FIB_geneset[i], 
                 plot_data)
})

# Extract plots
plots_list <- lapply(plots_results, function(x) x$plot)

# Arrange plots in a grid
# Calculate grid dimensions
n_cols <- length(FIB_clean)
n_rows <- length(CRC_selected_clean)

# Create a grid of plots
grid_plot <- do.call(grid.arrange, c(plots_list, ncol = n_cols))

print(grid_plot)

# Save the plot
output_plot_dir <- file.path(getwd(), 'Visium_output')
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save as PDF (5:4 aspect ratio)
plot_filename_pdf <- file.path(output_plot_dir, "FIB_vs_Tumor_genesets_scatter_plots.pdf")
pdf(plot_filename_pdf, width = n_cols * 5, height = n_rows * 4)
do.call(grid.arrange, c(plots_list, ncol = n_cols))
dev.off()
cat("PDF plot saved to:", plot_filename_pdf, "\n")

# Save as PNG (5:4 aspect ratio)
plot_filename_png <- file.path(output_plot_dir, "FIB_vs_Tumor_genesets_scatter_plots.png")
png(plot_filename_png, width = n_cols * 5 * 300, height = n_rows * 4 * 300, res = 300)
do.call(grid.arrange, c(plots_list, ncol = n_cols))
dev.off()
cat("PNG plot saved to:", plot_filename_png, "\n")

# Create correlation and p-value matrices
cat("\n=== Correlation Matrix ===\n")
cor_matrix <- matrix(NA, nrow = length(CRC_selected_clean), ncol = length(FIB_clean))
rownames(cor_matrix) <- CRC_selected_clean
colnames(cor_matrix) <- FIB_clean

pval_matrix <- matrix(NA, nrow = length(CRC_selected_clean), ncol = length(FIB_clean))
rownames(pval_matrix) <- CRC_selected_clean
colnames(pval_matrix) <- FIB_clean

for(i in 1:length(CRC_selected_clean)) {
  for(j in 1:length(FIB_clean)) {
    cor_test <- cor.test(plot_data[[CRC_selected_clean[i]]], 
                        plot_data[[FIB_clean[j]]], 
                        method = "spearman")
    cor_matrix[i, j] <- cor_test$estimate
    pval_matrix[i, j] <- cor_test$p.value
  }
}

cat("\nCorrelation coefficients:\n")
print(round(cor_matrix, 3))

cat("\nP-values:\n")
print(round(pval_matrix, 4))

# Create a combined matrix with significance stars
cat("\n=== Correlation Matrix with Significance ===\n")
cor_sig_matrix <- matrix("", nrow = length(CRC_selected_clean), ncol = length(FIB_clean))
rownames(cor_sig_matrix) <- CRC_selected_clean
colnames(cor_sig_matrix) <- FIB_clean

for(i in 1:length(CRC_selected_clean)) {
  for(j in 1:length(FIB_clean)) {
    cor_val <- cor_matrix[i, j]
    p_val <- pval_matrix[i, j]
    
    sig <- ""
    if (p_val < 0.001) {
      sig <- "***"
    } else if (p_val < 0.01) {
      sig <- "**"
    } else if (p_val < 0.05) {
      sig <- "*"
    }
    
    cor_sig_matrix[i, j] <- paste0(round(cor_val, 3), sig)
  }
}

print(cor_sig_matrix)
cat("\nSignificance: *** p<0.001, ** p<0.01, * p<0.05\n")

cat("\n=== Analysis Complete ===\n")

