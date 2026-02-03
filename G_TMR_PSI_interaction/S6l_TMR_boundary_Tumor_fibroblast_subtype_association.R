# TMR Boundary Tumor-Fibroblast Subtype Association Analysis (v2)
# This script analyzes correlations between tumor and fibroblast subtypes
# grouped by TMR boundary regions (tmr_group)
# Version 2: Filters with tumor_count > 100 & fibroblast_count > 100

# Load required libraries
library(tidyverse, quietly = TRUE)
library(ComplexHeatmap, quietly = TRUE)
library(circlize, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(gridExtra, quietly = TRUE)

# Load data
xenium_anno = '/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/mCRC_N26_Xenium_banky_celltype_metadata_20250713.csv'
xenium_metadata = read.csv(xenium_anno, header = TRUE)
rownames(xenium_metadata) <- xenium_metadata$barcode

tn_out_filter = read.csv('/diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_annotations/tn_out_filter.csv', 
                         header = TRUE, 
                         row.names = 1) |>
                         rename(tn_out_group = variable, tn_out_distance = value) |>
                         select(barcode, tn_out_group, tn_out_distance)
head(tn_out_filter, 3)

xenium_metadata <- left_join(xenium_metadata, tn_out_filter, by="barcode")
head(xenium_metadata, 3)

# Check column names
colnames(xenium_metadata)

# Create tumor-fibroblast metadata with TMR groups
tumor_fib_metadata = xenium_metadata |> 
                     filter(Broad_cell_type1 %in% c('Tumor', 'Fibroblast') & (tn_inward < 4 | tn_outward < 4)) |>
                     mutate(tn_group = paste0(Tissue_ID, '_tn', tn_label),
                            tn_out_group = if_else(is.na(tn_out_group), 
                                                   NA_character_, 
                                                   paste0(Tissue_ID, '_', tn_out_group)),
                            tmr_group = case_when(!is.na(tn_out_group) ~ tn_out_group,
                                                  !is.na(tn_group) ~ tn_group, 
                                                  TRUE ~ NA))
dim(tumor_fib_metadata)   
length(unique(tumor_fib_metadata$tmr_group))

# Filter out NA tmr_group
tumor_fib_metadata <- tumor_fib_metadata %>%
  filter(!is.na(tmr_group))
length(unique(tumor_fib_metadata$tmr_group))

# ============================================================================
# Colon Analysis
# ============================================================================

colon_tumor_fib_count <- tumor_fib_metadata %>%
  filter(Organ == 'Colon') %>%
  group_by(tmr_group) %>%
  summarise(
    all_cell_count = n(),
    tumor_count = sum(Broad_cell_type1 == "Tumor", na.rm = TRUE),
    fibroblast_count = sum(Broad_cell_type1 == "Fibroblast", na.rm = TRUE),
    proinvasive_tumor_count = sum(All_cell_type1 == "Non-canonical", na.rm = TRUE),
    stem_tumor_count = sum(All_cell_type1 == "Stem-like", na.rm = TRUE),
    prolif_tumor_count = sum(All_cell_type1 == "Proliferative-like", na.rm = TRUE),
    intestinal_tumor_count = sum(All_cell_type1 == "Intestine-like", na.rm = TRUE),
    mCAF_count = sum(All_cell_type1 == "mCAF", na.rm = TRUE),
    WNT5A_BMP_count = sum(All_cell_type1 == "WNT5A_BMP", na.rm = TRUE),
    WNT5A_infl_count = sum(All_cell_type1 == "WNT5A_infl", na.rm = TRUE),
    iCAF_count = sum(All_cell_type1 == "iCAF", na.rm = TRUE), 
    Organ = first(Organ),
    Tissue_ID = first(Tissue_ID),
    TMR_region = first(tmr_region),
    .groups = "drop"
  )

colon_tumor_fib_prop <- colon_tumor_fib_count %>%
    mutate(
        tumor_prop = round(tumor_count / all_cell_count, 3) ,
        fibroblast_prop = round(fibroblast_count / all_cell_count, 3),
        proinvasive_tumor_prop = round(proinvasive_tumor_count / tumor_count, 3),
        stem_tumor_prop = round(stem_tumor_count / tumor_count, 3),
        prolif_tumor_prop = round(prolif_tumor_count / tumor_count, 3),
        intestinal_tumor_prop = round(intestinal_tumor_count / tumor_count, 3),
        mCAF_prop = round(mCAF_count / fibroblast_count, 3),
        WNT5A_BMP_prop = round(WNT5A_BMP_count / fibroblast_count, 3),
        WNT5A_infl_prop = round(WNT5A_infl_count / fibroblast_count, 3),
        iCAF_prop = round(iCAF_count / fibroblast_count, 3),
        .groups = "drop"
    )

colon_tumor_fib_prop_filtered <- colon_tumor_fib_prop %>%
    filter(tumor_count > 100 & fibroblast_count > 100)
summary(colon_tumor_fib_prop_filtered$all_cell_count)
summary(colon_tumor_fib_prop_filtered$tumor_count)
summary(colon_tumor_fib_prop_filtered$fibroblast_count)

# Define tumor and fibroblast subtypes
# Order fibroblasts: mCAF, WNT5A_infl, WNT5A_BMP, iCAF
tumor_subtypes <- c("proinvasive_tumor_prop", "stem_tumor_prop", "prolif_tumor_prop", "intestinal_tumor_prop")
fibroblast_subtypes <- c("mCAF_prop", "WNT5A_infl_prop", "WNT5A_BMP_prop", "iCAF_prop")

cat("Tumor subtypes:", paste(tumor_subtypes, collapse=", "), "\n")
cat("Fibroblast subtypes:", paste(fibroblast_subtypes, collapse=", "), "\n")
cat("Total combinations:", length(tumor_subtypes) * length(fibroblast_subtypes), "\n")

# Create all combinations
combinations <- expand.grid(
  Tumor_subtype = tumor_subtypes,
  Fibroblast_subtype = fibroblast_subtypes,
  stringsAsFactors = FALSE
)

# Function to create a scatter plot for one combination
create_scatter_subtype <- function(tumor_st, fib_st, data) {
  # Calculate correlation and p-value
  cor_test <- cor.test(data[[tumor_st]], data[[fib_st]], use = "complete.obs", method = "spearman")
  cor_val <- cor_test$estimate
  p_val <- cor_test$p.value
  
  # Format p-value
  if (p_val < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- paste0("p = ", round(p_val, 3))
  }
  
  # Clean up subtype names for display
  tumor_name <- gsub("_prop", "", tumor_st)
  tumor_name <- gsub("_", " ", tumor_name)
  fib_name <- gsub("_prop", "", fib_st)
  fib_name <- gsub("_", " ", fib_name)
  
  p <- ggplot(data, aes_string(x = fib_st, y = tumor_st)) +
    geom_point(size = 2.5, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed", linewidth = 0.5) +
    labs(
      x = paste0("Fibroblast: ", fib_name),
      y = paste0("Tumor: ", tumor_name),
      title = paste0("Tumor: ", tumor_name, " vs Fibroblast: ", fib_name)
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

# Create plots for all combinations (Colon)
plots_results <- lapply(1:nrow(combinations), function(i) {
  create_scatter_subtype(combinations$Tumor_subtype[i], 
                        combinations$Fibroblast_subtype[i], 
                        colon_tumor_fib_prop_filtered)
})

# Extract plots
plots_list <- lapply(plots_results, function(x) x$plot)

# Arrange plots in a grid
n_cols <- length(fibroblast_subtypes)
n_rows <- length(tumor_subtypes)

# Create a grid of plots
grid_plot <- do.call(grid.arrange, c(plots_list, ncol = n_cols))

print(grid_plot)

# Save the plots
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Colon")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save the raw data table used for scatter plots
raw_data_filename <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_raw_data_Colon_v2.csv")
tryCatch({
  write.csv(colon_tumor_fib_prop_filtered, raw_data_filename, row.names = FALSE)
  cat("Raw data table saved to:", raw_data_filename, "\n")
}, error = function(e) {
  cat("Error saving raw data table:", e$message, "\n")
})

# Save as PDF (5:4 aspect ratio)
plot_filename_pdf <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_boundary_v2.pdf")
pdf(plot_filename_pdf, width = n_cols * 5, height = n_rows * 4)
do.call(grid.arrange, c(plots_list, ncol = n_cols))
dev.off()
cat("PDF plot saved to:", plot_filename_pdf, "\n")

# Create correlation and p-value matrices (Colon)
cat("\n=== Correlation Matrix (Colon) ===\n")
cor_matrix <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_matrix) <- tumor_subtypes
colnames(cor_matrix) <- fibroblast_subtypes

pval_matrix <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(pval_matrix) <- tumor_subtypes
colnames(pval_matrix) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
    cor_test <- cor.test(colon_tumor_fib_prop_filtered[[tumor_subtypes[i]]], 
                        colon_tumor_fib_prop_filtered[[fibroblast_subtypes[j]]], 
                        use = "complete.obs", method = "spearman")
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
cor_sig_matrix <- matrix("", nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_sig_matrix) <- tumor_subtypes
colnames(cor_sig_matrix) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
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

# Clean up names for display
clean_tumor_names <- c(
  "proinvasive_tumor_prop" = "Proinvasive",
  "stem_tumor_prop" = "Stem",
  "prolif_tumor_prop" = "Prolif",
  "intestinal_tumor_prop" = "Intestinal"
)

clean_fib_names <- c(
  "mCAF_prop" = "mCAF",
  "WNT5A_infl_prop" = "WNT5A infl.",
  "WNT5A_BMP_prop" = "WNT5A BMP",
  "iCAF_prop" = "iCAF"
)

# Create annotation matrix for significance markers
# Reorder pval_matrix to match the desired column order before creating sig_matrix
pval_matrix_ordered <- pval_matrix[, fibroblast_subtypes, drop = FALSE]
sig_matrix <- matrix("", nrow = nrow(pval_matrix_ordered), ncol = ncol(pval_matrix_ordered))
rownames(sig_matrix) <- clean_tumor_names[rownames(pval_matrix_ordered)]
colnames(sig_matrix) <- clean_fib_names[colnames(pval_matrix_ordered)]

for(i in 1:nrow(pval_matrix_ordered)) {
  for(j in 1:ncol(pval_matrix_ordered)) {
    p_val <- pval_matrix_ordered[i, j]
    if (p_val < 0.01) {
      sig_matrix[i, j] <- "**"
    } else if (p_val < 0.05) {
      sig_matrix[i, j] <- "*"
    }
  }
}

# Set row and column names for correlation matrix
# Explicitly reorder columns to: mCAF, WNT5A_infl, WNT5A_BMP, iCAF (before renaming)
cor_matrix_display <- cor_matrix[, fibroblast_subtypes, drop = FALSE]
rownames(cor_matrix_display) <- clean_tumor_names[rownames(cor_matrix_display)]
colnames(cor_matrix_display) <- clean_fib_names[colnames(cor_matrix_display)]

# Create color mapping (blue for negative, orange/red for positive)
col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#2166ac", "#67a9cf", "white", "#f4a582", "#ca0020")
)

# Create the heatmap
ht <- Heatmap(
  cor_matrix_display,
  name = "Correlation",
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add significance markers
    if(sig_matrix[i, j] != "") {
      grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(
    title = "Correlation",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

# Draw the heatmap with title
draw(ht, heatmap_legend_side = "bottom", 
     column_title = "Colon (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))

# Save the heatmap
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Colon")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save as PDF
heatmap_filename_pdf <- file.path(output_plot_dir, "Tumor_Fibroblast_correlation_heatmap_boundary_v2.pdf")
pdf(heatmap_filename_pdf, width = 4, height = 5)
draw(ht, heatmap_legend_side = "bottom", 
     column_title = "Colon (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
cat("PDF heatmap saved to:", heatmap_filename_pdf, "\n")

# ============================================================================
# Liver Analysis
# ============================================================================

liver_tumor_fib_count <- tumor_fib_metadata %>%
  filter(Organ == 'Liver') %>%
  group_by(tmr_group) %>%
  summarise(
    all_cell_count = n(),
    tumor_count = sum(Broad_cell_type1 == "Tumor", na.rm = TRUE),
    fibroblast_count = sum(Broad_cell_type1 == "Fibroblast", na.rm = TRUE),
    proinvasive_tumor_count = sum(All_cell_type1 == "Non-canonical", na.rm = TRUE),
    stem_tumor_count = sum(All_cell_type1 == "Stem-like", na.rm = TRUE),
    prolif_tumor_count = sum(All_cell_type1 == "Proliferative-like", na.rm = TRUE),
    intestinal_tumor_count = sum(All_cell_type1 == "Intestine-like", na.rm = TRUE),
    mCAF_count = sum(All_cell_type1 == "mCAF", na.rm = TRUE),
    WNT5A_BMP_count = sum(All_cell_type1 == "WNT5A_BMP", na.rm = TRUE),
    WNT5A_infl_count = sum(All_cell_type1 == "WNT5A_infl", na.rm = TRUE),
    iCAF_count = sum(All_cell_type1 == "iCAF", na.rm = TRUE), 
    Organ = first(Organ),
    Tissue_ID = first(Tissue_ID),
    TMR_region = first(tmr_region),
    .groups = "drop"
  )

liver_tumor_fib_prop <- liver_tumor_fib_count %>%
    mutate(
        tumor_prop = round(tumor_count / all_cell_count, 3) ,
        fibroblast_prop = round(fibroblast_count / all_cell_count, 3),
        proinvasive_tumor_prop = round(proinvasive_tumor_count / tumor_count, 3),
        stem_tumor_prop = round(stem_tumor_count / tumor_count, 3),
        prolif_tumor_prop = round(prolif_tumor_count / tumor_count, 3),
        intestinal_tumor_prop = round(intestinal_tumor_count / tumor_count, 3),
        mCAF_prop = round(mCAF_count / fibroblast_count, 3),
        WNT5A_BMP_prop = round(WNT5A_BMP_count / fibroblast_count, 3),
        WNT5A_infl_prop = round(WNT5A_infl_count / fibroblast_count, 3),
        iCAF_prop = round(iCAF_count / fibroblast_count, 3),
        .groups = "drop"
    )

liver_tumor_fib_prop_filtered <- liver_tumor_fib_prop %>%
    filter(tumor_count > 100 & fibroblast_count > 100)
summary(liver_tumor_fib_prop_filtered$all_cell_count)
summary(liver_tumor_fib_prop_filtered$tumor_count)
summary(liver_tumor_fib_prop_filtered$fibroblast_count)

# Create plots for all combinations (Liver)
plots_results_liver <- lapply(1:nrow(combinations), function(i) {
  create_scatter_subtype(combinations$Tumor_subtype[i], 
                        combinations$Fibroblast_subtype[i], 
                        liver_tumor_fib_prop_filtered)
})

# Extract plots
plots_list_liver <- lapply(plots_results_liver, function(x) x$plot)

# Arrange plots in a grid
n_cols <- length(fibroblast_subtypes)
n_rows <- length(tumor_subtypes)

# Create a grid of plots
grid_plot_liver <- do.call(grid.arrange, c(plots_list_liver, ncol = n_cols))

print(grid_plot_liver)

# Save the plots
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Liver")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save the raw data table used for scatter plots
raw_data_filename <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_raw_data_Liver_v2.csv")
tryCatch({
  write.csv(liver_tumor_fib_prop_filtered, raw_data_filename, row.names = FALSE)
  cat("Raw data table saved to:", raw_data_filename, "\n")
}, error = function(e) {
  cat("Error saving raw data table:", e$message, "\n")
})

# Save as PDF (5:4 aspect ratio)
plot_filename_pdf <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_boundary_v2.pdf")
pdf(plot_filename_pdf, width = n_cols * 5, height = n_rows * 4)
do.call(grid.arrange, c(plots_list_liver, ncol = n_cols))
dev.off()
cat("PDF plot saved to:", plot_filename_pdf, "\n")

# Create correlation and p-value matrices (Liver)
cat("\n=== Correlation Matrix (Liver) ===\n")
cor_matrix_liver <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_matrix_liver) <- tumor_subtypes
colnames(cor_matrix_liver) <- fibroblast_subtypes

pval_matrix_liver <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(pval_matrix_liver) <- tumor_subtypes
colnames(pval_matrix_liver) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
    cor_test <- cor.test(liver_tumor_fib_prop_filtered[[tumor_subtypes[i]]], 
                        liver_tumor_fib_prop_filtered[[fibroblast_subtypes[j]]], 
                        use = "complete.obs", method = "spearman")
    cor_matrix_liver[i, j] <- cor_test$estimate
    pval_matrix_liver[i, j] <- cor_test$p.value
  }
}

cat("\nCorrelation coefficients:\n")
print(round(cor_matrix_liver, 3))

cat("\nP-values:\n")
print(round(pval_matrix_liver, 4))

# Create a combined matrix with significance stars
cat("\n=== Correlation Matrix with Significance ===\n")
cor_sig_matrix_liver <- matrix("", nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_sig_matrix_liver) <- tumor_subtypes
colnames(cor_sig_matrix_liver) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
    cor_val <- cor_matrix_liver[i, j]
    p_val <- pval_matrix_liver[i, j]
    
    sig <- ""
    if (p_val < 0.001) {
      sig <- "***"
    } else if (p_val < 0.01) {
      sig <- "**"
    } else if (p_val < 0.05) {
      sig <- "*"
    }
    
    cor_sig_matrix_liver[i, j] <- paste0(round(cor_val, 3), sig)
  }
}

print(cor_sig_matrix_liver)
cat("\nSignificance: *** p<0.001, ** p<0.01, * p<0.05\n")

# Create annotation matrix for significance markers
# Reorder pval_matrix_liver to match the desired column order before creating sig_matrix_liver
pval_matrix_liver_ordered <- pval_matrix_liver[, fibroblast_subtypes, drop = FALSE]
sig_matrix_liver <- matrix("", nrow = nrow(pval_matrix_liver_ordered), ncol = ncol(pval_matrix_liver_ordered))
rownames(sig_matrix_liver) <- clean_tumor_names[rownames(pval_matrix_liver_ordered)]
colnames(sig_matrix_liver) <- clean_fib_names[colnames(pval_matrix_liver_ordered)]

for(i in 1:nrow(pval_matrix_liver_ordered)) {
  for(j in 1:ncol(pval_matrix_liver_ordered)) {
    p_val <- pval_matrix_liver_ordered[i, j]
    if (p_val < 0.01) {
      sig_matrix_liver[i, j] <- "**"
    } else if (p_val < 0.05) {
      sig_matrix_liver[i, j] <- "*"
    }
  }
}

# Set row and column names for correlation matrix
# Explicitly reorder columns to: mCAF, WNT5A_infl, WNT5A_BMP, iCAF (before renaming)
cor_matrix_display_liver <- cor_matrix_liver[, fibroblast_subtypes, drop = FALSE]
rownames(cor_matrix_display_liver) <- clean_tumor_names[rownames(cor_matrix_display_liver)]
colnames(cor_matrix_display_liver) <- clean_fib_names[colnames(cor_matrix_display_liver)]

# Create color mapping (blue for negative, orange/red for positive)
col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#2166ac", "#67a9cf", "white", "#f4a582", "#ca0020")
)

# Create the heatmap
ht_liver <- Heatmap(
  cor_matrix_display_liver,
  name = "Correlation",
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add significance markers
    if(sig_matrix_liver[i, j] != "") {
      grid.text(sig_matrix_liver[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(
    title = "Correlation",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

# Draw the heatmap with title
draw(ht_liver, heatmap_legend_side = "bottom", 
     column_title = "Liver (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))

# Save the heatmap
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Liver")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save as PDF
heatmap_filename_pdf <- file.path(output_plot_dir, "Tumor_Fibroblast_correlation_heatmap_boundary_v2.pdf")
pdf(heatmap_filename_pdf, width = 4, height = 5)
draw(ht_liver, heatmap_legend_side = "bottom", 
     column_title = "Liver (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
cat("PDF heatmap saved to:", heatmap_filename_pdf, "\n")

# ============================================================================
# Lung Analysis
# ============================================================================

lung_tumor_fib_count <- tumor_fib_metadata %>%
  filter(Organ == 'Lung') %>%
  group_by(tmr_group) %>%
  summarise(
    all_cell_count = n(),
    tumor_count = sum(Broad_cell_type1 == "Tumor", na.rm = TRUE),
    fibroblast_count = sum(Broad_cell_type1 == "Fibroblast", na.rm = TRUE),
    proinvasive_tumor_count = sum(All_cell_type1 == "Non-canonical", na.rm = TRUE),
    stem_tumor_count = sum(All_cell_type1 == "Stem-like", na.rm = TRUE),
    prolif_tumor_count = sum(All_cell_type1 == "Proliferative-like", na.rm = TRUE),
    intestinal_tumor_count = sum(All_cell_type1 == "Intestine-like", na.rm = TRUE),
    mCAF_count = sum(All_cell_type1 == "mCAF", na.rm = TRUE),
    WNT5A_BMP_count = sum(All_cell_type1 == "WNT5A_BMP", na.rm = TRUE),
    WNT5A_infl_count = sum(All_cell_type1 == "WNT5A_infl", na.rm = TRUE),
    iCAF_count = sum(All_cell_type1 == "iCAF", na.rm = TRUE), 
    Organ = first(Organ),
    Tissue_ID = first(Tissue_ID),
    TMR_region = first(tmr_region),
    .groups = "drop"
  )

lung_tumor_fib_prop <- lung_tumor_fib_count %>%
    mutate(
        tumor_prop = round(tumor_count / all_cell_count, 3) ,
        fibroblast_prop = round(fibroblast_count / all_cell_count, 3),
        proinvasive_tumor_prop = round(proinvasive_tumor_count / tumor_count, 3),
        stem_tumor_prop = round(stem_tumor_count / tumor_count, 3),
        prolif_tumor_prop = round(prolif_tumor_count / tumor_count, 3),
        intestinal_tumor_prop = round(intestinal_tumor_count / tumor_count, 3),
        mCAF_prop = round(mCAF_count / fibroblast_count, 3),
        WNT5A_BMP_prop = round(WNT5A_BMP_count / fibroblast_count, 3),
        WNT5A_infl_prop = round(WNT5A_infl_count / fibroblast_count, 3),
        iCAF_prop = round(iCAF_count / fibroblast_count, 3),
        .groups = "drop"
    )

lung_tumor_fib_prop_filtered <- lung_tumor_fib_prop %>%
    filter(tumor_count > 100 & fibroblast_count > 100)
summary(lung_tumor_fib_prop_filtered$all_cell_count)
summary(lung_tumor_fib_prop_filtered$tumor_count)
summary(lung_tumor_fib_prop_filtered$fibroblast_count)

# Create plots for all combinations (Lung)
plots_results_lung <- lapply(1:nrow(combinations), function(i) {
  create_scatter_subtype(combinations$Tumor_subtype[i], 
                        combinations$Fibroblast_subtype[i], 
                        lung_tumor_fib_prop_filtered)
})

# Extract plots
plots_list_lung <- lapply(plots_results_lung, function(x) x$plot)

# Arrange plots in a grid
n_cols <- length(fibroblast_subtypes)
n_rows <- length(tumor_subtypes)

# Create a grid of plots
grid_plot_lung <- do.call(grid.arrange, c(plots_list_lung, ncol = n_cols))

print(grid_plot_lung)

# Save the plots
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Lung")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save the raw data table used for scatter plots
raw_data_filename <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_raw_data_Lung_v2.csv")
tryCatch({
  write.csv(lung_tumor_fib_prop_filtered, raw_data_filename, row.names = FALSE)
  cat("Raw data table saved to:", raw_data_filename, "\n")
}, error = function(e) {
  cat("Error saving raw data table:", e$message, "\n")
})

# Save as PDF (5:4 aspect ratio)
plot_filename_pdf <- file.path(output_plot_dir, "Tumor_vs_Fibroblast_subtypes_scatter_plots_boundary_v2.pdf")
pdf(plot_filename_pdf, width = n_cols * 5, height = n_rows * 4)
do.call(grid.arrange, c(plots_list_lung, ncol = n_cols))
dev.off()
cat("PDF plot saved to:", plot_filename_pdf, "\n")

# Create correlation and p-value matrices (Lung)
cat("\n=== Correlation Matrix (Lung) ===\n")
cor_matrix_lung <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_matrix_lung) <- tumor_subtypes
colnames(cor_matrix_lung) <- fibroblast_subtypes

pval_matrix_lung <- matrix(NA, nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(pval_matrix_lung) <- tumor_subtypes
colnames(pval_matrix_lung) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
    cor_test <- cor.test(lung_tumor_fib_prop_filtered[[tumor_subtypes[i]]], 
                        lung_tumor_fib_prop_filtered[[fibroblast_subtypes[j]]], 
                        use = "complete.obs")
    cor_matrix_lung[i, j] <- cor_test$estimate
    pval_matrix_lung[i, j] <- cor_test$p.value
  }
}

cat("\nCorrelation coefficients:\n")
print(round(cor_matrix_lung, 3))

cat("\nP-values:\n")
print(round(pval_matrix_lung, 4))

# Create a combined matrix with significance stars
cat("\n=== Correlation Matrix with Significance ===\n")
cor_sig_matrix_lung <- matrix("", nrow = length(tumor_subtypes), ncol = length(fibroblast_subtypes))
rownames(cor_sig_matrix_lung) <- tumor_subtypes
colnames(cor_sig_matrix_lung) <- fibroblast_subtypes

for(i in 1:length(tumor_subtypes)) {
  for(j in 1:length(fibroblast_subtypes)) {
    cor_val <- cor_matrix_lung[i, j]
    p_val <- pval_matrix_lung[i, j]
    
    sig <- ""
    if (p_val < 0.001) {
      sig <- "***"
    } else if (p_val < 0.01) {
      sig <- "**"
    } else if (p_val < 0.05) {
      sig <- "*"
    }
    
    cor_sig_matrix_lung[i, j] <- paste0(round(cor_val, 3), sig)
  }
}

print(cor_sig_matrix_lung)
cat("\nSignificance: *** p<0.001, ** p<0.01, * p<0.05\n")

# Create annotation matrix for significance markers
# Reorder pval_matrix_lung to match the desired column order before creating sig_matrix_lung
pval_matrix_lung_ordered <- pval_matrix_lung[, fibroblast_subtypes, drop = FALSE]
sig_matrix_lung <- matrix("", nrow = nrow(pval_matrix_lung_ordered), ncol = ncol(pval_matrix_lung_ordered))
rownames(sig_matrix_lung) <- clean_tumor_names[rownames(pval_matrix_lung_ordered)]
colnames(sig_matrix_lung) <- clean_fib_names[colnames(pval_matrix_lung_ordered)]

for(i in 1:nrow(pval_matrix_lung_ordered)) {
  for(j in 1:ncol(pval_matrix_lung_ordered)) {
    p_val <- pval_matrix_lung_ordered[i, j]
    if (p_val < 0.01) {
      sig_matrix_lung[i, j] <- "**"
    } else if (p_val < 0.05) {
      sig_matrix_lung[i, j] <- "*"
    }
  }
}

# Set row and column names for correlation matrix
# Explicitly reorder columns to: mCAF, WNT5A_infl, WNT5A_BMP, iCAF (before renaming)
cor_matrix_display_lung <- cor_matrix_lung[, fibroblast_subtypes, drop = FALSE]
rownames(cor_matrix_display_lung) <- clean_tumor_names[rownames(cor_matrix_display_lung)]
colnames(cor_matrix_display_lung) <- clean_fib_names[colnames(cor_matrix_display_lung)]

# Create color mapping (blue for negative, orange/red for positive)
col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  c("#2166ac", "#67a9cf", "white", "#f4a582", "#ca0020")
)

# Create the heatmap
ht_lung <- Heatmap(
  cor_matrix_display_lung,
  name = "Correlation",
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Add significance markers
    if(sig_matrix_lung[i, j] != "") {
      grid.text(sig_matrix_lung[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11),
  heatmap_legend_param = list(
    title = "Correlation",
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(4, "cm"),
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1", "-0.5", "0", "0.5", "1")
  )
)

# Draw the heatmap with title
draw(ht_lung, heatmap_legend_side = "bottom", 
     column_title = "Lung (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))

# Save the heatmap
output_plot_dir <- "/diskmnt/Projects/MetNet_analysis/Colorectal/analysis/figure_4/supp"
output_plot_dir <- file.path(output_plot_dir, "TMR_tumor_fiborblast_association", "Lung")
dir.create(output_plot_dir, showWarnings = FALSE, recursive = TRUE)

# Save as PDF
heatmap_filename_pdf <- file.path(output_plot_dir, "Tumor_Fibroblast_correlation_heatmap_boundary_v2.pdf")
pdf(heatmap_filename_pdf, width = 4, height = 5)
draw(ht_lung, heatmap_legend_side = "bottom", 
     column_title = "Lung (TMR Boundary) v2", 
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
cat("PDF heatmap saved to:", heatmap_filename_pdf, "\n")

cat("\n=== Analysis Complete ===\n")
