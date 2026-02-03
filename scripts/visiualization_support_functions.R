Highlight_Cluster_UMAP <- function(seurat_object, metadata_column, cell_type, umap_reduction, highlight_color = "#B56727", background_color = "lightgray", background_alpha = 0.3) {
  Idents(seurat_object) <- metadata_column
  
  umap_data <- Embeddings(seurat_object, reduction = umap_reduction)
  
  umap_data <- as.data.frame(umap_data)
  umap_data$cluster <- Idents(seurat_object)
  
  umap_dims <- colnames(umap_data)[1:2]
  
  ggplot(umap_data, aes_string(x = umap_dims[1], y = umap_dims[2])) +
    geom_point(aes(color = ifelse(cluster == cell_type, "highlight", "background"), 
                   alpha = ifelse(cluster == cell_type, 1, background_alpha)), 
               size = 1) +
    scale_color_manual(values = c("highlight" = highlight_color, "background" = background_color)) +
    scale_alpha_continuous(range = c(background_alpha, 1)) +
    theme_minimal() +
    labs(title = paste("UMAP with Highlighted", cell_type, "Cluster")) +
    theme(legend.position = "none")
}

CustomDotPlot <- function(
  seurat_obj,
  features,
  cell_type_col,
  gene_annotation_table = NULL,
  gene_annotation_col = NULL,
  FDA_approved_list = NULL,  
  assay = "SCT",
  scale_min = -2.5,
  scale_max = 2.5,
  dot_min = 0,
  legend_pct_values = c(25, 50, 75, 100),
  max_radius_size = 5, # in mm
  color_palette = "RdYlBu",
  color_palette_length = 10
) {
  # Load required libraries
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
  
  # Define PercentAbove function
  PercentAbove <- function(x, threshold) {
    return(sum(x > threshold, na.rm = TRUE) / length(x))
  }
  
  # Set the DefaultAssay
  DefaultAssay(seurat_obj) <- assay
  
  # Set identities to your cell types
  Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_col]]
  
  # Define your cell types
  cell_types <- unique(Idents(seurat_obj))
  
  # Fetch data for the specified features
  data.features <- FetchData(object = seurat_obj, vars = features)
  
  # Add identities to the data
  data.features$id <- Idents(object = seurat_obj)
  
  # Initialize data frames to store results
  avg_exp_df <- data.frame(row.names = features)
  pct_exp_df <- data.frame(row.names = features)
  
  # Loop over each cell type
  for (cell_type in cell_types) {
    # Get cells in the current cell type
    cells_in_type <- WhichCells(object = seurat_obj, idents = cell_type)
    
    # Subset data for current cell type
    cluster_data <- data.features[cells_in_type, features, drop = FALSE]
    
    # Calculate average expression
    avg_exp <- apply(cluster_data, 2, function(x) mean(x, na.rm = TRUE))
    
    # Calculate percent expression
    pct_exp <- apply(cluster_data, 2, function(x) PercentAbove(x, threshold = 0)) * 100
    
    # Add results to data frames
    avg_exp_df <- cbind(avg_exp_df, avg_exp)
    pct_exp_df <- cbind(pct_exp_df, pct_exp)
  }
  
  # Set column names to cell types
  colnames(avg_exp_df) <- cell_types
  colnames(pct_exp_df) <- cell_types
  
  # Convert data frames to matrices
  avg_exp_mat <- as.matrix(avg_exp_df)
  pct_exp_mat <- as.matrix(pct_exp_df)
  
  # Transpose the matrices to have cell types as rows and features as columns
  avg_exp_mat <- t(avg_exp_mat)
  pct_exp_mat <- t(pct_exp_mat)
  
  # Adjust row and column names after transposition
  rownames(avg_exp_mat) <- cell_types
  colnames(avg_exp_mat) <- features
  rownames(pct_exp_mat) <- cell_types
  colnames(pct_exp_mat) <- features
  
  # Log-transform the data
  avg_exp_mat_log <- log1p(avg_exp_mat)
  avg_exp_mat_scaled <- scale(avg_exp_mat_log)
  avg_exp_mat_scaled <- pmax(pmin(avg_exp_mat_scaled, scale_max), scale_min)
  
  # Define color function based on average expression values
  rdylbu_colors <- rev(colorRampPalette(brewer.pal(color_palette_length, color_palette))(color_palette_length))
  color_breaks <- seq(scale_min, scale_max, length.out = length(rdylbu_colors))
  col_fun <- colorRamp2(color_breaks, rdylbu_colors)
  
  # Prepare gene annotation if provided
  if (!is.null(gene_annotation_table) && !is.null(gene_annotation_col)) {
    # Remove duplicates, keeping the first occurrence
    gene_annotation <- gene_annotation_table[!duplicated(gene_annotation_table$gene), ]
    
    # Get the order of genes (features) in the heatmap
    heatmap_genes <- colnames(avg_exp_mat_scaled)
    
    # Reorder the annotation data to match the heatmap
    gene_annotation <- gene_annotation[match(heatmap_genes, gene_annotation$gene), ]
    
    # Check for any mismatches or missing annotations
    if (any(is.na(gene_annotation[[gene_annotation_col]]))) {
      warning("Some genes in the heatmap do not have corresponding annotations.")
    }
    
    # Define colors for the annotations
    annotation_categories <- unique(gene_annotation[[gene_annotation_col]])
    
    # Use a color palette with enough colors
    num_categories <- length(annotation_categories)
    if (num_categories <= 8) {
      color_palette_annotation <- brewer.pal(n = num_categories, name = "Dark2")
    } else {
      color_palette_annotation <- colorRampPalette(brewer.pal(8, "Dark2"))(num_categories)
    }
    annotation_colors <- setNames(color_palette_annotation, annotation_categories)
    
    # Create the bottom annotation
    bottom_annotation <- HeatmapAnnotation(
      Annotation = gene_annotation[[gene_annotation_col]],
      col = list(Annotation = annotation_colors),
      annotation_label = gene_annotation_col,
      show_annotation_name = TRUE
    )
  } else {
    bottom_annotation <- NULL
  }

  # Prepare FDA approval annotation if provided
  if (!is.null(FDA_approved_list)) {
    FDA_approved <- FDA_approved_list[features]  # Match FDA status with features
    
    if (any(is.na(FDA_approved))) {
      warning("Some genes in the heatmap do not have corresponding FDA approval status.")
    }
    
    # Assign colors for FDA approval (TRUE = green, FALSE = red)
    FDA_colors <- c("TRUE" = "seagreen2", "FALSE" = "gray60")
    
    FDA_annotation <- HeatmapAnnotation(
      FDA_Approved = FDA_approved,
      col = list(FDA_Approved = FDA_colors),
      annotation_label = "FDA Approved",
      show_annotation_name = TRUE
    )
    
    # Combine annotations
    if (!is.null(bottom_annotation)) {
      bottom_annotation <- c(bottom_annotation, FDA_annotation)
    } else {
      bottom_annotation <- FDA_annotation
    }
  }
  
  # Define the maximum radius for dots in the heatmap
  max_radius <- unit(max_radius_size, "mm")  # Adjust as needed
  
  # Compute the corresponding radii for the legend dots
  legend_dots <- sqrt(legend_pct_values / 100) * unit(max_radius_size, "mm")
  
  # Create the legend object
  percent_expression_legend <- Legend(
    at = legend_pct_values,
    labels = paste0(legend_pct_values, "%"),
    type = "points",
    legend_gp = gpar(col = "black", fill = "black"),
    background = NA,
    size = legend_dots,
    grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm"),
    title = "Percent Expressed",
    direction = "horizontal",
    nrow = 1,
    labels_gp = gpar(fontsize = 10)
  )
  
  # Define the spacing for grid.rect
  spacing <- unit(0.5, "pt")  # Adjust as needed
  
  # Update your cell_fun to use consistent max_radius
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Draw grid rectangle (cell border)
    grid.rect(
      x = x, y = y, width = width - spacing, height = height - spacing,
      gp = gpar(col = "grey", fill = NA)
    )
    
    # Get scaled average and percent expression values
    avg_exp_value <- avg_exp_mat_scaled[i, j]
    pct_exp_value <- pct_exp_mat[i, j]
    
    # Apply minimum percent expression threshold if needed
    if (pct_exp_value < dot_min) {
      return()  # Skip drawing if below threshold
    }
    
    # Convert max_radius to "npc" units
    max_radius_npc <- convertUnit(max_radius, unitTo = "npc", valueOnly = TRUE)
    
    # Calculate the radius
    radius <- sqrt(pct_exp_value / 100) * max_radius_npc
    
    # Draw circle with color based on scaled average expression
    grid.circle(
      x = x, y = y, r = unit(radius, "npc"),
      gp = gpar(fill = col_fun(avg_exp_value), col = NA)
    )
  }
  
  # Create the heatmap object
  heatmap_obj <- Heatmap(
    matrix = avg_exp_mat_scaled,
    name = "Average Expression",
    col = col_fun,
    rect_gp = gpar(type = "none"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_side = "bottom",
    bottom_annotation = bottom_annotation,
    heatmap_legend_param = list(
      title = "Scaled Avg. Expression",
      at = c(scale_min, 0, scale_max),
      labels = c(scale_min, "0", scale_max)
    ),
    cell_fun = cell_fun
  )
  
  return(list(
    heatmap = heatmap_obj,
    percent_expression_legend = percent_expression_legend
  ))
}
