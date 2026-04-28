#!/usr/bin/env Rscript

# ==============================================================================
# mCRC epithelial-cell annotation
# Clean reproducible script for final manuscript-level epithelial labels
#
# This script keeps only the steps that determine the final epithelial annotation.
# Exploratory marker checks, UMAP highlight loops,
# unused subcluster branches, clinical boxplots, and downstream proportion
# analyses are not included.
#
# Final annotation columns:
#   epi_cell_type   : first-pass epithelial/tumor annotation
#   epi_cell_type3  : canonical tumor-state refinement before final collapse
#   epi_cell_type6  : final manuscript-level annotation
# 
# Note: angiogenic tumor was finally named as CMETS 
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
})

set.seed(1)
options(stringsAsFactors = FALSE)

# ------------------------------- User settings --------------------------------

graph_name <- "RNA_snn"

paths <- list(
  # Initial objects
  tumor_reint_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/tumor/57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_reINT.rds",

  epithelial_reint_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_reINT.rds",

  # Re-integrated epithelial object after the first cleanup step.
  # This object is generated upstream from epithelial_clean1.rds.
  epithelial_clean1_reint_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1_reINT.rds",

  # Re-integrated epithelial object used for final label assignment.
  # In the original workflow this object contains epithelial_cluster_sub6.
  # Change this path if your current final pre-annotation object has another name.
  epithelial_final_reint_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5.rds",

  output_dir =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/Analysis/snRNA/01_Cohort/Metadata"
)

outputs <- list(
  tumor_clean1_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/tumor/57_Integrated_normalized_mCRC_snRNA_noDB_v7_tumor_clean1.rds",

  epithelial_clean1_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean1.rds",

  epithelial_clean3_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean3.rds",

  epithelial_clean4_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean4.rds",

  epithelial_final_rds =
    "/diskmnt/Projects/MetNet_analysis_2/Colorectal/snMultiome/rds_objects/cell_types_subset/RNA/epithelial/57_Integrated_normalized_mCRC_snRNA_noDB_v7_epithelial_clean5_final_annotated.rds",

  metadata_csv =
    file.path(paths$output_dir, "57_Integrated_normalized_mCRC_snRNA_noDB_v7_final_epithelial_cell_type.csv"),

  final_umap_pdf =
    file.path(paths$output_dir, "Dimplot_mCRC_final_epithelial_cell_type6.pdf"),

  final_dotplot_pdf =
    file.path(paths$output_dir, "Dotplot_mCRC_final_epithelial_cell_type6.pdf")
)

dir.create(paths$output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- Helper functions -----------------------------

check_cols <- function(obj, cols, object_name = deparse(substitute(obj))) {
  missing_cols <- setdiff(cols, colnames(obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(
      object_name, " is missing metadata columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

check_graph <- function(obj, graph_name, object_name = deparse(substitute(obj))) {
  if (!graph_name %in% names(obj@graphs)) {
    stop(
      object_name, " does not contain graph '", graph_name, "'. ",
      "Run FindNeighbors upstream or change graph_name.",
      call. = FALSE
    )
  }
}

keep_by_metadata <- function(obj, column, values) {
  check_cols(obj, column)
  cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[column]]) %in% as.character(values)]
  subset(obj, cells = cells)
}

drop_by_metadata <- function(obj, column, values) {
  check_cols(obj, column)
  cells <- rownames(obj@meta.data)[!as.character(obj@meta.data[[column]]) %in% as.character(values)]
  subset(obj, cells = cells)
}

plot_dim <- function(obj, group_by, reduction, out_pdf, width = 7, height = 6, colors = NULL) {
  p <- DimPlot(
    obj,
    reduction = reduction,
    group.by = group_by,
    label = FALSE,
    raster = FALSE
  )

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors, name = group_by)
  }

  pdf(out_pdf, width = width, height = height)
  print(p)
  dev.off()

  invisible(p)
}

# ==============================================================================
# 1. Clean tumor epithelial states and remove low-confidence tumor subclusters
# ==============================================================================

tumor <- readRDS(paths$tumor_reint_rds)
DefaultAssay(tumor) <- "SCT"
check_graph(tumor, graph_name, "tumor")

tumor <- FindClusters(
  tumor,
  resolution = 0.1,
  graph.name = graph_name,
  cluster.name = "tumor_clusters_0.1"
)

tumor <- keep_by_metadata(tumor, "tumor_clusters_0.1", 0:4)

Idents(tumor) <- "tumor_clusters_0.1"
tumor <- FindSubCluster(
  object = tumor,
  cluster = "4",
  graph.name = graph_name,
  subcluster.name = "tumor_clusters_sub4",
  resolution = 0.1,
  algorithm = 1
)

Idents(tumor) <- "tumor_clusters_sub4"
tumor <- FindSubCluster(
  object = tumor,
  cluster = "2",
  graph.name = graph_name,
  subcluster.name = "tumor_clusters_sub4_2",
  resolution = 0.1,
  algorithm = 1
)

Idents(tumor) <- "tumor_clusters_sub4_2"
tumor <- FindSubCluster(
  object = tumor,
  cluster = "0",
  graph.name = graph_name,
  subcluster.name = "tumor_clusters_sub4_2_0",
  resolution = 0.1,
  algorithm = 1
)

Idents(tumor) <- "tumor_clusters_sub4_2_0"
tumor <- FindSubCluster(
  object = tumor,
  cluster = "1",
  graph.name = graph_name,
  subcluster.name = "tumor_clusters_0.1_sub",
  resolution = 0.1,
  algorithm = 1
)

tumor@meta.data <- tumor@meta.data %>%
  mutate(
    tumor_cluster = case_when(
      tumor_clusters_0.1_sub %in% c("0_0", "0_2", "0_3", "0_4") ~ "0",
      tumor_clusters_0.1_sub == "0_1" ~ "1",
      tumor_clusters_0.1_sub %in% c("1_0", "1_2", "1_3") ~ "2",
      tumor_clusters_0.1_sub == "1_1" ~ "3",
      tumor_clusters_0.1_sub %in% c("2_0", "2_1") ~ "4",
      tumor_clusters_0.1_sub == "2_2" ~ "5",
      tumor_clusters_0.1_sub == "3" ~ "6",
      tumor_clusters_0.1_sub %in% c("4_0", "4_2") ~ "7",
      tumor_clusters_0.1_sub == "4_1" ~ "8",
      TRUE ~ as.character(tumor_clusters_0.1_sub)
    )
  )

# Clusters 7 and 8 were excluded from final annotation:
#   7: T-cell-like/doublet signal
#   8: high mitochondrial/dead-cell-like signal
#   Other removed clusters were only have one cell
tumor_clean <- tumor %>%
  drop_by_metadata("tumor_cluster", c("7", "8"))

saveRDS(tumor_clean, outputs$tumor_clean1_rds)

# ==============================================================================
# 2. First epithelial cleanup using tumor-derived labels
# ==============================================================================

epithelial <- readRDS(paths$epithelial_reint_rds)
DefaultAssay(epithelial) <- "SCT"
check_graph(epithelial, graph_name, "epithelial")
check_cols(epithelial, "cell_type_simple3", "epithelial")

epithelial <- FindClusters(
  epithelial,
  resolution = 0.1,
  graph.name = graph_name,
  cluster.name = "epithelial_clusters_0.1"
)

# Add tumor_cluster labels back onto the full epithelial object.
tumor_cluster_meta <- tumor@meta.data %>%
  select(tumor_cluster)

epithelial <- AddMetaData(
  object = epithelial,
  metadata = tumor_cluster_meta
)

epithelial@meta.data <- epithelial@meta.data %>%
  mutate(
    epithelial_cluster = if_else(
      !is.na(tumor_cluster),
      as.character(tumor_cluster),
      as.character(cell_type_simple3)
    )
  )

# Remove low-confidence tumor clusters (doublet like).
epithelial_clean1 <- epithelial %>%
  drop_by_metadata("epithelial_cluster", c("7", "8"))

epithelial_clean1 <- keep_by_metadata(
  epithelial_clean1,
  "epithelial_cluster",
  c(
    "0", "1", "2", "3", "4", "5", "6",
    "Enterocyte", "Goblet"
  )
)

saveRDS(epithelial_clean1, outputs$epithelial_clean1_rds)

# ==============================================================================
# 3. Second epithelial cleanup after re-integration of epithelial_clean1
# ==============================================================================

# This step uses the re-integrated clean1 object generated upstream.
epithelial_reint <- readRDS(paths$epithelial_clean1_reint_rds)
DefaultAssay(epithelial_reint) <- "SCT"
check_graph(epithelial_reint, graph_name, "epithelial_reint")

epithelial_reint <- FindClusters(
  epithelial_reint,
  resolution = 0.1,
  graph.name = graph_name,
  cluster.name = "epithelial_clusters_0.1"
)

# Remove the microglia-like/doublet cluster and retain epithelial-enriched clusters.
epithelial_clean2 <- keep_by_metadata(
  epithelial_reint,
  "epithelial_clusters_0.1",
  c("0", "1", "2", "3", "4", "9")
)

epithelial_clean2 <- FindClusters(
  epithelial_clean2,
  resolution = 0.1,
  graph.name = graph_name,
  cluster.name = "epithelial_clusters_0.1"
)

epithelial_clean2 <- keep_by_metadata(
  epithelial_clean2,
  "epithelial_clusters_0.1",
  c("0", "1", "2", "3", "4", "5")
)

Idents(epithelial_clean2) <- "epithelial_clusters_0.1"
epithelial_clean2 <- FindSubCluster(
  object = epithelial_clean2,
  cluster = "3",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub1",
  resolution = 0.1,
  algorithm = 1
)

Idents(epithelial_clean2) <- "epithelial_cluster_sub1"
epithelial_clean2 <- FindSubCluster(
  object = epithelial_clean2,
  cluster = "3_1",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub2",
  resolution = 0.3,
  algorithm = 1
)

epithelial_clean2@meta.data <- epithelial_clean2@meta.data %>%
  mutate(
    epithelial_cluster_sub2 = case_when(
      epithelial_cluster_sub2 %in% c("3_1_0", "3_1_2", "3_1_4") ~ "3_1",
      TRUE ~ as.character(epithelial_cluster_sub2)
    )
  )

Idents(epithelial_clean2) <- "epithelial_cluster_sub2"
epithelial_clean2 <- FindSubCluster(
  object = epithelial_clean2,
  cluster = "0",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub3",
  resolution = 0.1,
  algorithm = 1
)

Idents(epithelial_clean2) <- "epithelial_cluster_sub3"
epithelial_clean2 <- FindSubCluster(
  object = epithelial_clean2,
  cluster = "1",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub4",
  resolution = 0.1,
  algorithm = 1
)

epithelial_clean2 <- epithelial_clean2 %>%
  drop_by_metadata("epithelial_cluster_sub4", "1_3")

epithelial_clean2@meta.data <- epithelial_clean2@meta.data %>%
  mutate(
    epithelial_cluster_sub4 = case_when(
      epithelial_cluster_sub4 == "1_2" ~ "1_0",
      TRUE ~ as.character(epithelial_cluster_sub4)
    )
  )

# Remove low-confidence subcluster 1_1 and rerun UMAP for the cleaned object.
epithelial_clean3 <- epithelial_clean2 %>%
  drop_by_metadata("epithelial_cluster_sub4", "1_1")

epithelial_clean3 <- RunUMAP(
  epithelial_clean3,
  reduction = "integrated.scvi",
  dims = 1:30,
  reduction.name = "epithelial_umap2.scvi"
)

Idents(epithelial_clean3) <- "epithelial_cluster_sub4"
epithelial_clean3 <- FindSubCluster(
  object = epithelial_clean3,
  cluster = "2",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub5",
  resolution = 0.2,
  algorithm = 1
)

epithelial_clean3 <- epithelial_clean3 %>%
  drop_by_metadata("epithelial_cluster_sub5", "2_4")

Idents(epithelial_clean3) <- "epithelial_cluster_sub5"
epithelial_clean3 <- FindSubCluster(
  object = epithelial_clean3,
  cluster = "2_3",
  graph.name = graph_name,
  subcluster.name = "epithelial_cluster_sub6",
  resolution = 0.1,
  algorithm = 1
)

epithelial_clean3@meta.data <- epithelial_clean3@meta.data %>%
  mutate(
    epithelial_cluster_sub6 = case_when(
      epithelial_cluster_sub6 %in% c("2_2", "2_1") ~ "2_0",
      epithelial_cluster_sub6 %in% c("2_3_2", "2_3_0") ~ "2_3",
      TRUE ~ as.character(epithelial_cluster_sub6)
    )
  )

saveRDS(epithelial_clean3, outputs$epithelial_clean3_rds)

# ==============================================================================
# 4. First-pass epithelial/tumor cell-type annotation
# ==============================================================================

# Use the final re-integrated epithelial object used for manuscript labels.
# This object should contain epithelial_clusters_0.1 and/or epithelial_cluster_sub6
# after the upstream cleanup above.
epithelial_final <- readRDS(paths$epithelial_final_reint_rds)
DefaultAssay(epithelial_final) <- "SCT"
check_graph(epithelial_final, graph_name, "epithelial_final")

# If epithelial_cluster_sub3 is not already present, regenerate the final cluster
# labels from the final re-integrated epithelial object.
# Remove doublet-like clusters
if (!"epithelial_cluster_sub3" %in% colnames(epithelial_final@meta.data)) {
  check_cols(epithelial_final, "epithelial_clusters_0.1", "epithelial_final")

  epithelial_final <- keep_by_metadata(
    epithelial_final,
    "epithelial_clusters_0.1",
    c("0", "1", "10", "2", "3", "4")
  )

  epithelial_final <- FindClusters(
    epithelial_final,
    resolution = 0.1,
    graph.name = graph_name,
    cluster.name = "epi_clusters_0.1"
  )

  epithelial_final <- keep_by_metadata(
    epithelial_final,
    "epi_clusters_0.1",
    c("0", "1", "2", "3", "4")
  )

  Idents(epithelial_final) <- "epi_clusters_0.1"
  epithelial_final <- FindSubCluster(
    object = epithelial_final,
    cluster = "0",
    graph.name = graph_name,
    subcluster.name = "epithelial_cluster_sub0",
    resolution = 0.1,
    algorithm = 1
  )

  epithelial_final@meta.data <- epithelial_final@meta.data %>%
    mutate(
      epithelial_cluster_sub0 = case_when(
        epithelial_cluster_sub0 == "0_3" ~ "0_1",
        TRUE ~ as.character(epithelial_cluster_sub0)
      )
    )

  Idents(epithelial_final) <- "epithelial_cluster_sub0"
  epithelial_final <- FindSubCluster(
    object = epithelial_final,
    cluster = "2",
    graph.name = graph_name,
    subcluster.name = "epithelial_cluster_sub2",
    resolution = 0.1,
    algorithm = 1
  )

  Idents(epithelial_final) <- "epithelial_cluster_sub2"
  epithelial_final <- FindSubCluster(
    object = epithelial_final,
    cluster = "2_2",
    graph.name = graph_name,
    subcluster.name = "epithelial_cluster_sub2",
    resolution = 0.1,
    algorithm = 1
  )

  epithelial_final@meta.data <- epithelial_final@meta.data %>%
    mutate(
      epithelial_cluster_sub2 = case_when(
        epithelial_cluster_sub2 %in% c("2_1", "2_3") ~ "2_0",
        epithelial_cluster_sub2 == "2_2_2" ~ "2_2_0",
        TRUE ~ as.character(epithelial_cluster_sub2)
      )
    )

  Idents(epithelial_final) <- "epithelial_cluster_sub2"
  epithelial_final <- FindSubCluster(
    object = epithelial_final,
    cluster = "3",
    graph.name = graph_name,
    subcluster.name = "epithelial_cluster_sub3",
    resolution = 0.1,
    algorithm = 1
  )

  Idents(epithelial_final) <- "epithelial_cluster_sub3"
  epithelial_final <- FindSubCluster(
    object = epithelial_final,
    cluster = "3_1",
    graph.name = graph_name,
    subcluster.name = "epithelial_cluster_sub3",
    resolution = 0.2,
    algorithm = 1
  )

  epithelial_final@meta.data <- epithelial_final@meta.data %>%
    mutate(
      epithelial_cluster_sub3 = case_when(
        epithelial_cluster_sub3 %in% c("3_1_1", "3_1_3") ~ "3_1_0",
        TRUE ~ as.character(epithelial_cluster_sub3)
      )
    )
}

epithelial_final@meta.data <- epithelial_final@meta.data %>%
  mutate(
    epi_cell_type = case_when(
      epithelial_cluster_sub3 == "0_0" ~ "Stem-like tumor",
      epithelial_cluster_sub3 == "0_1" ~ "Angiogenic tumor",
      epithelial_cluster_sub3 == "0_2" ~ "Proliferative stem-like tumor",
      epithelial_cluster_sub3 == "1" ~ "Proliferative tumor",
      epithelial_cluster_sub3 == "2_0" ~ "APCDD1+ tumor",
      epithelial_cluster_sub3 == "2_2_0" ~ "Enteroendocrine-like cells",
      epithelial_cluster_sub3 == "2_2_1" ~ "Tuft cells",
      epithelial_cluster_sub3 == "3_0" ~ "Enterocytes",
      epithelial_cluster_sub3 == "3_1_0" ~ "Transit-amplifying cells",
      epithelial_cluster_sub3 == "3_1_2" ~ "Stem cells",
      epithelial_cluster_sub3 == "4" ~ "Goblet cells",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(epithelial_final$epi_cell_type))) {
  print(table(epithelial_final$epithelial_cluster_sub3, useNA = "ifany"))
  stop("Some epithelial_cluster_sub3 values were not mapped to epi_cell_type.", call. = FALSE)
}

epithelial_final$epi_cell_type <- factor(
  epithelial_final$epi_cell_type,
  levels = c(
    "APCDD1+ tumor", "Angiogenic tumor", "Proliferative tumor",
    "Proliferative stem-like tumor", "Stem-like tumor",
    "Stem cells", "Enterocytes", "Transit-amplifying cells",
    "Goblet cells", "Tuft cells", "Enteroendocrine-like cells"
  )
)

saveRDS(epithelial_final, outputs$epithelial_clean4_rds)

# ==============================================================================
# 5. Refine canonical tumor-state labels and generate final manuscript labels
# ==============================================================================

# Only the stem-like/canonical epithelial refinements are retained here.
# Exploratory subclustering of Angiogenic tumor and APCDD1+ tumor is omitted
# because those subclusters were not used for the final manuscript labels.
#
# Final non-canonical mapping:
#   Angiogenic tumor -> Non_Canonical_CRC_1 / CMETS
#   APCDD1+ tumor    -> Non_Canonical_CRC_2

# 5a. Split stem-like tumor into canonical intestine-like vs tumor-ISC-like states.
Idents(epithelial_final) <- "epi_cell_type"
epithelial_final <- FindSubCluster(
  object = epithelial_final,
  cluster = "Stem-like tumor",
  graph.name = graph_name,
  subcluster.name = "epi_cell_type_sub0",
  resolution = 0.3,
  algorithm = 1
)

epithelial_final@meta.data <- epithelial_final@meta.data %>%
  mutate(
    epi_cell_type2 = case_when(
      epi_cell_type_sub0 %in% c(
        "Stem-like tumor_0",
        "Stem-like tumor_4",
        "Stem-like tumor_5",
        "Stem-like tumor_6",
        "Stem-like tumor_7",
        "Stem-like tumor_8"
      ) ~ "Intestine-like",
      epi_cell_type_sub0 %in% c(
        "Stem-like tumor_1",
        "Stem-like tumor_2",
        "Stem-like tumor_3"
      ) ~ "Tumor-ISC-like"
      TRUE ~ as.character(epi_cell_type)
    )
  )

# 5b. Collapse retained labels into final manuscript-level epithelial labels.
# Angiogenic tumor and APCDD1+ tumor are mapped directly, without unused
# exploratory subcluster labels.
epithelial_final@meta.data <- epithelial_final@meta.data %>%
  mutate(
    epi_cell_type6 = case_when(
      epi_cell_type3 == "Tumor-ISC-like" ~ "Canonical_CRC_Stem",
      epi_cell_type3 == "Intestine-like" ~ "Canonical_CRC_Intestine",
      epi_cell_type3 == "Proliferative tumor" ~ "Canonical_CRC_Intestine_Proliferation",
      epi_cell_type3 == "Proliferative stem-like tumor" ~ "Canonical_CRC_Stem_Proliferation",
      epi_cell_type3 == "Angiogenic tumor" ~ "Non_Canonical_CRC_1",
      epi_cell_type3 == "APCDD1+ tumor" ~ "Non_Canonical_CRC_2",
      TRUE ~ as.character(epi_cell_type3)
    )
  )

epithelial_final$epi_cell_type6 <- factor(
  epithelial_final$epi_cell_type6,
  levels = rev(c(
    "Non_Canonical_CRC_1", "Non_Canonical_CRC_2",
    "Canonical_CRC_Stem", "Canonical_CRC_Intestine",
    "Canonical_CRC_Stem_Proliferation",
    "Canonical_CRC_Intestine_Proliferation",
    "Stem cells", "Transit-amplifying cells", "Enterocytes",
    "Goblet cells", "Tuft cells", "Enteroendocrine-like cells"
  ))
)

# Optional convenience alias used in the manuscript narrative.
epithelial_final@meta.data <- epithelial_final@meta.data %>%
  mutate(
    CMETS_status = case_when(
      epi_cell_type6 == "Non_Canonical_CRC_1" ~ "CMETS",
      epi_cell_type6 %in% c(
        "Canonical_CRC_Stem",
        "Canonical_CRC_Intestine",
        "Canonical_CRC_Stem_Proliferation",
        "Canonical_CRC_Intestine_Proliferation",
        "Non_Canonical_CRC_2"
      ) ~ "Other_CRC_tumor_state",
      TRUE ~ "Non_tumor_epithelial"
    )
  )

# ==============================================================================
# 6. Save final object, metadata, and compact QC plots
# ==============================================================================

saveRDS(epithelial_final, outputs$epithelial_final_rds)

metadata_dt <- as.data.table(epithelial_final@meta.data, keep.rownames = "cell_id")

metadata_cols <- intersect(
  c(
    "cell_id", "orig.ident", "Patient_ID", "Site_of_Origin", "Tissue_Type",
    "Primary_Side", "MSI", "Tx_in_6mo",
    "epithelial_cluster_sub3",
    "epi_cell_type", "epi_cell_type2", "epi_cell_type3", "epi_cell_type6", "CMETS_status"
  ),
  colnames(metadata_dt)
)

fwrite(metadata_dt[, ..metadata_cols], outputs$metadata_csv)

epi_cell_type6_colors <- c(
  "Non_Canonical_CRC_1" = "#be0032",
  "Non_Canonical_CRC_2" = "#008856",
  "Canonical_CRC_Stem" = "#e68fac",
  "Canonical_CRC_Intestine" = "#f99379",
  "Canonical_CRC_Stem_Proliferation" = "#8db600",
  "Canonical_CRC_Intestine_Proliferation" = "#f3c300",
  "Stem cells" = "khaki",
  "Transit-amplifying cells" = "#a1caf1",
  "Enterocytes" = "#882d17",
  "Goblet cells" = "#c2b280",
  "Tuft cells" = "#f99379",
  "Enteroendocrine-like cells" = "#2b3d26"
)

# Use the final UMAP if available; otherwise fall back to the earlier epithelial UMAP.
plot_reduction <- if ("epithelial_umap.scvi" %in% names(epithelial_final@reductions)) {
  "epithelial_umap.scvi"
} else if ("epithelial_umap2.scvi" %in% names(epithelial_final@reductions)) {
  "epithelial_umap2.scvi"
} else {
  stop("No epithelial UMAP reduction found.", call. = FALSE)
}

plot_dim(
  epithelial_final,
  group_by = "epi_cell_type6",
  reduction = plot_reduction,
  out_pdf = outputs$final_umap_pdf,
  width = 7,
  height = 6,
  colors = epi_cell_type6_colors
)

DefaultAssay(epithelial_final) <- "RNA"

marker_genes <- c(
  "EMP1", "KRT20", "VEGFA", "APCDD1", "PROX1", "LGR5", "RGMB",
  "MKI67", "TOP2A", "CLCA4", "SLC26A3", "MUC2", "CLCA1",
  "SH2D6", "TRPM5", "CHGA", "CHGB"
)

marker_genes <- intersect(marker_genes, rownames(epithelial_final))

if (length(marker_genes) > 0) {
  rdylbu_colors <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(10))

  p_dot <- DotPlot(
    epithelial_final,
    group.by = "epi_cell_type6",
    features = marker_genes
  ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    coord_fixed() +
    scale_color_gradientn(
      colors = rdylbu_colors,
      limits = c(-3, 3),
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    scale_size_area(limits = c(0, 100), oob = scales::squish)

  pdf(outputs$final_dotplot_pdf, width = 8, height = 8)
  print(p_dot)
  dev.off()
}

message("Final epithelial annotation complete.")
message("Final RDS: ", outputs$epithelial_final_rds)
message("Final metadata: ", outputs$metadata_csv)
