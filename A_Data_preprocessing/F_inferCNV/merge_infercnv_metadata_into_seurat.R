# Merge infercnv_cell_metadata*.csv into a (merged) Seurat object.
# Cell names in Seurat must match the `cell` column (InferCNV naming: LIBRARY_BARCODE-1).
# Wide tables: infercnv_cell_metadata.csv (20241001) or build_infercnv_cell_metadata_20260403.sh (+ infercnv_cnv_call).

infercnv_meta <- read.csv(
  "/diskmnt/Users2/epeng/Projects/mCRC/Revesion/infercnv_cell_metadata.csv",
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
rownames(infercnv_meta) <- infercnv_meta$cell

# Example (uncomment and adapt):
# common <- intersect(Cells(seurat_merged), rownames(infercnv_meta))
# seurat_merged <- AddMetaData(
#   seurat_merged,
#   metadata = infercnv_meta[common, setdiff(names(infercnv_meta), "cell")]
# )
# message("Cells with InferCNV metadata: ", length(common), " / ", ncol(seurat_merged))
