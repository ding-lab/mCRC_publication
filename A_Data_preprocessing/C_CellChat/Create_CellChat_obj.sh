Rscript /diskmnt/Users2/epeng/Projects/mCRC/scripts/spatial_omics/Xenium/CellChat/Create_CellChat_obj.R \
  --seurat /diskmnt/Projects/MetNet_analysis/Colorectal/data/Xenium_objects/mCRC_Xenium_N26_liver_TSI2.qs \
  --annotation All_cell_type3 \
  --subset_annotation Broad_cell_type1 \
  --cell_types Tumor T_NK_cell Myeloid_cells B_cell Fibroblast Endothelial_cell Cholangiocyte Hepatocyte \
  --filterDB FALSE \
  --signal "Secreted Signaling" \
  --output_dir /diskmnt/Projects/MetNet_analysis/Colorectal/data/cellchat_objects/by_tmr_regions/ \
  --output_file liver_TSI_cellchat_ver2