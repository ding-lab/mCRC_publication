#!/bin/bash

# Define input directory, output file, and plot file paths
INPUT_DIR="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_anndata"
OUTPUT_FILE="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_merged_anndata/mCRC_Xenium_5K_21"
LOG_FILE="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_merged_anndata/merge_and_stagger.log"
GENE_FILE="/diskmnt/Projects/Users/Evan.p/scripts/Python/Banksy/5K_gene_to_remove.tsv"

# Activate the Python environment
eval "$(conda shell.bash hook)"
conda activate scanpy_env

# Run the Python script
python ./src/Xenium_anndata_merge_v2.py \
    --input_dir "$INPUT_DIR" \
    --output_prefix "$OUTPUT_FILE" \
    --remove_genes_tsv "$GENE_FILE" \
    --samples_per_row 6 \
    --grid_width 15000 \
    --grid_height 15000 \
    &> "$LOG_FILE"

# Deactivate the Python environment
conda deactivate
