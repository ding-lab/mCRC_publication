#!/bin/bash

# Set the input arguments
CSV_FILE="/diskmnt/Projects/MetNet_analysis_2/Colorectal/Xenium/Tracking_sheet/mCRC_tracking_Xenium_20250304.csv"
XENIUM_PATH_COLUMN="Output file path"
OUTPUT_DIR="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_anndata"
FILTER_COLUMN="New_processing"
FILTER_CELLS="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium_python/cells_to_remove/mCRC_Xenium_5K_cell_to_remove.tsv"

mkdir -p "$OUTPUT_DIR"

# Activate the Python environment
eval "$(conda shell.bash hook)"
conda activate scanpy_env

# Run the Python script
python ./src/Xenium_anndata_generator_v2.py \
    --csv_file "$CSV_FILE" \
    --Xenium_path "$XENIUM_PATH_COLUMN" \
    --output_dir "$OUTPUT_DIR" \
    --filter_column "$FILTER_COLUMN" \
    --cell_removal_tsv "$FILTER_CELLS"

# Deactivate the Python environment
conda deactivate
