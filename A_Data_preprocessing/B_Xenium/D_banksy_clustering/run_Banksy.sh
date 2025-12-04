#!/bin/bash

# Define paths and parameters
INPUT_FILE="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_merged_anndata/mCRC_Xenium_5K_21_merged.h5ad"  # Path to input merged AnnData file
OUTPUT_DIR="/diskmnt/Projects/MetNet_analysis_2/Colorectal/Xenium/Banksy/mCRC_Xenium_5K_N21/VG2000_K15M0L8PC20leR5"          # Directory to save output files
OUTPUT_PREFIX="mCRC_Xenium_5K_N21_VG2000_K15M0L8PC20leR5"                 # Prefix for all output files
LOG_FILE="${OUTPUT_DIR}/banksy_pipeline.log"    # Path to log file
CONDA_ENV="scanpy_env"                          # Name of the Conda environment to use

# Define Banksy parameters
N_TOP_GENES=2000
K_GEOM=15
MAX_M=0
NBR_WEIGHT_DECAY="scaled_gaussian"
PCA_DIMS=20
LAMBDA_VALUE=0.8
HARMONY_BATCH_KEY="dataset"
RUN_CLUSTERING="both"                           # Options: "leiden", "secuer", or "both"
LEIDEN_RESOLUTION=0.5
SECUER_RESOLUTION=1.2
SAMPLE_ID_COLUMN="dataset"
PLOT="true"                                     # Set to "true" to generate a spatial scatter plot

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# Run the Python script
python ./src/Xenium_Banksy_Harmony_v1.py \
    --input_merged_anndata "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --output_prefix "$OUTPUT_PREFIX" \
    --n_top_genes "$N_TOP_GENES" \
    --k_geom "$K_GEOM" \
    --max_m "$MAX_M" \
    --nbr_weight_decay "$NBR_WEIGHT_DECAY" \
    --pca_dims "$PCA_DIMS" \
    --lambda_list "$LAMBDA_VALUE" \
    --harmony_batch_key "$HARMONY_BATCH_KEY" \
    --run_clustering "$RUN_CLUSTERING" \
    --leiden_resolution "$LEIDEN_RESOLUTION" \
    --secuer_resolution "$SECUER_RESOLUTION" \
    --sample_id_column "$SAMPLE_ID_COLUMN" \
    $( [ "$PLOT" = "true" ] && echo "--plot" ) \
    &> "$LOG_FILE"

# Deactivate Conda environment
conda deactivate

# Print completion message
echo "Banksy pipeline completed. Check the log file at: $LOG_FILE"
