# This script creates a slim seurat object from Xenium outputs
eval "$(conda shell.bash hook)"
conda activate "/diskmnt/Projects/Users/simonmo/Softwares/miniconda3/envs/seurat521"

# Setup parameters
## -------------------- Modify BELOW per analysis ------------------------##
# Parameters - Paths
SAMPLE_TABLE="./mCRC_Xenium_sample_table.tsv"
OUTPUT_DIR="./merged_objects"
ANALYSIS_NAME="mCRC_Xenium_N26"

# Parameters - 1-create object
shared_genes_only=FALSE # TRUE or FALSE. If FALSE, all genes are kept
remove_microbiome=TRUE # TRUE or FALSE. If TRUE, remove microbiome genes
min_transcripts=10 # minimum number of transcripts per gene. Default to 5. Set 10 here
## -------------------- Modify ABOVE per analysis ------------------------##

# Validate parameters
if [ ! -f "$SAMPLE_TABLE" ]; then
    echo "ERROR: SAMPLE_TABLE not found: $SAMPLE_TABLE"
    exit 1
fi
mkdir -p $OUTPUT_DIR

# create folder structure
step1_out="${OUTPUT_DIR}/${ANALYSIS_NAME}/1_createObject"
mkdir -p $step1_out

# scripts
script1_create_object="./create_slim_object_layer_v4.R"

# step1 output
merged_obj_path="${step1_out}/1_xenium_slim_${ANALYSIS_NAME}_normalized.qs" 

# 1. create object
# check if the object already exists
if [ -f "$merged_obj_path" ]; then
    echo "INFO: Object already exists: $merged_obj_path"
    echo "Skip creating object"
else
    echo "INFO: Creating object: $merged_obj_path"
    Rscript $script1_create_object \
        --sample_table "$SAMPLE_TABLE" \
        --out_dir "$step1_out" \
        --analysis_name "$ANALYSIS_NAME" \
        --shared_genes_only "$shared_genes_only" \
        --remove_microbiome "$remove_microbiome" \
        --min_transcripts "$min_transcripts" \
        &> "${step1_out}/run.log"
fi

echo "Object creation complete: $merged_obj_path"
