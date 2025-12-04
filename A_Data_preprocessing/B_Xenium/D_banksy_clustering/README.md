# Python Pipeline for Harmonized Banksy Clustering

This repository provides a step-by-step pipeline for processing Xenium data, merging datasets, and performing clustering using Banksy and Harmony. This pipeline is tailored for the computing environment in Ding Lab at WashU.

For the latest one-click solution, please see Simon’s repository:
https://github.com/simoncmo/PyBanksy-Harmony-master/tree/main

## Environment
/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy \
Python version of BANKSY: https://github.com/prabhakarlab/Banksy_py

---

## **Pipeline Overview**

### **Step 1: Generate Anndata for Each Xenium Sample**
Use the provided script to generate anndata objects from individual Xenium samples.

Run the following command:
```bash
bash run_making_anndata.sh
```
### **Step 2: Merge anndata and Stagger the Coordinates**
Merge the generated anndata objects and stagger their spatial coordinates. Optionally, remove unneeded genes by specifying a gene list.

Run the following command:
```bash
bash run_merge_anndata.sh
```

#### **Arguments**

##### **1. Required Arguments**
- `--input_dir`  
  - **Type**: String  
  - **Description**: Path to the folder containing `.h5ad` files that need to be merged.  
  - **Required**: Yes  

- `--output_prefix`  
  - **Type**: String  
  - **Description**: Prefix for the output files. All resulting files will use this prefix.  
  - **Required**: Yes  

---

##### **2. Optional Arguments**
- `--remove_genes_tsv`  
  - **Type**: String  
  - **Description**: Path to a `.tsv` file containing a list of genes to remove from the data. If not provided, no genes are removed.  
  - **Required**: No  

---

##### **3. Configuration Arguments**
- `--samples_per_row`  
  - **Type**: Integer  
  - **Default**: `4`  
  - **Description**: The number of samples to display per row in the spatial plot. Helps in organizing the visualization layout.

- `--grid_width`  
  - **Type**: Integer  
  - **Default**: `5000`  
  - **Description**: The width of each grid block when staggering spatial coordinates. Controls the spacing between samples.

- `--grid_height`  
  - **Type**: Integer  
  - **Default**: `5000`  
  - **Description**: The height of each grid block when staggering spatial coordinates. Controls the vertical spacing between samples.



### **Step 3: Run Banksy and Harmony, Then Perform Clustering**
Run Banksy and Harmony for data integration, followed by clustering.

Run the following command:
```bash
bash run_Banksy.sh
```

#### Banksy Harmonized Pipeline: Argument Descriptions

The following arguments are used in the `Banksy Harmonized Pipeline` script:

##### **Required Arguments**
- `--input_merged_anndata`:  
  Path to the input merged AnnData file. This is the primary input for the pipeline.
  - **Type**: File path
  - **Required**: Yes

- `--output_dir`:  
  Directory where the results will be saved.
  - **Type**: Directory path
  - **Required**: Yes

---

##### **Optional Arguments**
- `--output_prefix`:  
  Prefix for naming all output files. Default is `banksy`.
  - **Type**: String
  - **Default**: `"banksy"`

- `--n_top_genes`:  
  Number of top highly variable genes to retain for downstream analysis.
  - **Type**: Integer
  - **Default**: `2000`

- `--k_geom`:  
  Specifies the `K` parameter for Banksy initialization geometry.
  - **Type**: Integer
  - **Default**: `15`

- `--max_m`:  
  Maximum order for azimuthal transform (m-th order). Default is 1.
  - **Type**: Integer
  - **Default**: `1`

- `--nbr_weight_decay`:  
  Method for neighbor weight decay. Choose from:
  - `"scaled_gaussian"`: Scaled Gaussian decay
  - `"reciprocal"`: Reciprocal decay
  - `"uniform"`: Uniform weights
  - `"ranked"`: Ranked decay
  - **Type**: String
  - **Choices**: `"scaled_gaussian"`, `"reciprocal"`, `"uniform"`, `"ranked"`
  - **Default**: `"scaled_gaussian"`

- `--pca_dims`:  
  Dimensionality for PCA reduction. Can specify multiple dimensions.
  - **Type**: List of integers
  - **Default**: `[20]`

- `--lambda_list`:  
  List of lambda parameters for Banksy optimization.
  - **Type**: List of floats
  - **Default**: `[0.8]`

- `--harmony_batch_key`:  
  Column name in the `AnnData` object used for Harmony batch correction.
  - **Type**: String
  - **Default**: `"dataset"`

---

##### **Clustering Arguments**
- `--run_clustering`:  
  Leiden clustering is very slow for a dataset with million cells. Secuer clustering (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010753) is an alternative ultrafast algorithm to save time. \
  Specifies clustering method(s) to use. Options:
  - `"leiden"`: Run Leiden clustering.
  - `"secuer"`: Run Secuer clustering.
  - `"both"`: Run both methods.
  - **Type**: String
  - **Choices**: `"leiden"`, `"secuer"`, `"both"`
  - **Default**: `"both"`

- `--leiden_resolution`:  
  Resolution parameter for Leiden clustering.
  - **Type**: Float
  - **Default**: `0.5`

- `--secuer_resolution`:  
  Resolution parameter for Secuer clustering.
  - **Type**: Float
  - **Default**: `1.0`

---

##### **Additional Options**
- `--sample_id_column`:  
  Column name in `adata.obs` that represents the sample ID. Used for data grouping.
  - **Type**: String
  - **Default**: `"dataset"`

- `--plot`:  
  Generate spatial scatter plots. Include this flag to enable plotting.
  - **Type**: Boolean flag
  - **Default**: `False`

