import argparse
import pandas as pd
import spatialdata as sd
from spatialdata_io import xenium
import anndata as ad
import scanpy as sc
import logging

def setup_logging(output_dir):
    logging.basicConfig(
        filename=f"{output_dir}/xenium_anndata.log",
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )
    logging.info("Logging initialized.")

def generate_anndata(csv_file_path, xenium_path_column, output_dir, filter_column=None):
    setup_logging(output_dir)

    try:
        # Load the CSV file into a DataFrame
        data = pd.read_csv(csv_file_path)
        logging.info("CSV file loaded successfully.")
    except Exception as e:
        logging.error(f"Error loading CSV file: {e}")
        return

    # Ensure the CSV is not empty
    if data.empty:
        logging.error("The CSV file is empty. Please check the file content.")
        return

    # Strip whitespace from column names and values
    data.columns = data.columns.str.strip()

    # Check if the specified column exists
    if xenium_path_column not in data.columns:
        logging.error(f"Column '{xenium_path_column}' not found in the CSV file. Available columns: {data.columns.tolist()}")
        return

    # Apply filtering if filter_column is specified
    if filter_column:
        if filter_column not in data.columns:
            logging.error(f"Filter column '{filter_column}' not found in the CSV file. Available columns: {data.columns.tolist()}")
            return
        data = data[data[filter_column].str.strip() == "Yes"]
        if data.empty:
            logging.warning("No rows match the filter criteria.")
            return

    # Iterate through each sample in the specified column
    for index, row in data.iterrows():
        try:
            xenium_path = row[xenium_path_column]
            sample_name = row["Section ID"].strip()  # Use 'Section ID' as the file name
            Casette_id = row["Casette ID"].strip()

            # Derive the output Zarr path
            zarr_path = f"{output_dir}/{sample_name}.zarr"
            h5ad_path = f"{output_dir}/{sample_name}.h5ad"

            logging.info(f"Processing sample: {sample_name}")

            # Generate SpatialData object from Xenium path
            sdata = xenium(xenium_path)
            
            # Write the Zarr output
            sdata.write(zarr_path)

            # Extract the AnnData table
            adata = sdata.tables["table"]
            adata.obs["dataset"] = sample_name
            adata.obs["casette"] = Casette_id
            logging.info(f"Sample {sample_name} AnnData shape: {adata.shape}")

            # Perform QC and filter cells and genes
            sc.pp.filter_cells(adata, min_counts=20, min_genes=10)
            logging.info(f"Sample {sample_name} post-QC AnnData shape: {adata.shape}")

            # Save the AnnData table to H5AD format
            adata.write(h5ad_path)
            logging.info(f"Sample {sample_name} saved as {h5ad_path}")

        except Exception as e:
            logging.error(f"Error processing sample {row['Section ID']}: {e}")

        # Additional processing if needed
        # ...

def main():
    parser = argparse.ArgumentParser(description="Generate Xenium AnnData from paths in a CSV file.")
    parser.add_argument("--csv_file", required=True, help="Path to the CSV file.")
    parser.add_argument("--Xenium_path", required=True, help="Name of the column containing Xenium output paths.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the outputs.")
    parser.add_argument("--filter_column", help="Name of the column to filter rows where the value is 'Yes'.")

    args = parser.parse_args()

    # Call the function with the parsed arguments
    generate_anndata(args.csv_file, args.Xenium_path, args.output_dir, args.filter_column)

if __name__ == "__main__":
    main()