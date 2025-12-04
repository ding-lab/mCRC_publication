import argparse
import pandas as pd
import spatialdata as sd
from spatialdata_io import xenium
import anndata as ad
import scanpy as sc
import logging
import os

def setup_logging(output_dir):
    log_file = f"{output_dir}/xenium_anndata.log"
    logging.basicConfig(
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Console output
        ],
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )
    logging.info("Logging initialized.")


def remove_cells(adata, cell_removal_csv):
    """Remove specific cells from an AnnData object based on a CSV file."""
    try:
        # Read the CSV file and skip the first three rows
        removal_data = pd.read_csv(cell_removal_csv, skiprows=3, header=None)
        
        # The first column contains the cell IDs to remove
        cells_to_remove = removal_data.iloc[:, 0].str.strip().tolist()
        
        # Ensure valid cell IDs exist in AnnData
        valid_cells_to_remove = [cell for cell in cells_to_remove if cell in adata.obs_names]

        if len(valid_cells_to_remove) > 0:
            # Remove the cells from AnnData
            adata = adata[~adata.obs_names.isin(valid_cells_to_remove)].copy()
            logging.info(f"Removed {len(valid_cells_to_remove)} cells from the AnnData object.")
        else:
            logging.info("No matching cells found to remove.")

    except Exception as e:
        logging.error(f"Error removing cells: {e}")

    return adata

def generate_anndata(csv_file_path, xenium_path_column, output_dir, filter_column=None, cell_removal_tsv=None):
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

    # Load the cell removal information from the TSV file
    removal_info = None
    if cell_removal_tsv:
        try:
            removal_info = pd.read_csv(cell_removal_tsv, sep="\t", header=None, names=["sample", "removal_csv"])
            logging.info("Cell removal TSV file loaded successfully.")
        except Exception as e:
            logging.error(f"Error loading cell removal TSV file: {e}")

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

            if not os.path.exists(zarr_path):
                logging.info(f"Writing Zarr output to {zarr_path}.")
                sdata.write(zarr_path)
            else:
                logging.info(f"Zarr output already exists at {zarr_path}. Skipping.")

            # Extract the AnnData table
            adata = sdata.tables["table"]
            adata.obs["dataset"] = sample_name
            adata.obs["casette"] = Casette_id
            adata.obs["barcode"] = adata.obs["cell_id"]
            adata.obs.set_index('barcode', inplace=True)
            logging.info(f"Sample {sample_name} AnnData shape: {adata.shape}")

            # Perform QC and filter cells and genes
            sc.pp.filter_cells(adata, min_counts=20)
            sc.pp.filter_cells(adata, min_genes=10)
            logging.info(f"Sample {sample_name} post-QC AnnData shape: {adata.shape}")

            # Check for cell removal
            if removal_info is not None and sample_name in removal_info["sample"].values:
                cell_removal_csv = removal_info.loc[removal_info["sample"] == sample_name, "removal_csv"].values[0]
                if os.path.exists(cell_removal_csv):
                    logging.info(f"Removing cells for sample {sample_name} using {cell_removal_csv}.")
                    adata = remove_cells(adata, cell_removal_csv)
            else:
                logging.warning(f"Cell removal CSV file not found for sample {sample_name}.")

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
    parser.add_argument("--cell_removal_tsv", help="Path to the TSV file specifying cells to remove.")

    args = parser.parse_args()

    print("Input Arguments:")
    print(f"CSV File: {args.csv_file}")
    print(f"Xenium Path Column: {args.Xenium_path}")
    print(f"Output Directory: {args.output_dir}")
    print(f"Filter Column: {args.filter_column}")
    print(f"Cell Removal TSV: {args.cell_removal_tsv}")

    

    # Call the function with the parsed arguments
    generate_anndata(args.csv_file, args.Xenium_path, args.output_dir, args.filter_column, args.cell_removal_tsv)

if __name__ == "__main__":
    main()