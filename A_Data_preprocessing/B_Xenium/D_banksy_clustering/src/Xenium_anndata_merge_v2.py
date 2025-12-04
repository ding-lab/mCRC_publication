import os
import argparse
import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def stagger_spatial_coordinates_grid(adatas, samples_per_row=4, grid_width=10000, grid_height=10000):
    """
    Stagger spatial coordinates for multiple AnnData objects into a fixed grid layout.

    Parameters:
        adatas (list of AnnData): List of AnnData objects to stagger spatial coordinates.
        samples_per_row (int): Number of samples to arrange in a row.
        grid_width (float): Width of each grid block assigned to a sample.
        grid_height (float): Height of each grid block assigned to a sample.

    Returns:
        list of AnnData: List of AnnData objects with staggered spatial coordinates.
    """
    staggered_adatas = []

    for i, adata in enumerate(adatas):
        # Determine the grid position
        row, col = divmod(i, samples_per_row)
        x_offset = col * grid_width
        y_offset = row * grid_height

        # Center the sample within its grid block
        spatial = adata.obsm['spatial']
        x_center = spatial[:, 0].mean()
        y_center = spatial[:, 1].mean()

        # Calculate shifts to center the sample in the grid block
        x_shift = x_offset + grid_width / 2 - x_center
        y_shift = y_offset + grid_height / 2 - y_center

        # Apply the calculated shifts
        staggered_adata = stagger_spatial_coordinates(adata, x_offset=x_shift, y_offset=y_shift)
        staggered_adatas.append(staggered_adata)

    print(f"Samples arranged in a grid with {samples_per_row} samples per row.")
    return staggered_adatas


def stagger_spatial_coordinates(adata, x_offset=0, y_offset=0):
    """
    Stagger spatial coordinates stored in `adata.obsm['spatial']`.

    Parameters:
        adata (AnnData): The AnnData object with spatial coordinates.
        x_offset (float): Offset to add to the x-coordinate.
        y_offset (float): Offset to add to the y-coordinate.

    Returns:
        AnnData: A new AnnData object with staggered spatial coordinates.
    """
    adata_copy = adata.copy()
    if 'spatial' in adata_copy.obsm:
        spatial = adata_copy.obsm['spatial'].copy()
        spatial[:, 0] += x_offset
        spatial[:, 1] += y_offset
        adata_copy.obsm['spatial'] = spatial
    else:
        raise KeyError("'spatial' not found in obsm. Ensure spatial data is available.")
    return adata_copy


def load_all_anndata(folder_path):
    """
    Load all AnnData objects from a folder and add file name as dataset.

    Parameters:
        folder_path (str): Path to the folder containing .h5ad files.

    Returns:
        list of tuple: List of (file name, AnnData) tuples.
    """
    adatas = []
    for file in os.listdir(folder_path):
        if file.endswith(".h5ad"):
            file_path = os.path.join(folder_path, file)
            adata = sc.read_h5ad(file_path)
            dataset_name = file.replace(".h5ad", "")
            adata.obs["dataset"] = dataset_name  # Add dataset name as a column
            adatas.append((dataset_name, adata))
    return adatas


def remove_genes(adata, genes_to_remove):
    """
    Remove specified genes from an AnnData object.

    Parameters:
        adata (AnnData): The AnnData object.
        genes_to_remove (list): List of gene names to remove.

    Returns:
        AnnData: AnnData object with the specified genes removed.
    """
    genes_to_keep = [gene for gene in adata.var_names if gene not in genes_to_remove]
    return adata[:, genes_to_keep].copy()


def merge_anndata(adatas):
    """
    Merge multiple AnnData objects.

    Parameters:
        adatas (list of tuple): List of (dataset name, AnnData) tuples.

    Returns:
        AnnData: Merged AnnData object.
    """
    keys, adata_list = zip(*adatas)
    merged_adata = ad.concat(
        adata_list,
        join="outer",
        label="dataset",
        keys=keys,
        index_unique="_"
    )

    # Add x and y coordinates to obs
    if 'spatial' in merged_adata.obsm:
        spatial_coords = merged_adata.obsm['spatial']
        merged_adata.obs['xcoord'] = spatial_coords[:, 0]
        merged_adata.obs['ycoord'] = spatial_coords[:, 1]
    else:
        raise KeyError("'spatial' not found in obsm of merged AnnData.")

    return merged_adata


def main():
    parser = argparse.ArgumentParser(description="Merge and stagger spatial AnnData objects.")
    parser.add_argument("--input_dir", required=True, help="Path to the folder containing .h5ad files.")
    parser.add_argument("--output_prefix", required=True, help="Prefix for the output files.")
    parser.add_argument("--remove_genes_tsv", required=False, help="Path to the TSV file listing genes to remove (optional).")
    parser.add_argument("--samples_per_row", type=int, default=4, help="Number of samples per row in the spatial plot.")
    parser.add_argument("--grid_width", type=int, default=5000, help="Width of each grid block.")
    parser.add_argument("--grid_height", type=int, default=5000, help="Height of each grid block.")

    args = parser.parse_args()

    # Define output file paths
    output_file = f"{args.output_prefix}_merged.h5ad"
    plot_file_png = f"{args.output_prefix}_spatial_plot.png"

    # Load genes to remove if provided
    genes_to_remove = []
    if args.remove_genes_tsv:
        # Read the TSV file to get a list of genes to remove
        if not os.path.exists(args.remove_genes_tsv):
            raise FileNotFoundError(f"Gene list file not found: {args.remove_genes_tsv}")
        genes_to_remove = pd.read_csv(args.remove_genes_tsv, header=None)[0].tolist()
        print(f"Genes to remove: {genes_to_remove}")

    # Load and process AnnData objects
    adatas = load_all_anndata(args.input_dir)
    if genes_to_remove:
        # Remove specified genes from each AnnData object
        adatas = [(name, remove_genes(adata, genes_to_remove)) for name, adata in adatas]

    # Stagger spatial coordinates and merge AnnData objects
    staggered_adatas = stagger_spatial_coordinates_grid(
        [adata for _, adata in adatas],
        samples_per_row=args.samples_per_row,
        grid_width=args.grid_width,
        grid_height=args.grid_height
    )
    merged_adata = merge_anndata(list(zip([key for key, _ in adatas], staggered_adatas)))

    # Save the merged AnnData
    merged_adata.write(output_file)
    print(f"Merged AnnData saved to {output_file}")

    # Plot the merged spatial coordinates and save to PNG
    plt.figure()
    sq.pl.spatial_scatter(
        merged_adata,
        library_id="spatial",
        shape=None,
        color=["dataset"],
        wspace=0.4,
        save=None  # Disable automatic saving to avoid appending to the path
    )
    plt.savefig(plot_file_png, dpi=300, bbox_inches="tight")  # Save explicitly with high resolution
    plt.close()
    print(f"Spatial plot saved to {plot_file_png}")


if __name__ == "__main__":
    main()
