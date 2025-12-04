import os
import argparse
import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def stagger_spatial_coordinates_dynamic(adatas):
    """
    Stagger spatial coordinates for multiple AnnData objects.

    Parameters:
        adatas (list of AnnData): List of AnnData objects to stagger spatial coordinates.
    
    Returns:
        list of AnnData: List of AnnData objects with staggered spatial coordinates.
    """
    staggered_adatas = []
    x_offset = 0

    for adata in adatas:
        # Get spatial dimensions
        if 'spatial' in adata.obsm:
            spatial = adata.obsm['spatial']

        # Apply staggered offsets
        staggered_adata = stagger_spatial_coordinates(adata, x_offset=x_offset, y_offset=0)
        staggered_adatas.append(staggered_adata)

        # Update offsets for the next dataset
        x_offset += spatial[:, 0].max() + 1000  # Add a buffer between samples

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

    args = parser.parse_args()

    # Define output file paths
    output_file = f"{args.output_prefix}_merged.h5ad"
    plot_file_png = f"{args.output_prefix}_spatial_plot.png"

    # Load genes to remove if provided
    genes_to_remove = []
    if args.remove_genes_tsv:
        genes_to_remove = pd.read_csv(args.remove_genes_tsv, header=None)[0].tolist()
        print(f'genes to remove: {genes_to_remove}')

    # Load and process AnnData objects
    adatas = load_all_anndata(args.input_dir)
    if genes_to_remove:
        adatas = [(name, remove_genes(adata, genes_to_remove)) for name, adata in adatas]
    staggered_adatas = stagger_spatial_coordinates_dynamic([adata for _, adata in adatas])
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

