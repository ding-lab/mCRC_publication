import os
import argparse
import scanpy as sc
import secuer as sr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append('/diskmnt/Projects/Users/Evan.p/tools/Banksy_py')
from banksy_utils.filter_utils import normalize_total, filter_hvg
from banksy.main import median_dist_to_nearest_neighbour
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy_utils.pca_harmony_umap import pca_harmony_umap

def main(args):
    try:
        # Create output directory
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Step 1: Load data
        adata = sc.read_h5ad(args.input_merged_anndata)
        print("Loaded AnnData object:", adata)
        
        # Step 2: Normalize and filter
        adata = normalize_total(adata)
        adata, adata_allgenes = filter_hvg(adata, n_top_genes=args.n_top_genes, flavor="seurat")
        
        # Save intermediate result
        adata.write(os.path.join(args.output_dir, f"{args.output_prefix}_normalized_filtered.h5ad"))
        
        # Step 3: Calculate median distance and initialize Banksy
        nbrs = median_dist_to_nearest_neighbour(adata, key='spatial')
        coord_keys=('xcoord', 'ycoord', 'spatial')
        k_geom=args.k_geom
        banksy_dict = initialize_banksy(
            adata,
            coord_keys,
            k_geom,
            nbr_weight_decay=args.nbr_weight_decay,
            max_m=args.max_m,
            plt_edge_hist=False,
            plt_nbr_weights=False,
            plt_agf_angles=False,
            plt_theta=False
            )
        
        # Step 4: Generate Banksy matrix
        harmony_checkpoint = os.path.join(args.output_dir, f"{args.output_prefix}_harmony.pkl")
        if os.path.exists(harmony_checkpoint):
            print(f"Loading Banksy dictionary with Harmony applied from {harmony_checkpoint}")
            banksy_dict = pd.read_pickle(harmony_checkpoint)
        else:
            banksy_dict, banksy_matrix = generate_banksy_matrix(
                adata, banksy_dict, args.lambda_list, args.max_m
            )
            pca_harmony_umap(
                banksy_dict=banksy_dict,
                pca_dims=args.pca_dims,
                plt_remaining_var=False,
                add_umap=True,
                harmony_batch_key=args.harmony_batch_key
            )
            # pd.to_pickle(banksy_dict, harmony_checkpoint)
            # print(f"Saved Banksy dictionary with Harmony applied to {harmony_checkpoint}")

        # Step 5: Extract and save Banksy matrix
        banksy_matrix2 = banksy_dict['scaled_gaussian'][args.lambda_list[0]]['adata']
        banksy_matrix2.obsm['spatial'] = adata.obsm['spatial']
        banksy_matrix2.obs['dataset__cell_id'] = banksy_matrix2.obs['dataset'].astype(str) + '__' + banksy_matrix2.obs['cell_id'].astype(str)

        banksy_matrix_file = os.path.join(args.output_dir, f"{args.output_prefix}_banksy_matrix.h5ad")
        if not os.path.exists(banksy_matrix_file):
            banksy_matrix2.write(banksy_matrix_file)
            print(f"Saved Banksy matrix to {banksy_matrix_file}")

        num_components = banksy_matrix2.obsm['X_pca_harmony'].shape[1]
        column_names = [f"Banksy_harmony_PC{i+1}" for i in range(num_components)]
        print(column_names)

        pca_harmony_df = pd.DataFrame(
            banksy_matrix2.obsm['X_pca_harmony'], 
            index=banksy_matrix2.obs['dataset__cell_id'], 
            columns=column_names)

        banksy_reduction_file = os.path.join(args.output_dir, f"{args.output_prefix}_banksy_reduction.csv")
        pca_harmony_df.to_csv(banksy_reduction_file)
        print(f"Banksy reduction file saved to {banksy_reduction_file}")

        num_umap_components = banksy_matrix2.obsm['X_pca_harmony_umap'].shape[1]
        umap_column_names = [f"UMAP{i+1}" for i in range(num_umap_components)]
        print(umap_column_names)

        pca_harmony_umap_df = pd.DataFrame(
        banksy_matrix2.obsm['X_pca_harmony_umap'], 
        index=banksy_matrix2.obs['dataset__cell_id'], 
        columns=umap_column_names)

        banksy_umap_file = os.path.join(args.output_dir, f"{args.output_prefix}_banksy_umap_reduction.csv")
        pca_harmony_umap_df.to_csv(banksy_umap_file)
        print(f"Banksy umap reduction file saved to {banksy_umap_file}")

        # Save final annotated matrix
        final_output = os.path.join(args.output_dir, f"{args.output_prefix}_final_banksy_matrix.h5ad")
        banksy_matrix2.write(final_output)
        print(f"Saved final Banksy matrix to {final_output}")
            
        
        # Step 6: Clustering
        if args.run_clustering in ["secuer", "both"]:
            fea = banksy_matrix2.obsm['X_pca_harmony']
            res = sr.secuer(fea=fea, Knn=7, eskResolution=args.secuer_resolution, num_multiProcesses=20)
            banksy_matrix2.obs['secuer_cluster'] = res

            secuer_dir = os.path.join(args.output_dir, "secuer_clusters")
            os.makedirs(secuer_dir, exist_ok=True)

            for sample_id in banksy_matrix2.obs[args.sample_id_column].unique():
                # Filter the obs dataframe by the sample ID
                subset_obs = banksy_matrix2.obs[banksy_matrix2.obs[args.sample_id_column] == sample_id].copy()
                
                # Rename the 'secuer_cluster' column to 'group'
                subset_obs.rename(columns={'secuer_cluster': 'group'}, inplace=True)
                
                # Define the output file name
                output_file = os.path.join(secuer_dir, f"{sample_id}_secuer_obs.csv")
                
                # Save the subset to a CSV file
                subset_obs.to_csv(output_file, index=True)
                print(f"Saved Secuer clustering results: {output_file}")

            # Optional: Generate spatial scatter plot
            if args.plot:
                plot_spatial_scatter(banksy_matrix2, args.output_dir, args.output_prefix, "secuer_cluster")

        cleanup_intermediate_files(args.output_dir, args.output_prefix)

        if args.run_clustering in ["leiden", "both"]:
            n_pcs = args.pca_dims[0]
            sc.pp.neighbors(banksy_matrix2, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=n_pcs)
            sc.tl.leiden(banksy_matrix2, resolution=args.leiden_resolution)

            leiden_dir = os.path.join(args.output_dir, "leiden_clusters")
            os.makedirs(leiden_dir, exist_ok=True)

            for sample_id in banksy_matrix2.obs[args.sample_id_column].unique():
                # Filter the obs dataframe by the sample ID
                subset_obs = banksy_matrix2.obs[banksy_matrix2.obs[args.sample_id_column] == sample_id].copy()
            
                # Rename the 'leiden' column to 'group'
                subset_obs.rename(columns={'leiden': 'group'}, inplace=True)
            
                #Define the output file name
                output_file = os.path.join(leiden_dir, f"{sample_id}_leiden_obs.csv")
            
                # Save the subset to a CSV file
                subset_obs.to_csv(output_file, index=True)
                print(f"Saved Leiden clustering results: {output_file}")

            if args.plot:
                plot_spatial_scatter(banksy_matrix2, args.output_dir, args.output_prefix, "leiden")

        
        # Save final annotated matrix
        final_output = os.path.join(args.output_dir, f"{args.output_prefix}_final_banksy_matrix.h5ad")
        banksy_matrix2.write(final_output)
        print(f"Saved final Banksy matrix to {final_output}")
        


    except Exception as e:
        print(f"An error occurred: {e}")
        raise

def plot_spatial_scatter(adata, output_dir, prefix, cluster):
    spatial_coords = adata.obsm['spatial']
    group = adata.obs.get(cluster, None)
    if group is not None:
        group = group.astype(int)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    scatter = ax.scatter(spatial_coords[:, 0], spatial_coords[:, 1], c=group, cmap='tab20', s=5, alpha=0.8)
    plt.colorbar(scatter, ax=ax, label="Clusters")
    plt.xlabel("Spatial Coordinate 1")
    plt.ylabel("Spatial Coordinate 2")
    plt.title(f"Spatial Scatter Plot ({cluster})")
    plt.tight_layout()
    plot_path = os.path.join(output_dir, f"{prefix}_{cluster}_spatial_scatterplot.png")
    plt.savefig(plot_path)
    print(f"Saved {cluster} spatial scatter plot to {plot_path}")

def cleanup_intermediate_files(output_dir, output_prefix):
    """
    Remove intermediate files if the pipeline completes successfully.
    """
    intermediate_files = [
        os.path.join(output_dir, f"{output_prefix}_normalized_filtered.h5ad"),
        os.path.join(output_dir, f"{output_prefix}_harmony.pkl"),
        os.path.join(output_dir, f"{output_prefix}_banksy_matrix.h5ad")
    ]
    
    for file_path in intermediate_files:
        if os.path.exists(file_path):
            os.remove(file_path)
            print(f"Deleted intermediate file: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Banksy Harmonized Pipeline")
    parser.add_argument("--input_merged_anndata", required=True, help="Path to input merged AnnData file.")
    parser.add_argument("--output_dir", required=True, help="Output directory for results.")
    parser.add_argument("--output_prefix", default="banksy", help="Prefix for all output files.")
    parser.add_argument("--n_top_genes", type=int, default=2000, help="Number of top highly variable genes.")
    parser.add_argument("--k_geom", type=int, default=15, help="K geometry for Banksy initialization.")
    parser.add_argument("--max_m", type=int, default=1, help="Azumithal transform up to kth order.")
    parser.add_argument("--nbr_weight_decay", default="scaled_gaussian", choices=["scaled_gaussian", "reciprocal", "uniform", "ranked"], help="Neighbor weight decay method.")
    parser.add_argument("--pca_dims", type=int, nargs="+", default=[20], help="Dimensionality reduction for PCA.")
    parser.add_argument("--lambda_list", type=float, nargs="+", default=[0.8], help="List of lambda parameters for Banksy.")
    parser.add_argument("--harmony_batch_key", default="dataset", help="Key for Harmony batch correction.")
    parser.add_argument("--run_clustering", choices=["leiden", "secuer", "both"], default="both", help="Run Leiden, Secuer clustering, or both.")
    parser.add_argument("--leiden_resolution", type=float, default=0.5, help="Resolution for Leiden clustering.")
    parser.add_argument("--secuer_resolution", type=float, default=1.0, help="Resolution for Secuer clustering.")
    parser.add_argument("--sample_id_column", default="dataset", help="Column name in adata.obs representing sample ID.")
    parser.add_argument("--plot", action="store_true", help="Generate spatial scatter plot.")
    args = parser.parse_args()
    
    main(args)


