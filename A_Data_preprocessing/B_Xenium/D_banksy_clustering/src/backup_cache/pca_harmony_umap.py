"""
This module performs PCA, Harmony batch correction, and optional UMAP
embedding on the banksy matrices.

Author: Adapted by ChatGPT, June 2025
"""

from sklearn.decomposition import PCA
import scanpy as sc
import harmonypy as hm
import umap
import numpy as np
from scipy.sparse import issparse
from banksy_utils.pca import plot_remaining_variance
from typing import List


def pca_harmony_umap(
    banksy_dict: dict,
    pca_dims: List[int] = [20],
    plt_remaining_var: bool = True,
    add_umap: bool = True,
    harmony_batch_key: str = "dataset",
    harmony_seed: int = 42,
) -> None:
    """
    Performs PCA and Harmony on banksy_dict entries, and optionally adds UMAP.

    Args:
        banksy_dict (dict): Dictionary of Banksy matrices.
        pca_dims (List[int]): List of PCA dimensions to reduce to.
        plt_remaining_var (bool): Whether to plot explained variance.
        add_umap (bool): Whether to compute and store UMAP embeddings.
        harmony_batch_key (str): Key in .obs to use for batch correction.
        harmony_seed (int): Seed for Harmony.
    """

    for decay in banksy_dict:
        for lambda_val in banksy_dict[decay]:
            if isinstance(lambda_val, str):
                continue  # skip non-lambda entries

            print(f"\n>>> Processing {decay} | lambda = {lambda_val} <<<\n")
            adata = banksy_dict[decay][lambda_val]["adata"]

            # Extract matrix
            X = adata.X.todense() if issparse(adata.X) else adata.X
            X[np.isnan(X)] = 0

            for pca_dim in pca_dims:
                print(f"  → Running PCA with {pca_dim} components")
                pca = PCA(n_components=pca_dim)
                X_pca = pca.fit_transform(X)
                adata.obsm[f"X_pca"] = X_pca

                if plt_remaining_var:
                    plot_remaining_variance(pca, title=f"{decay}, λ={lambda_val}")

                print(f"  → Running Harmony on 'X_pca' with key '{harmony_batch_key}'")
                ho = hm.run_harmony(X_pca, adata.obs, vars_use=[harmony_batch_key], random_state=harmony_seed)
                X_harmony = ho.Z_corr.T
                adata.obsm["X_pca_harmony"] = X_harmony

                if add_umap:
                    print(f"  → Running UMAP on Harmony output")
                    reducer = umap.UMAP(random_state=42)
                    umap_result = reducer.fit_transform(X_harmony)
                    adata.obsm["X_pca_harmony_umap"] = umap_result

                print(f"✓ Finished processing: {decay} λ={lambda_val}")