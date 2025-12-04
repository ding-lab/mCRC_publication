for sample_id in adata.obs["dataset"].unique():
    subset_obs = adata.obs[adata.obs["dataset"] == sample_id].copy()
    subset_obs.rename(columns={'leiden': 'group'}, inplace=True)
    output_file = os.path.join(leiden_dir, f"{sample_id}_leiden_obs.csv")
    subset_obs.to_csv(output_file, index=True)
    print(f"Saved Leiden clustering results: {output_file}")