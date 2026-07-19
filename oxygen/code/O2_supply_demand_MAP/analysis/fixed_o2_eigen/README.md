# Fixed-O2 eigen-attractor analysis

## File registry

### `analyze_fixo2_eigen_attractor_embeddings.R`

- Reads only materialized fixed-O2 eigen wide feature tables.
- Produces z-score metadata, PCA/UMAP/t-SNE coordinate tables, PCA variance,
  k-means cluster assignments, and silhouette tables.
- Writes no figures or reports and never invokes the simulator.
- Called by `runner/fixed_o2_eigen/run_fixo2_eigen_attractor_pipeline.R`.
