# Multi-warmup analysis

All canonical files in this directory are table/statistics/embedding producers.
They do not draw figures, assemble HTML, invoke numerical simulators, or create
downstream `figures/` and `report/` directories.

## Canonical files

### `build_multi_warmup_seed_plan_tables.R`

- Reads separate-fit seed summaries and parameter tables.
- Produces paired/single-side 18D profiles, UMAP coordinate tables, cluster
  summaries, representatives, and the warm-up manifest.

### `build_multi_warmup_landscape_tables.R`

- Consumes materialized parameter-landscape seed tables.
- Produces PCA/UMAP/t-SNE coordinates, cluster/subcluster tables, seed groups,
  pair manifests, seed-space task tables, and validation comparisons.

### `collect_multi_warmup_tables.R`

- Collects completed joint-fit summaries and parameter ratios.
- Produces objective, DEoptim-iteration, basin, and initial-to-final distance tables.

### `analyze_warm_up_joint_results.R`

- Produces joint initial/best in-vivo parameter tables, embeddings, curve-class
  summary tables, association tests, and the joint master table.
- Dense-grid curve execution is deliberately owned by the runner.

## Deprecated compatibility entrypoints

- `build_multi_warmup_seed_plan.R` delegates to the seed-plan runner.
- `build_multi_warmup_pairs_from_landscape_subclusters.R` delegates to the
  landscape-pair runner.
- `collect_multi_warmup_results.R` delegates to the collection runner.
- `multi_warmup_results_report.R` delegates to the visualization/report runner.
