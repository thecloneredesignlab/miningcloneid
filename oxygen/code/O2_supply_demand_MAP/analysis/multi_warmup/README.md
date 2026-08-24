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
- For the active joint workflow, produces a pooled t-SNE, clusters only in-vivo
  best points at the primary level, selects each cluster's objective-minimum
  seed, and pairs those representatives with one globally best in-vitro seed.
- Does not run second-level clustering or curve filtering.

### `build_joint_primary_cluster_pairs.R`

- Canonical command-line entry for the only supported joint pair-selection
  workflow.

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
  landscape-pair runner for archived callers; unified fitting does not call it.
- `collect_multi_warmup_results.R` delegates to the collection runner.
- `multi_warmup_results_report.R` delegates to the visualization/report runner.
