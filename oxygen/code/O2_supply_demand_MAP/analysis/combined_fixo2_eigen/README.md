# Combined FixO2 eigen-attractor analysis

This directory contains table-only consumers. It does not simulate trajectories, draw figures, or render reports.

- `calculate_regression_curve_average_slope.R` reads the dense-grid regression curves and optional by-seed annotations, computes one average-slope row per seed, and writes `fixed_o2_ploidy_regression_curve_average_slope_by_seed.tsv`.
- `prepare_fixo2_eigen_curve_class_tables.R` reads materialized FixO2 eigen-attractor coordinate CSVs, the dense-grid curve-class table, and the slope table. It writes annotated coordinate CSVs, legacy best-point/count tables, and `pooled_embedding_curve_class_analysis_manifest.tsv` for visualization.
- `MIGRATION.md` records the legacy-to-canonical entrypoint mapping and preserved output contract.

Run the complete workflow through `runner/combined_fixo2_eigen/run_fixo2_eigen_curve_class_pipeline.R`.
