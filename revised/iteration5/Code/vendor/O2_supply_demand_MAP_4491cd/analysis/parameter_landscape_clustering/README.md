# Parameter-landscape analysis

This directory consumes materialized simulation tables and produces reductions,
clusters, contribution statistics, and analysis manifests. It does not run the
model or draw graphics.

| File | Responsibility |
|---|---|
| `analyze_parameter_landscape.R` | Canonical clustering/reduction analysis entrypoint. |
| `parameter_contribution_analysis.R` | Canonical mode and dominant-ploidy contribution statistics. |
| `parameter_landscape_analysis_utils.R` | Analysis-only transforms, scaling, clustering, summaries, and compatibility adapters. |
| `parameter_landscape_utils.R` | Deprecated analysis utility loader. |
| `clustering_analysis.R`, `full_data_in_vivo_clustring.R` | Deprecated clustering entrypoints forwarding to the canonical analysis. |
| `mode_parameter_contribution_analysis.R`, `dominant_ploidy_parameter_contribution_analysis.R` | Deprecated contribution entrypoints. |
| `clustering_report.R`, `mode_parameter_contribution_report.R`, `dominant_ploidy_parameter_contribution_report.R` | Deprecated report entrypoints forwarding to `report/parameter_landscape/`. |
| `clustering_runner.R`, `mode_parameter_contribution_runner.R`, `dominant_ploidy_parameter_contribution_runner.R` | Deprecated orchestration entrypoints forwarding to `runner/parameter_landscape/`. |
