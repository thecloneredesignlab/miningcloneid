# Parameter-landscape simulation

This directory materializes fitted-parameter numerical features. It does not
classify clusters, draw figures, or assemble reports.

| File | Responsibility |
|---|---|
| `generate_parameter_landscape_simulation_tables.R` | CLI producer for all/in-vivo/in-vitro materialized feature tables and simulation manifests. |
| `parameter_landscape_invivo_feature_simulation.R` | Numerical in-vivo probes derived from each fitted parameter set. |
| `parameter_landscape_simulation_utils.R` | Simulation-only feature builders and fitted-result adapters. |

Use `../../runner/parameter_landscape/run_parameter_landscape.R`; direct
outputs are consumed by `analysis/parameter_landscape_clustering/`.
