# Visualization layer

`vis/` is a pure consumer layer. It reads materialized simulation or analysis
tables and writes figures plus visualization manifests. It does not read
`fit_result.rds`, load best-parameter files, source model/optimizer code, or
reconstruct numerical trajectories.

## Visualization folders

| Folder | Figures produced from existing tables |
|---|---|
| `invivo/` | Reusable in-vivo plot constructors used by the standard in-vivo consumer. |
| `invitro/` | In-vitro plot and diagnostic constructors. |
| `joint/` | In-vivo/in-vitro comparisons and joint-parameter ratio figures. |
| `joint_soft_coupling_stability/` | Within-/between-pair ClassA/B/C stability, sensitivity, direction, and biological-process figures. |
| `joint_ploidy_coupling_association/` | CatA/B/C/U trajectories, classifier diagnostics, Cat×parameter, and Cat×ratio-class association figures. |
| `joint_fixed_o2_ploidy_classification/` | Pair-specific and all-pair regression-smoothed fixed-O2 steady-state ploidy curves, faceted by curve class and colored only by temporal Cat. |
| `fit_results/` | Cross-seed objective, boundary, paired-sigma, and sigma-burden figures. |
| `fixed_o2/` | Fixed-O2 sweep and reliability figures. |
| `fixed_o2_eigen/` | Fixed-O2 eigen-attractor embedding and clustering figures. |
| `combined_fixo2_eigen/` | Fixed-O2 eigen-attractor embeddings annotated from analysis-prepared curve-class and slope tables. |
| `combined_parameter_landscape/` | Pooled parameter-landscape embeddings annotated from existing classification tables. |
| `dense_grid_monotonicity/` | Dense-grid monotonicity and curve-class figures from analysis tables. |
| `interactions/` | Factorial-interaction figures. |
| `multi_warmup/` | Warm-up landscape, seed-plan, collection, and report figures. |
| `parameter_landscape/` | Parameter-landscape embeddings and contribution figures. |
| `perturbation/` | Mixed-ploidy perturbation figures. |
| `process_fingerprints/` | Process fingerprints, ploidy regimes, medium-O2 windows, and O2-ploidy coupling figures. |
| `profile_likelihood/` | Live-cell effective `p_ms` comparison figures. |

The historical top-level consumers
`viz_invivo_model_O2_supply_demand_MAP_results.R` and
`viz_invitro_model_O2_supply_demand_MAP_results.R` now require already
materialized tables. Numerical production is handled by
`../simulation/invivo/generate_invivo_simulation_outputs.R` and
`../simulation/invitro/run_invitro_simulation_outputs.R`.

Shared plot-scale constructors used by more than one visualization consumer
live in `o2_supply_demand_map_common_plot_utils.R`; domain-specific plot
constructors remain in their functional subfolders.

## Recommended execution

For standard completed seeds, use `../runner/run_postfit_pipeline.R`; it invokes
simulation, analysis, visualization, and report entrypoints in dependency
order. For a visualization-only rerender, call the applicable script only after
its upstream manifest and required tables exist.

The global architecture gate enforces the consume-only boundary:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
```

Concrete purpose, required inputs, outputs, functions, and direct tests for
every visualization file are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`.
