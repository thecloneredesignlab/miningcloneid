# Fixed-O2 simulation

This directory owns numerical production under fixed oxygen. It may load fitted
parameters and the model, construct the fixed-O2 generator, solve or simulate
trajectories, and write numerical tables. It does not draw figures or assemble
reports.

## Files

- `run_fixed_o2_simulation.R`: canonical command-line simulator. It resolves
  fitted parameters, initializes 2N/4N states, runs requested replicates, and
  writes population, state, rate, parameter-audit, and metadata tables. Batch
  mode expands seed/O2/initial-ploidy/replicate combinations and skips complete
  tasks unless forced.
- `fixed_o2_simulation_utils.R`: compatibility loader for shared fixed-O2 path,
  argument, parameter, and model-loading utilities in `util/`.
- `fixed_o2_numerical_producers.R`: numerical functions extracted from the
  former `analysis/**/FixO2_invivo.R` monoliths. It owns attractor/eigen/expm/
  Euler calculations, counterfactual trajectories, missing-task generation,
  analytical trajectories, and readers that convert materialized state tables
  into numerical trajectory metrics. This is a function module, not a CLI.

## Output contract

Numerical producers must write tables before downstream stages run. The staged
analysis entry accepts, at minimum:

- `attractors/tables/fixed_o2_attractors_by_seed.tsv`
- `attractors/tables/parameter_values_long.tsv` when parameter correlations are
  requested
- `counterfactual_trajectories/tables/fixed_o2_counterfactual_summary_by_seed.tsv`
- `counterfactual_trajectories/tables/fixed_o2_counterfactual_trajectories.tsv`
- validation and analytical-agreement tables under their respective `tables/`
  directories

Shared helpers belong in `util/`; cross-stage compatibility orchestration lives
in `runner/fixed_o2/`.
