# Fixed-O2 analysis

This directory owns table-to-table analysis only. It does not load fitted
objects, source the model, estimate trajectories, draw figures, or assemble
HTML.

## Files

- `run_fixed_o2_analysis.R`: staged CLI. It requires `--simulation_dir`, reads
  materialized numerical tables, and writes derived statistics under
  `--analysis_dir` (default: `SIMULATION_DIR/analysis`). `--parts` supports
  `attractors`, `counterfactual`, and `agreement`.
- `fixed_o2_analysis_utils.R`: pure classification and analysis functions:
  reference-O2 mode assignment, regime summaries/tests, parameter correlations,
  spectral-gap reliability tables, counterfactual summaries, replicate
  aggregation, and analytical-versus-simulation agreement metrics.
- `FixO2_invivo.R`: deprecated 53-line compatibility wrapper. It preserves the
  historical CLI/default output location by delegating to
  `runner/fixed_o2/fixed_o2_pipeline_loader.R`. New jobs should call the staged
  simulation, analysis, visualization, and report entrypoints explicitly.

## Main outputs

The attractor stage writes mode tables, regime summaries/tests, parameter
correlations, and spectral-gap reliability tables under `attractors/tables/`.
The counterfactual stage writes regime summaries/tests and parameter
correlations under `counterfactual_trajectories/tables/`. The agreement stage
writes residual-augmented matched data and Bland-Altman summaries under
`analytical_simulation_agreement/tables/`.

All inputs are ordinary TSV files, so this stage can run from any working
directory without model compilation or access to optimizer state.
