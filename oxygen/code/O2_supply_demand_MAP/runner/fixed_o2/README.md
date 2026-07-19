# Fixed-O2 compatibility runner

This directory exists only to preserve the two historical all-in-one fixed-O2
CLIs while callers migrate to staged entrypoints.

## Files

- `fixed_o2_pipeline_loader.R`: loads the numerical, analysis, visualization,
  report, and legacy orchestration modules into one compatibility environment.
  It parses the canonical simulator to import required function definitions
  without executing its top-level model load or CLI `main()`.
- `fixed_o2_legacy_pipeline_functions.R`: the old cross-stage orchestration
  functions and `fixo2_main()`. It preserves historical arguments, output names,
  and execution order. It is intentionally isolated here because it calls more
  than one functional stage.

The compatibility wrappers are:

- `analysis/fixed_o2/FixO2_invivo.R`, defaulting to
  `oxygen/results/analysis/FixO2_invivo_500seed`
- `analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo.R`, defaulting
  to `oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed`

No new workflow should source this directory for a single-stage task.
