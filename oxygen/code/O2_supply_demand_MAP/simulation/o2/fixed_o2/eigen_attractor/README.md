# Fixed-O2 eigen-attractor simulation

This directory owns numerical feature production from fitted or optimizer-initial
parameter vectors. It may load the fixed-O2 simulator and protected model, but it
does not perform embeddings, clustering, plotting, or report assembly.

## File registry

### `build_fixo2_eigen_attractor_features.R`

- Reads selected seed directories, `fit_config.rds` or embedded fit-result
  configuration, parameter-landscape best/initial parameter tables, and an O2 grid.
- Evaluates the dominant eigen-attractor and spectral gap at every requested O2.
- Writes long/wide feature tables, status summaries, HPC task tables, and merged
  task-row products under `FixO2EigenAttractors/`.
- Preserves configuration precedence: `fit_config.rds`, embedded fit-result cfg,
  then the fixed-O2 fallback configuration.
- Called by the fixed-O2 eigen runner and the worker in `hpc/fixo2_eigen_attractor/`.
