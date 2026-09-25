# Fixed-O2 eigen-attractor runners

## File registry

### `run_fixo2_eigen_attractor_pipeline.R`

- Full cross-stage CLI for `features`, `reductions`/`analysis`, `figures`, and
  `report` parts.
- Sequences simulation, analysis, visualization, and report modules through
  materialized files while preserving the historical CLI defaults.
- The historical `05_.../fixo2_eigen_attractor_clustering_runner.R` delegates here.
