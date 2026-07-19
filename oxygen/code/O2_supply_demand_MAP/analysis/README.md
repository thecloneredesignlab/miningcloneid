# Analysis layer

`analysis/` consumes fitted summaries and already materialized simulation
tables. It produces statistical summaries, diagnostics, rankings, clustering
results, and plot-ready tables. It does not run the numerical model, launch a
simulation, or draw figures.

The required dependency direction is:

```text
fit outputs -> simulation -> analysis -> vis -> report
```

Only `runner/` and `hpc/` may sequence more than one executable layer. Shared
functions used by multiple workflows belong in `../util/`.

## Feature folders

| Folder | Analysis responsibility |
|---|---|
| `fit_diagnostics/` | In-vitro and joint-parameter diagnostics from materialized post-fit data. |
| `fit_results/` | Cross-seed ranking, boundary summaries, paired-run comparisons, and seed selection. Historical executable paths in this folder are compatibility orchestrators where noted in their headers. |
| `fixed_o2/` | Fixed-O2 summaries and analytical tables from fixed-O2 simulation products. |
| `fixed_o2_eigen/` | Clustering and embedding summaries from materialized fixed-O2 eigen-attractor features. |
| `interactions/` | Factorial-interaction effect tables from perturbation simulations. |
| `multi_warmup/` | Warm-up seed-plan, landscape, pair-selection, and collection tables. |
| `parameter_landscape_clustering/` | Parameter-landscape clustering and contribution tables. |
| `perturbation/` | Mixed-ploidy perturbation comparisons. |
| `process_fingerprints/` | Process fingerprints, ploidy regimes, medium-O2 windows, and O2-ploidy event-coupling analyses. `process_fingerprint_utils.R` and `ploidy_regime_utils.R` are source-compatible analysis loaders for canonical shared modules in `../util/`; they contain no duplicated simulation or plotting implementation. |
| `profile_likelihood/` | Profile-likelihood collection and comparisons of materialized live-cell effective `p_ms` values. |
| `warm_up_joint_fitting_results_extra/` | Joint warm-up fit-result summaries; HPC launchers for this workflow live in `../hpc/warm_up_joint_fitting_results_extra/`. |
| `best_fit_parameter_feature/` | Deprecated path-compatible forwards for the staged best-fit feature workflows. Canonical producers, analyses, visualizations, reports, and runners live in their corresponding layers. |
| `dense-grid_monotonicity_classification/` | Historical compatibility surface for the staged dense-grid workflow. |

## How to run analyses

For the standard post-fit workflow, use the layer orchestrator:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/run_postfit_pipeline.R \
  --fit_dir=/absolute/path/to/seed \
  --scope=invivo
```

Specialized pipelines are under `../runner/fit_results/`,
`../runner/fixed_o2/`, `../runner/fixed_o2_eigen/`,
`../runner/multi_warmup/`, `../runner/parameter_landscape/`, and
`../runner/profile_likelihood/`. Run an analysis script directly only when its
documented upstream simulation tables and manifest already exist.

## Boundary checks and file details

The architecture gate rejects analysis code that sources `model/`,
`simulation/`, or `vis/`, invokes a simulator, or writes analytical figures:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
```

Concrete purpose, inputs, outputs, functions, and direct tests for every source
file are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`. Regenerate them after any structural change
with `../runner/documentation/generate_code_file_registry.R`.
