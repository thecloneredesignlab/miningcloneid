# Simulation layer

`simulation/` uses completed fitting results and best parameters to materialize
concrete biological and numerical quantities. Examples include O2 trajectories,
population sizes, ploidy states, CIN/missegregation quantities, functional
responses, perturbation trajectories, and fitted-parameter predictions needed
by downstream analysis.

Simulation writes data and manifests. It does not rank or statistically compare
runs, draw figures, or assemble reports.

## Functional organization

| Folder | Numerical products |
|---|---|
| `invivo/` | Standard in-vivo post-fit materialization. Domain modules under `o2/`, `population/`, `ploidy/`, `cin/`, and `functional_response/` keep biological responsibilities separate. |
| `invitro/` | Standard in-vitro post-fit materialization. Domain modules under `o2/`, `population/`, `ploidy/`, and `cin/` provide the corresponding tables. |
| `o2/fixed_o2/` | Fixed-O2 sweeps, attractors, reliability products, and eigen-attractor numerical features. |
| `fit_results/` | Predictions and long-ploidy metrics required before cross-seed fit-result analysis. |
| `parameter_landscape/` | Per-seed numerical features used by parameter-landscape analyses. |
| `perturbation/` | Mixed-ploidy and factorial perturbation simulations. |
| `process_fingerprints/` | Numerical process trajectories and event inputs used to build process/ploidy fingerprints. `process_fingerprint_simulation_utils.R` handles materialized simulation contracts; `process_fingerprint_simulation_legacy_utils.R` retains simulation-only fitted-result readers and numerical feature production for the compatibility surface. Shared parsing/transformation/distance/clustering and ploidy-regime primitives come from `../util/o2_supply_demand_map_process_fingerprint_utils.R` and `../util/o2_supply_demand_map_ploidy_regime_utils.R`. |

The top-level `fix_o2_simulation*.R` and `simulate_invivo_*.R` files are
deprecated compatibility entrypoints. Their headers identify the canonical
staged producers they forward to.

## Standard entrypoints

Materialize a completed in-vivo seed:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/simulation/invivo/generate_invivo_simulation_outputs.R \
  --fit_dir=/absolute/path/to/seed
```

Materialize a completed in-vitro seed:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/simulation/invitro/run_invitro_simulation_outputs.R \
  --fit_dir=/absolute/path/to/seed
```

For simulation followed by analysis, visualization, and reporting, call
`../runner/run_postfit_pipeline.R` instead of manually chaining layers.

## Boundary checks and file details

Simulation may read protected fitting/model interfaces only where needed to
evaluate an already fitted parameter set. It must not source downstream layers
or write `figures/` or `report/` products. Verify this with:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
```

Concrete purpose, inputs, outputs, functions, and direct tests for every source
file are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`.
