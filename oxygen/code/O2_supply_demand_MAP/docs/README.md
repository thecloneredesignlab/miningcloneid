# O2 supply-demand MAP workflow guide

This guide describes the current fitting and post-fit architecture under
`oxygen/code/O2_supply_demand_MAP/`. Detailed fitting options remain in
`oxygen/README.md`; concrete purpose, inputs, outputs, defined functions, and
direct tests for every source file are generated in `CODE_FILE_REGISTRY.md` and
`code_file_registry.tsv`.

## Protected core

The complete contents of these directories are immutable:

- `model/`
- `optimizer/`

Do not add, remove, rename, or edit source, cache, metadata, or documentation in
either directory. Verify their tracked and complete-tree manifests with:

```bash
Rscript oxygen/tests/check_immutable_o2_core.R
Rscript oxygen/tests/check_immutable_o2_core.R --full
```

The protected fitting dispatcher relies on stable backends under `util/`; those
paths and fitting-facing interfaces are documented in `../util/README.md`.

## Architecture

The canonical dependency direction is:

```text
fit outputs / best parameters
              |
              v
          simulation
              |
              v
           analysis
              |
              v
              vis
              |
              v
            report
```

Each executable layer may use shared libraries from `util/`. Only `runner/`
and `hpc/` may sequence multiple executable layers.

| Layer | Contract |
|---|---|
| `simulation/` | Evaluate completed fitting results/best parameters and materialize numerical O2, population, ploidy, CIN, response, perturbation, or feature tables. |
| `analysis/` | Consume fit summaries or simulation artifacts and write statistics, diagnostics, clustering, ranking, and plot-ready tables. |
| `vis/` | Consume simulation/analysis tables and write figures plus visualization manifests. |
| `report/` | Assemble existing tables and figures into HTML/PDF presentation artifacts. |
| `util/` | Provide source-safe shared functions and stable fitting backends. |
| `runner/` | Sequence local fitting or staged post-fit entrypoints. |
| `hpc/` | Configure and submit Slurm jobs that call canonical runner/layer entrypoints. |

Layer-specific folder maps and restrictions are in `../analysis/README.md`,
`../simulation/README.md`, `../vis/README.md`, `../report/README.md`,
`../util/README.md`, `../runner/README.md`, and `../hpc/README.md`.

## Fitting entrypoints

The unified local fitting launcher is:

```bash
bash oxygen/code/O2_supply_demand_MAP/runner/run_o2_fit.sh \
  --fitting_mode=invivo \
  --config_path=oxygen/config/O2_supply_demand.yaml \
  --invivo_total_seeds=1 \
  --n_cores=1
```

Supported fitting scopes are `invivo`, `invitro`, and `joint`. Joint execution
supports direct fitting, fitting from selected single-fit anchors, and
multi-warmup orchestration. Stable low-level shell paths remain under
`runner/run_fit_model_O2_supply_demand_MAP.sh` and
`runner/run_fit_joint_model_O2_supply_demand_MAP.sh` for existing callers.

For production arrays use:

```bash
bash oxygen/code/O2_supply_demand_MAP/hpc/submit/submit_o2_fit.sh \
  --fitting_mode=invivo \
  --config_path=/absolute/path/to/O2_supply_demand.yaml \
  --out_root=/absolute/path/to/results
```

Use `--dry_run=TRUE` on either orchestrator to inspect commands without running
local fits or submitting Slurm jobs.

## Standard completed-seed workflow

After fitting one seed, run the canonical post-fit orchestrator:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/run_postfit_pipeline.R \
  --fit_dir=/absolute/path/to/run/seed1 \
  --scope=invivo
```

Valid scopes are `invivo`, `invitro`, and `joint`. The runner invokes enabled
stages in simulation → analysis → visualization → report order. Per-layer
enable flags and explicit output-directory options support partial rerenders
without changing layer ownership.

Default seed-relative output roots are:

```text
seed1/
  simulation/
  analysis/
  viz/
  report/
```

Visualization and report consumers must continue to work when fit RDS objects
and best-parameter files are absent from a copied fixture, provided their
documented materialized tables/manifests are present.

## Specialized workflows

| Workflow | Canonical runner or producer |
|---|---|
| Cross-seed extra results | `runner/fit_results/run_extra_results.R` |
| Paired joint-sigma comparison | `runner/fit_results/run_joint_sigma_soft_coupled_paired_seeds.R` |
| Sigma-burden comparison | `runner/fit_results/run_sigma_burden_extra_results.R` |
| Long-ploidy seed selection | `runner/fit_results/select_invivo_best_long_ploidy_seed.R` |
| Fixed-O2 sweep/analysis | `runner/fixed_o2/` |
| Fixed-O2 eigen attractors | `runner/fixed_o2_eigen/run_fixo2_eigen_attractor_pipeline.R` |
| Parameter landscape | `runner/parameter_landscape/run_parameter_landscape.R` |
| Multi-warmup landscape/results | `runner/multi_warmup/` |
| Warm-start table construction | `runner/warm_start/` |
| Live-effective `p_ms` materialization | `simulation/invivo/cin/generate_live_effective_pms_outputs.R` |
| Live-effective `p_ms` cross-run comparison | `runner/profile_likelihood/run_live_effective_pms_comparison.R` |
| Profile-likelihood result collection | `analysis/profile_likelihood/collect_profile_likelihood_results.R` |

Feature-specific READMEs in these folders document required manifests and
output tables. Old public paths listed in `path_migration_table.tsv` are retained
only as labelled compatibility wrappers/orchestrators; new callers should use
the canonical paths above.

## In-vitro helper migration

The former `oxygen/code/in-vitro-utils/` directory is no longer a runtime or
compatibility path. Its reusable fitting helpers are loaded from:

```text
util/o2_supply_demand_map_invitro_utils.R
```

Pure in-vitro plot constructors live under `vis/invitro/`. The historical 58
`ivt_*` functions and their formal arguments are protected by migration tests.

## Output and manifest rules

- A simulation producer records its output tables and provenance in its own
  simulation directory.
- Analysis validates required upstream data and writes analysis-owned tables;
  it never backfills missing simulation by launching a producer.
- Visualization reads only materialized tables/manifests and writes figures
  plus a visualization manifest.
- Reports omit unavailable optional presentation assets or fail with a clear
  missing-required-input error; they do not regenerate scientific data or
  analytical figures.
- Compatibility wrappers preserve established CLI/output contracts while
  forwarding work to canonical staged entrypoints.

## Regression and maintenance

Stage evidence is recorded in:

- `stage0_regression_baseline.md`
- `stage1_regression.md`
- `stage2_regression.md`
- `stage3_regression.md`
- `stage4_regression.md`
- `stage5_regression.md`

Run the complete organization regression after structural changes:

```bash
Rscript oxygen/tests/run_o2_reorganization_regression.R
```

The aggregate gate parses every R file, checks shell/Slurm and Python syntax,
runs the layer-boundary and immutable-core checks, runs the full unit suite, and
checks the Git diff for whitespace errors.

Regenerate the per-file registry after the source tree is stable:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/documentation/generate_code_file_registry.R \
  --workflow_root=oxygen/code/O2_supply_demand_MAP
```

The registry generation test verifies that every supported source file appears
once and that paths are current. `path_migration_table.tsv` remains the source
of truth for old-to-canonical entrypoint mappings.
