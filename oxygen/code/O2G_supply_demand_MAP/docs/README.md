# O2 Supply-Demand MAP Workflow README

This directory contains the fitting, visualization, post-processing, automated boundary calibration, and profile likelihood workflows for the `O2G_supply_demand_MAP` model, together with a few shared internal helper scripts.

This README focuses on:

- what each script does
- which scripts are intended as entrypoints
- how each script should be called
- the main arguments that matter in practice
- where outputs are written
- which files are internal dependencies and should not usually be called directly

## Script Overview

### Recommended entrypoint scripts

- `run_fit_model_O2G_supply_demand_MAP.sh`
- `fit_model_O2G_supply_demand_MAP.R`
- `viz_invivo_model_O2G_supply_demand_MAP_results.R`
- `render_fit_report.R`
- `extra_results.R`
- `extra_results_report.R`
- `auto_calibrate_boundary_params.R`
- `profile_likelihood_O2G_supply_demand_MAP.R`
- `collect_profile_likelihood_results.R`
- `estimate_live_effective_pms.R`
- `compare_sigma_burden_live_effective_pms.R`

### HPC submission scripts

- `submit_profile_likelihood_array.sh`
- `submit_profile_likelihood_array.sub`

### Internal dependency scripts that are not normally run directly

- `model_O2G_supply_demand_MAP.R`
- `model_O2G_supply_demand_MAP.cpp`
- `o2g_supply_demand_map_common_semantics.R`
- `o2g_supply_demand_map_fit_invivo_backend.R`
- `o2g_supply_demand_map_fit_invitro_backend.R`
- `o2g_supply_demand_map_shared.R`

### Historical or non-primary files

- `model_O2G_supply_demand_MAP copy.cpp`
  - This is a copy file and is not the standard entrypoint for the active workflow.

## Recommended Workflow Order

The most common workflow is:

1. Run fitting with `run_fit_model_O2G_supply_demand_MAP.sh --fit_invivo` or `fit_model_O2G_supply_demand_MAP.R --fit_invivo --mode=run`
2. Inspect each `seed` directory for `fit_summary.tsv`, `best_params.tsv`, `viz/`, and `reprot/`
3. Run `extra_results.R` on a completed run directory to rank seeds and assess boundary behavior
4. If boundary expansion is needed, run `auto_calibrate_boundary_params.R`
5. If single-parameter profile likelihood analysis is needed, run `profile_likelihood_O2G_supply_demand_MAP.R`
6. After all profile tasks finish, run `collect_profile_likelihood_results.R`
7. If an effective live-cell missegregation rate estimate is needed from an existing run or seed, run `estimate_live_effective_pms.R`
8. If multiple `sigma_burden` analyses need to be compared at the seed level for `p_misseg` and live-cell effective `p_ms`, run `compare_sigma_burden_live_effective_pms.R`

## Common Paths

This README assumes the current working directory is already the repository root:

- `oxygen/code/O2G_supply_demand_MAP/`
- `oxygen/config/O2G_supply_demand.yaml`
- `oxygen/data/O2G_supply_demand/parameter_table_O2.csv`
- `oxygen/data/O2G_supply_demand/parameter_table_O2G.csv`
- `oxygen/results/`

All command examples below therefore use repository-root-relative paths.

When `parameter_table` is left unset in the YAML or CLI:

- `glucose=FALSE` defaults to `oxygen/data/O2G_supply_demand/parameter_table_O2.csv`
- `glucose=TRUE` defaults to `oxygen/data/O2G_supply_demand/parameter_table_O2G.csv`

## Glucose Family Modes

The O2G workflow now has a top-level `glucose` family switch.

- `glucose=FALSE`
  - fully reverts to the legacy O2-only family
  - death uses only `h_O(O2)`
  - proliferation uses the legacy O2-only saturation formula with `k_o`
  - WGD uses the constant `p_wgd`
  - no glucose state or glucose-specific parameters are estimated
- `glucose=TRUE, glucose_dynamic=FALSE`
  - keeps the current O2G family without an explicit `G(t)` state
  - death and chromosome-instability modules use the coupled fallback `h_G := h_O`
  - proliferation still uses the legacy O2-only growth formula with `k_o`
  - WGD uses the oxygen-triggered `p_wgd_max / O2_wgd` parameterization
- `glucose=TRUE, glucose_dynamic=TRUE`
  - enables the full dynamic-glucose O2G family
  - proliferation uses the combined resource-availability term `R(O2,G)`
  - `k_o` is kept only as a compatibility field and is not estimated

The current natural-scale glucose parameter table defaults are:

- `G_S0`: init `100`, bounds `70-100`
- `kappa_G`: init `20`, bounds `1-30`
- `eta_G`: init `1.2`, bounds `1-2`
- `G_c`: init `30`, bounds `10-60`
- `tau_G`: init `0.1`, bounds `0.01-30`

The `k_o` parameter is therefore only estimated when the active family uses the
legacy O2-only growth law, namely when `glucose=FALSE` or when
`glucose=TRUE, glucose_dynamic=FALSE`.

## Harvest-Specific Initial-Size Multipliers

The in vivo workflow supports an optional harvest-level initial-size adjustment
controlled by `harvest_init_multiplier`.

- `harvest_init_multiplier=FALSE`
  - default behavior
  - no harvest-specific initial-size nuisance parameters are estimated
- `harvest_init_multiplier=TRUE`
  - adds one parameter per independent harvest sample:
    - `log_init_mult_<harvest_id>`
  - the natural-scale multiplier is:
    - `init_mult_<harvest_id> = exp(log_init_mult_<harvest_id>)`
  - each scenario starts from:
    - `N0_i = init_total_size * init_mult_<harvest_id>`

The associated config-controlled prior and bounds are:

- `prior_center_log_init_mult`
- `prior_sd_log_init_mult`
- `log_init_mult_lower`
- `log_init_mult_upper`

The current defaults are:

- `harvest_init_multiplier: FALSE`
- `prior_center_log_init_mult: 0.0`
- `prior_sd_log_init_mult: 0.35`
- `log_init_mult_lower: -1.0`
- `log_init_mult_upper: 1.0`

## 1. `run_fit_model_O2G_supply_demand_MAP.sh`

### Purpose

This is the thinnest wrapper around the unified fitter. It does not parse YAML itself and does not implement workflow logic. You must still choose a mode explicitly. For in vivo fitting it forwards all arguments to:

- `fit_model_O2G_supply_demand_MAP.R --fit_invivo --mode=run`

### Recommended use

- start a full local run
- preserve shell-based calling habits
- avoid typing `Rscript ... --fit_invivo --mode=run` manually

### Example call

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invivo \
  --config=oxygen/config/O2G_supply_demand.yaml
```

An example with YAML overrides:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invivo \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --run_prefix=smoke_local_20260403 \
  --seeds_csv=1,2,3 \
  --auto_viz=FALSE \
  --n_cores=1
```

### Common arguments

- `--config=...`
  - Required. Path to the YAML config file.
- `--run_prefix=...`
  - Overrides the run prefix from YAML.
- `--seeds_csv=...`
  - Explicit comma-separated seed list such as `1,2,3`.
- `--seeds_file=...`
  - Reads seeds from a file.
- `--auto_viz=TRUE|FALSE`
  - Controls whether the visualization script runs automatically after each seed fit.

### Main outputs

The output directory is usually:

- `out_root/run_prefix[_timestamp]/`

Common files in that run directory:

- `run_status.log`
- `config.input.yaml`
- `config.resolved.yaml`
- `parameter_table_input.csv`
- `seed1/`
- `seed2/`
- `...`

Each `seed` directory usually contains:

- `fit_summary.tsv`
- `best_params.tsv`
- `fit_status.log`
- `viz_status.log`
- `report_status.log`
- `single_stage_pass_summary.tsv`
- `viz/`
- `reprot/`

## 2. `fit_model_O2G_supply_demand_MAP.R`

### Purpose

This is the main unified fitter entrypoint. It requires exactly one fit selector and then supports the corresponding backend mode:

- `--fit_invivo`
  - Uses the in vivo workflow.
- `--fit_invitro`
  - Uses the in vitro anoxia workflow.

- `--mode=run`
  - Reads YAML config, creates the run directory, and launches a batch of seed fits.
- `--mode=fit_seed`
  - Executes a single low-level seed fit directly.

### Recommended use

- In almost all cases, use `--fit_invivo --mode=run` for the in vivo workflow
- `--mode=fit_seed` is mainly intended for internal automation, debugging, or controlled one-seed reproductions

### Mode 1: `--mode=run`

#### Recommended call

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/fit_model_O2G_supply_demand_MAP.R \
  --fit_invivo \
  --mode=run \
  --config=oxygen/config/O2G_supply_demand.yaml
```

If `--mode` is omitted but `--config=...` is supplied, the script automatically infers `run` mode:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/fit_model_O2G_supply_demand_MAP.R \
  --fit_invivo \
  --config=oxygen/config/O2G_supply_demand.yaml
```

#### What `run` mode does

- reads YAML config
- merges CLI overrides
- resolves `out_root`, `data_dir`, `parameter_table`, `seeds_file`, and `seeds_csv`
- creates the run directory
- writes:
  - `config.input.yaml`
  - `config.resolved.yaml`
  - `parameter_table_input.csv`
- launches one `fit_seed` run per seed
- optionally launches the visualization script after each seed if `auto_viz=TRUE`
- when `auto_viz=TRUE`, automatically builds the HTML/PDF fit report after successful visualization

#### Common override arguments

- `--run_prefix=...`
- `--out_root=...`
- `--data_dir=...`
- `--parameter_table=...`
- `--start_with=ploidy|chr_number`
  - Optional endpoint-observation mode override.
  - `ploidy` preserves the legacy workflow.
  - `chr_number` uses `total_chromosomes` for the terminal endpoint likelihood and uses chromosome-number scaling in the O2 demand term.
- `--seeds_csv=...`
- `--seeds_file=...`
- `--auto_viz=TRUE|FALSE`
- `--viz_report_dt=1`
- `--viz_top_n=6`
- plus many fitter-specific parameters already present in YAML

#### Seed-level fit report

When `auto_viz=TRUE`, each completed `seed` fit now produces:

- `reprot/fit_report.html`
- `reprot/fit_report.pdf`

The fit report is a report-style rendering of the seed result that includes:

- report metadata
- a 9-field fit summary table
- the best-parameter table
- figure sections for:
  - Burden
  - Ploidy
  - Missegregation
  - Oxygen / O2
  - Glucose / G

The `Glucose / G` section is adaptive:

- if dynamic glucose figures are present, it appears after `Oxygen / O2`
- if glucose is not modeled independently for that fit, the section is omitted

### Mode 2: `--mode=fit_seed`

#### Recommended note

This is a low-level interface. It expects a large number of arguments to already be fully resolved. It is not generally intended for routine manual use unless:

- a script-level debugging session is needed
- a single-seed reproduction is needed
- a new automation wrapper is being developed

#### Minimal schematic example

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/fit_model_O2G_supply_demand_MAP.R \
  --fit_invivo \
  --mode=fit_seed \
  --seed=1 \
  --out_dir=/abs/path/to/seed1 \
  --data_dir=/abs/path/to/data_dir \
  --parameter_table=/abs/path/to/parameter_table.csv \
  --n_cores=1 \
  --use_deoptim=TRUE \
  --deoptim_parallel=FALSE \
  --itermax=1 \
  --NP=12 \
  --n_starts=1 \
  --optim_maxit=20 \
  ...
```

#### Important note

`fit_seed` requires many explicit runtime arguments, including:

- optimization settings
- data filtering settings
- O2, death, and treatment settings
- prior settings
- cache settings
- population scale settings
- `seed`

For routine work, prefer:

- `run_fit_model_O2G_supply_demand_MAP.sh --fit_invivo`
- or `fit_model_O2G_supply_demand_MAP.R --fit_invivo --mode=run`

### Main outputs

A completed `seed` directory typically contains:

- `fit_summary.tsv`
  - total objective, data objective, prior, burden/ploidy objectives, raw `-2logL` diagnostics, and related metadata
- `best_params.tsv`
  - best-fit parameter values
- `single_stage_pass_summary.tsv`
  - summary of the internal fitting stages
- `fit_status.log`
  - seed-level fit log
- `fit_config.rds`
  - saved fit configuration used by downstream scripts

#### Note

- `fit_summary.tsv` now records `start_with`.
- When `start_with=chr_number`, the legacy field names are preserved for compatibility:
  - `objective_ploidy`
  - `objective_ploidy_neg2loglik_raw`
- In that mode, those fields no longer represent a ploidy-space endpoint likelihood. They represent the terminal chromosome-number endpoint likelihood evaluated from `total_chromosomes`.
- Burden volume is now independent of ploidy and chromosome number in all modes:
  - `V_pred = sum_n n_n * (1/rho_2N)`

## 3. `viz_invivo_model_O2G_supply_demand_MAP_results.R`

### Purpose

This script generates visualization outputs and prediction tables for an already completed `seed` directory.

### Recommended call

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R \
  --fit_dir=oxygen/results/your_run/seed1
```

A more explicit example:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R \
  --fit_dir=oxygen/results/your_run/seed1 \
  --data_dir=data/InVivoData_Gemcitabine \
  --report_dt=1 \
  --top_n=6 \
  --n_cores=1
```

### Key arguments

- `--fit_dir=...`
  - Strongly recommended. Path to a single seed directory.
- `--data_dir=...`
  - Optional. If omitted, the script tries to reconstruct it from stored config context.
- `--report_dt=...`
  - Output reporting interval.
- `--top_n=...`
  - Number of top chromosome-count states to highlight.
- `--n_cores=...`
  - Number of workers when multiple fit directories are processed.

### Behavior notes

- If `--fit_dir` is omitted, the script tries to locate the latest fit result directory under the results root.
- If `--out_dir` is supplied, it is ignored.
- Outputs are always written to `fit_dir/viz/`.

### Common outputs

Written under `seed_dir/viz/`:

- `burden_trend.pdf`
- `burden_trend_absolute.pdf`
- `burden_live_dead_decomposition.pdf`
- `predict_burden_vs_o2.tsv`
- `predict_burden_vs_o2.pdf`
- `ploidy_heatmap_over_time.pdf`
- `ploidy_top_states_over_time.pdf`
- `ploidy_weighted_mean_over_time.pdf`
- `ploidy_timecourse.tsv`
- `ploidy_weighted_mean_timecourse.tsv`
- `functional_curve_ploidy.tsv`
- `predict_burden_0_1000day.tsv`
- `overview_9panel.pdf`

### Note

- The visualization layer keeps historical file names for compatibility, including names that contain `ploidy`.
- When `start_with=ploidy`, these files keep their original meaning.
- When `start_with=chr_number`, the following outputs still keep their legacy names but switch to chromosome-number semantics:
  - `ploidy_weighted_mean_over_time.pdf`
  - `ploidy_weighted_mean_timecourse.tsv`
  - `predict_ploidy_weighted_mean_*.tsv`
  - `terminal_ploidy_observed_vs_predicted.tsv`
  - `terminal_ploidy_observed_vs_predicted_violin.pdf`
- In `chr_number` mode, these outputs summarize endpoint chromosome-number behavior or chromosome-number likelihood terms while cohort labels remain `2N` and `4N`.

## 4. `extra_results.R`

### Purpose

This script ranks seeds within a completed run, evaluates parameter boundary behavior, and produces summary tables and plots.

### Recommended call

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R \
  --run_dir=oxygen/results/your_run
```

An example with optional arguments:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R \
  --run_dir=oxygen/results/your_run \
  --out_dir=oxygen/results/your_run/extra_results_custom \
  --near_thresh=0.05
```

### Key arguments

- `--run_dir=...`
  - Required. Can be a multi-seed run directory or a single seed directory.
- `--out_dir=...`
  - Optional. Defaults to `run_dir/extra_results/`.
- `--near_thresh=...`
  - Optional. Boundary-nearness threshold. Default is `0.05`.
  - Must lie in `(0, 0.5)`.

### Main outputs

Written by default to `run_dir/extra_results/`:

- `seed_summary.tsv`
  - seed-level ranking summary
  - includes:
    - `objective`
    - `objective_data`
    - `objective_burden`
    - `objective_ploidy`
    - `objective_burden_neg2loglik_raw`
    - `objective_ploidy_neg2loglik_raw`
    - recommendation ranks
- `parameter_boundary_long.tsv`
  - long-form parameter boundary table
- `objective_components_long.tsv`
  - long-form table for the objective violin plot
- `objective_components_violin.pdf`
  - violin plot of objective, burden, and ploidy components across seeds
- `parameter_boundary_forest.pdf`
  - forest plot of boundary proximity
- `parameter_boundary_forest_pred1000_gt44_top3.pdf`
  - boundary forest plot restricted to the top-3 seeds among runs with 2N/4N 1000-day predictions above 44
- `objective_vs_boundary_risk.pdf`
  - objective vs boundary-risk plot
- `extra_results_report.html`
  - self-contained HTML report that embeds the main extra-results plots in this order:
    - `objective_vs_boundary_risk.pdf`
    - `objective_components_violin.pdf`
    - `parameter_boundary_forest.pdf`
    - `parameter_boundary_forest_pred1000_gt44_top3.pdf`

### Typical use cases

- choose the best seed from a run
- identify boundary-sticking parameters
- provide ranking input to the automated boundary-calibration workflow

## 5. `auto_calibrate_boundary_params.R`

### Purpose

This script automates boundary-sticking parameter calibration. It does not reimplement the model. Instead, it repeatedly orchestrates:

- the fitter
- the visualization script
- `extra_results.R`

while maintaining a working natural-scale parameter table.

### Current target parameters

The current script only calibrates these eight parameters, in this exact order:

1. `rho_2N`
2. `sigma_burden`
3. `eta_o2`
4. `mu_hp`
5. `gamma_mu`
6. `k_o_mis`
7. `p_misseg`
8. `p_wgd`

These parameters are intentionally excluded from the current automated boundary-calibration loop:

- `kappa_O`
- `O2_crit`

### Recommended call

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/auto_calibrate_boundary_params.R \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --output_root=oxygen/results/auto_boundary_local
```

### Configurable arguments

- `--parameter_table=...`
- `--max_rounds_per_parameter=3`
- `--seeds_per_round=10`
- `--boundary_expand_fraction=0.10`
- `--boundary_stick_lower_threshold=0.05`
- `--boundary_stick_upper_threshold=0.95`
- `--objective_threshold_total=9`
- `--objective_threshold_burden=1.5`
- `--objective_threshold_ploidy=7.5`
- `--day1000_min_burden_threshold=2`

### Important behavior

- Only one target parameter is allowed to remain `estimate=TRUE` in each round.
- All non-target parameters are locked.
- A full pre-parameter snapshot is saved before each parameter starts.
- If the parameter ultimately fails, rollback is applied to restore the pre-parameter working table.
- Only successfully completed parameters are allowed to propagate updated values and bounds downstream.

### Main outputs

Written under `output_root/`:

- `auto_calibration_round_log.tsv`
- `auto_calibration_parameter_summary.tsv`
- `auto_calibration_report.md`
- `parameter_table.updated.csv`
- `parameter_table_snapshots/`
- `rollback_records/`

Parameter-specific and round-specific subdirectories also contain:

- round-specific config files
- round-specific parameter tables
- per-seed outputs
- `extra_results` outputs

## 6. `profile_likelihood_O2G_supply_demand_MAP.R`

### Purpose

This script runs a supplement-style profile likelihood workflow for one parameter at a time.

It is not a simple fixed-grid sweep. Instead, it:

- starts from a baseline best fit
- profiles decreasing and increasing directions separately
- fixes the profiled parameter
- re-optimizes the remaining `estimate=TRUE` parameters
- uses warm starts
- uses transformed-space stepping, adaptive step control, threshold refinement, and boundary refinement to generate a denser profile path

### Recommended local test

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/profile_likelihood_O2G_supply_demand_MAP.R \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --baseline_seed_dir=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2 \
  --profile_bounds_table=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2/parameter_table_input.csv \
  --output_root=oxygen/results/profile_local_test \
  --param_name=lam_min \
  --max_steps_per_direction=5 \
  --seeds_per_step=2 \
  --n_cores=1
```

### Example by parameter index

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/profile_likelihood_O2G_supply_demand_MAP.R \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --baseline_seed_dir=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2 \
  --profile_bounds_table=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2/parameter_table_input.csv \
  --output_root=oxygen/results/profile_local_test \
  --param_index=1 \
  --max_steps_per_direction=20 \
  --seeds_per_step=20 \
  --n_cores=62
```

### Required arguments

- `--config=...`
- `--baseline_seed_dir=...`
- `--profile_bounds_table=...`
- `--output_root=...`
- `--param_index=...`
  - or `--param_name=...`

### Common optional arguments

- `--max_steps_per_direction=20`
- `--seeds_per_step=20`
- `--n_cores=62`
- `--target_delta_objective=0.2`
- `--ci_delta_threshold=1.92`
- `--step_fraction_initial=0.10`
- `--step_fraction_min=1e-6`
- `--step_fraction_max=0.30`
- `--boundary_start_tolerance=1e-8`
- `--max_attempts_per_step=5`
- `--min_interior_points_per_direction=4`
- `--ci_refine_tolerance=0.05`
- `--max_refine_steps_ci=6`
- `--max_refine_steps_boundary=4`
- `--boundary_refine_fraction_min=0.01`
- `--max_step_growth_factor=1.5`
- `--max_step_shrink_factor=0.67`
- `--use_soft_prior_for_profile=TRUE|FALSE`
- `--lambda_prior_for_profile=...`

### Baseline requirements

The current implementation expects an existing baseline seed directory. That directory should contain at least:

- `best_params.tsv`
- `fit_summary.tsv`
- `parameter_table_input.csv`

If the baseline is older and does not include the newer raw `-2logL` fields:

- total-objective profile outputs can still run
- raw `-2logL` delta plots and raw CI summaries may not be fully available

### Meaning of the profile plots

Each parameter directory now includes several plot types:

- `profile_curve.pdf`
  - all-seeds version
  - the center line is the mean across complete seeds at each accepted point
  - error bars represent the seed sample range, not confidence intervals
- `profile_delta_curve.pdf`
  - all-seeds delta-objective version
- `profile_delta_burden_neg2loglik_best.pdf`
- `profile_delta_ploidy_neg2loglik_best.pdf`
  - best-seed path
  - based on raw `-2logL` deltas
  - can overlay profile-CI threshold and crossing information
- `profile_delta_burden_neg2loglik_allseeds.pdf`
- `profile_delta_ploidy_neg2loglik_allseeds.pdf`
  - all-seeds version
  - center line is the mean
  - error bars are the seed sample range
  - CI interpretation is based on raw `-2logL` envelope threshold crossings, not on the pointwise bars themselves

### Main outputs

Each parameter usually gets its own directory, such as:

- `01_lam_min/`
- `02_lam_max/`
- `...`

Common files in each parameter directory:

- `profile_path_decreasing.tsv`
- `profile_path_increasing.tsv`
- `profile_path_combined.tsv`
- `profile_seed_results.tsv`
- `profile_point_summary.tsv`
- `profile_likelihood.tsv`
- `direction_summary.tsv`
- `parameter_relation_long.tsv`
- multiple profile PDFs

### Runtime cost

With `seeds_per_step=20` and `max_steps_per_direction=20`, a single parameter can be expensive. A practical pattern is:

- for a local smoke test:
  - `--seeds_per_step=1`
  - `--max_steps_per_direction=1`
  - `--n_cores=1`
- for production profile runs:
  - use the larger defaults on HPC

## 7. `collect_profile_likelihood_results.R`

### Purpose

This script collects all parameter-level outputs under a profile likelihood root directory and writes combined summary tables and a report.

### Recommended call

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/collect_profile_likelihood_results.R \
  --output_root=oxygen/results/profile_local_test
```

### Input expectation

`output_root/` should already contain a set of parameter subdirectories, each usually containing:

- `profile_path_combined.tsv`
- `direction_summary.tsv`
- `profile_point_summary.tsv`
- `parameter_relation_long.tsv`

### Main outputs

Written under `output_root/`:

- `profile_likelihood_all.tsv`
- `profile_direction_summary_all.tsv`
- `profile_point_summary_all.tsv`
- `profile_parameter_summary.tsv`
- `profile_parameter_relations_all.tsv`
- `profile_likelihood_report.md`

### Summary content

The collector summarizes, among other things:

- best-seed profile CI results
- all-seeds envelope raw `-2logL` CI results
- direction-level stopping status
- baseline and best raw `-2logL` summaries

## 8. `submit_profile_likelihood_array.sh`

### Purpose

This is a local wrapper that prepares and submits a SLURM array job for profile likelihood analysis. It is suitable when the shell script itself can be called directly.

### Recommended call

```bash
bash oxygen/code/O2G_supply_demand_MAP/hpc/submit_profile_likelihood_array.sh \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --baseline_seed_dir=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2 \
  --profile_bounds_table=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2/parameter_table_input.csv \
  --output_root=oxygen/results/profile_eq21 \
  --max_steps_per_direction=20 \
  --seeds_per_step=20 \
  --n_cores=62
```

### What it does automatically

- validates the config, baseline, and bounds table paths
- creates the output directory and `slurm_logs/`
- copies baseline metadata into the output directory
- counts `estimate=TRUE` parameters in `profile_bounds_table`
- writes `parameter_targets.tsv`
- writes `submission_manifest.tsv`
- constructs `sbatch --array=...` and submits the `.sub` file

### Key arguments

- `--array_concurrency=...`
- `--mem=128G`
- `--time=7-00:00:00`
- `--job_name=o2sd_profile`
- `--mail_user=...`
- `--mail_type=END`
- `--r_module=R/4.4`

### Current defaults

Current defaults include:

- `max_steps_per_direction = 20`
- `seeds_per_step = 20`
- `n_cores = 62`

## 9. `submit_profile_likelihood_array.sub`

### Purpose

This is the direct `sbatch` entrypoint for the profile likelihood SLURM array workflow.

It is the correct choice when:

- jobs must be submitted through `sbatch file.sub`
- the HPC environment expects the submission entrypoint to be a `.sub` file

### Recommended call

```bash
sbatch oxygen/code/O2G_supply_demand_MAP/hpc/submit_profile_likelihood_array.sub
```

Defaults can still be overridden with `--export`:

```bash
sbatch --export=ALL,OUTPUT_ROOT=oxygen/results/profile_run,MAX_STEPS_PER_DIRECTION=30,SEEDS_PER_STEP=10 \
  oxygen/code/O2G_supply_demand_MAP/hpc/submit_profile_likelihood_array.sub
```

### Important note

This `.sub` file is currently an HPC-specific version. Its default paths are hardcoded for a cluster environment, for example:

- `.../miningcloneid/...`

If the file is reused on another machine, the defaults usually need to be updated, or the paths must be fully overridden through `sbatch --export=...`.

### Current defaults

Current defaults include:

- `#SBATCH --cpus-per-task=62`
- `#SBATCH --array=1-20`
- `DEFAULT_MAX_STEPS_PER_DIRECTION=20`
- `DEFAULT_SEEDS_PER_STEP=20`

### Task granularity

- 1 parameter = 1 array task
- the array task count must cover all rows in `parameter_targets.tsv`

If the target parameter count exceeds the submitted array task count, the script stops with an explicit error instead of silently skipping parameters.

## 10. `estimate_live_effective_pms.R`

### Purpose

This script estimates the effective live-cell-view `p_ms` level from an existing run or seed, including:

- live-weighted effective `p_ms`
- harvest-only `p_ms`
- cohort-stratified `2N` and `4N` results

It reads existing fit results and visualization outputs. It does not refit the model.

### Recommended calls

#### Directly from a seed directory

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
  --seed_dir=oxygen/results/fit_invivo_o2_supply_demand_pmiss05_20260402_160038/seed7
```

#### From a run directory with an explicit seed

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
  --run_dir=oxygen/results/fit_invivo_o2_supply_demand_pmiss05_20260402_160038 \
  --seed=7
```

#### From a run directory with automatic seed selection

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/estimate_live_effective_pms.R \
  --run_dir=oxygen/results/fit_invivo_o2_supply_demand_pmiss05_20260402_160038
```

### Automatic seed selection logic

If only `--run_dir` is given, the script:

- reads `run_dir/extra_results/seed_summary.tsv`
- prefers:
  - `recommend_rank_burden_ploidy_boundary_first`
- falls back to other rank columns if needed
- falls back to objective ordering if no rank columns are available

### Optional arguments

- `--out_dir=...`
  - Defaults to `seed_dir/viz/live_effective_pms/`
- `--seed_id=seed7`
  - Alternative to `--seed=7`

### Main outputs

Default output directory:

- `seed_dir/viz/live_effective_pms/`

Files written there:

- `live_effective_pms_context.tsv`
- `live_effective_pms_sample_day.tsv`
- `live_effective_pms_overall.tsv`
- `live_effective_pms_harvest_only.tsv`
- `live_effective_pms_cohort_all_days.tsv`
- `live_effective_pms_cohort_harvest_only.tsv`

### Output interpretation

- `live_weighted_effective_p_ms`
  - instantaneously weighted by the live-cell composition
- `live_weighted_retained_p_ms_proxy`
  - a more conservative retained proxy that also incorporates `viability_after_ms`
- `harvest_only`
  - restricted to harvest timepoints
- `cohort_*`
  - stratified by `2N` and `4N`

## 11. `compare_sigma_burden_live_effective_pms.R`

### Purpose

This script compares seed-level fitted `p_misseg` values and live-cell effective `p_ms` values across multiple completed `sigma_burden` analyses.

It is intended for the situation where several completed run directories exist, for example:

- `oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_0.05`
- `oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_0.15`
- `oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_0.3`
- `oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_0.6`

For each `seed*` directory in each run, the script:

- reuses or runs `estimate_live_effective_pms.R`
- extracts:
  - fitted `p_misseg` parameter
  - overall live-cell effective `p_ms`
  - harvest-only live-cell effective `p_ms`
- writes a combined by-seed comparison table
- writes a per-`sigma_burden` summary table
- generates a violin-plus-boxplot figure comparing fitted `p_misseg` vs live-cell effective `p_ms`

### Recommended calls

#### Standard local or HPC comparison

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/compare_sigma_burden_live_effective_pms.R \
  --sigma_caps=0.05,0.15,0.3,0.6 \
  --run_dir_template=oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma} \
  --out_dir=oxygen/results/comp_live_effective_pms \
  --n_cores=8
```

#### Recompute all `estimate_live_effective_pms.R` outputs before summarizing

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/compare_sigma_burden_live_effective_pms.R \
  --sigma_caps=0.05,0.15,0.3,0.6 \
  --run_dir_template=oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma} \
  --out_dir=oxygen/results/comp_live_effective_pms \
  --n_cores=8 \
  --rerun_estimate=TRUE
```

#### Small smoke test

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/compare_sigma_burden_live_effective_pms.R \
  --sigma_caps=0.15 \
  --run_dir_template=oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma} \
  --out_dir=oxygen/results/comp_live_effective_pms_smoke \
  --n_cores=2 \
  --max_seeds=5
```

### Main arguments

- `--sigma_caps=...`
  - Comma-separated `sigma_burden` upper-bound values such as `0.05,0.15,0.3,0.6`.
- `--run_dir_template=...`
  - Run-directory template containing `{sigma}`.
  - Example:
    - `oxygen/results/fit_invivo_o2_supply_demand_MAP_pmiss_0.5_sigma_burden_{sigma}`
- `--out_dir=...`
  - Output directory for merged tables and comparison plots.
- `--n_cores=...`
  - Number of parallel workers.
  - On Unix-like systems this uses dynamic task scheduling via `mclapply(..., mc.preschedule=FALSE)`.
- `--rerun_estimate=TRUE|FALSE`
  - If `FALSE`, the script reuses existing:
    - `seed_dir/viz/live_effective_pms/`
  - If `TRUE`, it reruns `estimate_live_effective_pms.R` for every seed.
- `--max_seeds=...`
  - Optional cap for quick smoke tests.
- `--estimate_script=...`
  - Optional override path for `estimate_live_effective_pms.R`.

### Input requirements

Each `seed*` directory must either already contain:

- `viz/live_effective_pms/live_effective_pms_overall.tsv`
- `viz/live_effective_pms/live_effective_pms_harvest_only.tsv`

or, if `--rerun_estimate=TRUE` is used, it must contain the inputs required by `estimate_live_effective_pms.R`, namely:

- `best_params.tsv`
- `fit_config.rds`
- `viz/ploidy_timecourse.tsv`
- `viz/predict_burden_vs_o2.tsv`
- `viz/functional_curve_ploidy.tsv`

### Main outputs

The output directory contains:

- `sigma_burden_live_effective_pms_by_seed.tsv`
  - one row per seed
- `sigma_burden_live_effective_pms_summary.tsv`
  - one row per `sigma_burden` group
- `sigma_burden_live_effective_pms_task_manifest.tsv`
  - explicit `(sigma_burden, seed)` task list
- `sigma_burden_live_effective_pms_task_results.tsv`
  - per-task completion table, useful for debugging failures
- `sigma_burden_p_misseg_vs_live_cell_violin.pdf`
  - violin-plus-boxplot comparison of fitted `p_misseg` vs live-cell effective `p_ms`

### Output interpretation

- `p_misseg_parameter`
  - fitted model parameter from `best_params.tsv`
  - this is the `p_misseg` amplitude parameter, not directly the state-specific `p_ms`
- `live_cell_p_misseg`
  - overall live-cell effective `p_ms`
  - extracted from:
    - `live_weighted_effective_p_ms_mean`
    in `live_effective_pms_overall.tsv`
- `harvest_live_cell_p_misseg`
  - harvest-only live-cell effective `p_ms`
  - extracted from:
    - `live_weighted_effective_p_ms_mean`
    in `live_effective_pms_harvest_only.tsv`

### Failure behavior

If one or more seed tasks fail:

- the script writes `sigma_burden_live_effective_pms_task_results.tsv`
- then stops with an error message pointing to the first failed task

## 12. Internal Scripts

### `model_O2G_supply_demand_MAP.R`

- loads model logic and R/C++ interfaces
- usually sourced by the fitter or visualization scripts
- not intended to be used as a standalone entrypoint

### `model_O2G_supply_demand_MAP.cpp`

- implements the main C++ simulation and objective calculations
- includes burden, ploidy, and raw `-2logL` internals
- not run directly; called through the Rcpp interface from upper-layer scripts

### `o2g_supply_demand_map_common_semantics.R`

- stores shared semantics and configuration logic used by fit and simulation code
- internal shared layer

### `o2g_supply_demand_map_shared.R`

- provides shared utility functions such as:
  - CLI parsing
  - path resolution
  - type conversion
  - common helper functions
- sourced by nearly all entrypoint scripts in this directory

## Common Command Templates

### Lightweight one-seed smoke fit

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invivo \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --run_prefix=smoke_$(date +%Y%m%d_%H%M%S) \
  --seeds_csv=1 \
  --itermax=1 \
  --NP=12 \
  --n_cores=1 \
  --auto_viz=FALSE
```

### Re-run `extra_results` on a completed run

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R \
  --run_dir=oxygen/results/your_run
```

### Re-run visualization for one completed seed

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/vis/viz_invivo_model_O2G_supply_demand_MAP_results.R \
  --fit_dir=oxygen/results/your_run/seed1
```

### Local profile-likelihood smoke test

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/optimizer/profile_likelihood_O2G_supply_demand_MAP.R \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --baseline_seed_dir=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2 \
  --profile_bounds_table=oxygen/results/fit_invivo_o2_supply_demand_eq21_20260331_011709/seed2/parameter_table_input.csv \
  --output_root=oxygen/results/profile_smoke_local \
  --param_name=lam_min \
  --max_steps_per_direction=1 \
  --seeds_per_step=1 \
  --n_cores=1
```

### Collect profile-likelihood results

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/collect_profile_likelihood_results.R \
  --output_root=oxygen/results/profile_smoke_local
```

## Common Notes and Caveats

### Default parameter tables vs baseline bounds tables

Some profile and boundary-calibration workflows intentionally use:

- a specific baseline seed's `parameter_table_input.csv`

instead of the current repository-wide default tables:

- `oxygen/data/O2G_supply_demand/parameter_table_O2.csv`
- `oxygen/data/O2G_supply_demand/parameter_table_O2G.csv`

This is done to keep the workflow self-consistent with the selected baseline. These two tables should not be mixed casually.

### `objective` is not pure likelihood

In the current fitter:

- `objective`
  - is not pure likelihood
  - it is the current balanced objective plus optional prior penalty
- `objective_data`
  - is also not a raw likelihood sum
  - it is still a balanced or averaged data loss
- the newer diagnostics:
  - `objective_burden_neg2loglik_raw`
  - `objective_ploidy_neg2loglik_raw`
  are the separate raw `-2logL` reporting fields

### Error bars in all-seeds profile plots are not confidence intervals

For all-seeds profile plots:

- the error bars represent the seed sample range
- they are not confidence intervals

For raw `-2logL` profile CI interpretation:

- CI comes from threshold crossing logic
- not from the plotted seed-range error bars themselves

## Maintenance Suggestion

If new scripts are added later, the following four items should also be added to this README:

1. purpose
2. recommended invocation
3. key arguments
4. main outputs
