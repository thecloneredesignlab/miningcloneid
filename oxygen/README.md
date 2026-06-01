# O2G Supply-Demand MAP Fitting Workflow

This directory contains the in vivo, in vitro, and joint O2G supply-demand MAP
fitting workflow. The main implementation lives under
`code/O2G_supply_demand_MAP/`.

The current workflow supports:

- separate in vivo fitting against tumor-burden and terminal-ploidy data;
- separate in vitro fitting against passage growth, ploidy, karyotype, and flow
  data;
- joint in vivo/in vitro fitting;
- joint soft coupling for selected parameters that should be similar, but not
  forced to be identical, between in vivo and in vitro contexts.

## Main Paths

- Unified optimizer entry point:
  `code/O2G_supply_demand_MAP/optimizer/fit_model_O2G_supply_demand_MAP.R`
- Local shell runner:
  `code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh`
- Joint-only shell runner:
  `code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh`
- Unified HPC submitter:
  `code/O2G_supply_demand_MAP/hpc/submit_o2g_fit.sh`
- Default O2G config:
  `config/O2G_supply_demand.yaml`
- Warm-start joint config:
  `config/O2G_supply_demand_warmup_seed50_seed350.yaml`
- Labelled soft-coupling start tables:
  `data/O2G_supply_demand/joint_soft_coupling_parameters_table__<seed_label>.csv`

## Requirements

The workflow is R-based and uses packages such as `DEoptim`, `Matrix`, `dplyr`,
`ggplot2`, `data.table`, `yaml`, and the in-house helper code in this repo.

Run commands from the repository root unless an absolute path is used:

```bash
cd /Users/4482173/Documents/GitHub/soft_coupling
```

## Separate In Vivo Fitting

The in vivo fit is selected with `--fit_invivo`.

For normal use, run through YAML-driven runner mode. A one-seed local run is
just runner mode with `--seeds_csv=1`:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invivo \
  --mode=run \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --out_root=oxygen/results \
  --run_prefix=fit_invivo_example \
  --append_run_prefix_timestamp=FALSE \
  --seeds_csv=1 \
  --data_dir=data/InVivoData_Gemcitabine
```

The in vivo backend loads:

- the in vivo parameter table from the config or default O2G table;
- tumor burden workbook `dt_Gem_VT_20260209_v5.xlsx`;
- terminal ploidy data resolved from `data_dir`;
- config flags such as `glucose`, `O2_growth`, `fit_treatment`,
  `fit_tau_O2`, `truncate_at_treatment`, `ploidy_at_harvest`, and
  `paired_only`.

The optimizer runs DEoptim and then attempts a serial `L-BFGS-B` refinement from
the DEoptim best point. Outputs are written under the seed directory, including:

- `fit_summary.tsv`
- `best_params.tsv`
- `best_params_transformed.tsv` or stage/checkpoint parameter tables when
  available
- `fit_result.rds`
- optional visualization and HTML report outputs when auto-reporting is enabled

For a config-driven multi-seed local run, use the same mode with multiple seeds:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invivo \
  --mode=run \
  --config=oxygen/config/O2G_supply_demand.yaml \
  --seeds_csv=1,2,3
```

In run mode, the config or CLI overrides must provide `run_prefix`, `out_root`,
`data_dir`, and seed settings such as `seeds_csv` or `seeds_file`. The runner
creates:

```text
<out_root>/<run_prefix>/seed<seed>/
```

for each seed, snapshots the resolved config and parameter table, runs fitting,
then optionally runs visualization and report rendering.

## Separate In Vitro Fitting

The in vitro fit is selected with `--fit_invitro`.

For a single seed:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_invitro \
  --seed=1 \
  --out_dir=oxygen/results/fit_invitro_example/seed1 \
  --parameter_table=oxygen/data/O2G_supply_demand/parameter_table_invitro_buffering.csv \
  --fit_objects_dir=oxygen/ploidyOxygen/data/fit_objects \
  --flow_density_path=oxygen/data/g0g1_ploidy_density_grid.csv \
  --itermax=500 \
  --NP=80 \
  --n_cores=1
```

The in vitro backend validates and loads:

- the in vitro parameter table from `--parameter_table` or the default in vitro
  table. The HPC submitter accepts aliases such as `--invitro_parameter_table`
  and forwards them as `--parameter_table`;
- fit objects from `fit_objects_dir`;
- optional flow-density grid from `flow_density_path`;
- fixed-oxygen settings and in vitro numerical settings such as `dt`,
  `init_total_size`, and `o2_upper_bound`.

In vitro fitting keeps `glucose=FALSE`; if a glucose flag is passed, it is
ignored for this mode. The optimizer also runs DEoptim followed by optional
`L-BFGS-B` local refinement. Main outputs include:

- `fit_summary.tsv`
- `best_params.tsv`
- `best_params_transformed.tsv`
- `fit_result.rds`
- optional in vitro visualization/report files

For many in vitro seeds on HPC, use the unified submitter shown below instead of
manually launching each seed.

## Joint Fitting

Joint fitting is selected with `--fit_joint`. The joint-only shell wrapper adds
`--fit_joint` automatically:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh \
  --mode=run \
  --config=oxygen/config/O2G_supply_demand.yaml
```

The equivalent unified-runner call is:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_model_O2G_supply_demand_MAP.sh \
  --fit_joint \
  --mode=run \
  --config=oxygen/config/O2G_supply_demand.yaml
```

For one joint seed:

```bash
bash oxygen/code/O2G_supply_demand_MAP/runner/run_fit_joint_model_O2G_supply_demand_MAP.sh \
  --mode=fit_seed \
  --seed=1 \
  --out_dir=oxygen/results/fit_joint_example/seed1 \
  --data_dir=data/InVivoData_Gemcitabine \
  --config=oxygen/config/O2G_supply_demand.yaml
```

In `--mode=run`, the joint runner reads the config, resolves the seed plan, then
calls the optimizer once per seed. For each seed it runs:

1. joint fitting;
2. in vivo visualization;
3. in vitro visualization;
4. joint HTML report rendering.

The joint objective is:

```text
objective_total =
  joint_weight_invivo  * invivo_objective
  + joint_weight_invitro * invitro_objective
  + objective_soft_coupling
  + objective_constraints
```

The in vivo part is evaluated by the in vivo backend. The in vitro part is
evaluated by constructing the in vitro transformed parameter vector from the
joint state, then applying in vitro-only optimizer parameters and any
soft-coupled in vitro-specific values.

## HPC Submission

Use the unified HPC submitter for production multi-seed runs:

```bash
bash oxygen/code/O2G_supply_demand_MAP/hpc/submit_o2g_fit.sh \
  --fitting_mode=invivo \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2G_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --invivo_run_prefix=fit_invivo_O2G_buffering_500seed \
  --invivo_total_seeds=500 \
  --invivo_array_tasks=500 \
  --invivo_seeds_per_task=1
```

```bash
bash oxygen/code/O2G_supply_demand_MAP/hpc/submit_o2g_fit.sh \
  --fitting_mode=invitro \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2G_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --invitro_run_prefix=fit_invitro_O2G_buffering_500seed \
  --invitro_total_seeds=500 \
  --invitro_array_tasks=500 \
  --invitro_seeds_per_task=1
```

For joint fitting:

```bash
bash oxygen/code/O2G_supply_demand_MAP/hpc/submit_o2g_fit.sh \
  --fitting_mode=joint \
  --joint_fitting_mode=DIRECT \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2G_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --joint_run_prefix=fit_joint_O2G_buffering_500seed \
  --joint_job_name=o2g_joint_B \
  --invivo_best_seed_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invivo_O2G_buffering_500seed/seed50 \
  --invitro_best_seed_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invitro_O2G_buffering_500seed/seed350 \
  --joint_warmup_seed_label=invivo_seed50__invitro_seed350 \
  --joint_total_seeds=500 \
  --joint_array_tasks=500 \
  --joint_seeds_per_task=1
```

When both best-seed directories are provided, the submitter first runs
`make_joint_soft_coupling_parameters_table.R`, writes a labelled start table
under `data/O2G_supply_demand/`, passes that exact CSV path to the joint array,
and appends the same label to `joint_run_prefix` unless it is already present.
If `--joint_warmup_seed_label` is omitted, the label is inferred from the two
seed directory basenames, for example `invivo_seed50__invitro_seed350`.

`joint_fitting_mode` has these meanings:

- `OFF`: do not submit joint fitting; this is forced when `fitting_mode` is not
  `joint`.
- `DIRECT`: submit only the current joint fitter directly from the config.
- `JOINT`: submit in vivo and in vitro fits first, run extra-results
  postprocessing for each, then submit the current joint fitter.

After each fitting job finishes, a dependent postprocess job runs
`extra_results.R`. Existing extra-results outputs are skipped unless
`--force_extra_results=TRUE`.

Slurm stdout and stderr for jobs launched by the unified submitter are written
under `<out_root>/log` by default. Use `--log_root=/path/to/logs` to override
that location.

## Joint Soft-Coupled Parameters

Soft coupling is controlled by these config keys:

```yaml
joint_soft_coupling_enable: TRUE
joint_soft_coupling_params: O2_crit,mu_hp,p_misseg,k_o_mis,buffer_smax,buffer_beta,buffer_n_exp,n_O,lam_max,p_mis_base,p_wgd,gamma_mu
joint_soft_coupling_sigma_default: 0.35
joint_soft_coupling_delta_span_frac: 0.5
```

The currently soft-coupled parameters are:

| Natural parameter | Optimizer center name |
| --- | --- |
| `O2_crit` | `log10_O2_crit` |
| `mu_hp` | `log10_mu_hp` |
| `p_misseg` | `log10_p_misseg` |
| `k_o_mis` | `log10_k_o_mis` |
| `buffer_smax` | `buffer_smax` |
| `buffer_beta` | `log10_buffer_beta` |
| `buffer_n_exp` | `log10_buffer_n_exp` |
| `n_O` | `n_O` |
| `lam_max` | `log10_lam_max` |
| `p_mis_base` | `log10_p_mis_base` |
| `p_wgd` | `log10_p_wgd` |
| `gamma_mu` | `gamma_mu` |

`alpha_o2` and `gamma_growth` are still joint-shared parameters, but they are
not soft-coupled in the default config.

For each soft-coupled parameter, the optimizer uses two variables:

```text
center_name
delta__center_name
```

The context-specific transformed values are:

```text
in_vivo_transformed_raw  = center + delta / 2
in_vitro_transformed_raw = center - delta / 2
```

The soft-coupling penalty is:

```text
penalty = delta^2 / (2 * sigma_delta^2)
```

where `sigma_delta` comes from `joint_soft_coupling_sigma_<parameter>` when that
parameter-specific key exists, otherwise from
`joint_soft_coupling_sigma_default`.

## How Center Bounds Are Determined

Before soft coupling is applied, the joint backend computes merged natural
bounds for all shared parameters. For hard-shared parameters, it keeps the
existing merged behavior:

```text
joint_lower = min(invivo_lower, invitro_lower)
joint_upper = max(invivo_upper, invitro_upper)
```

For soft-coupled parameters, the center is not allowed to use the full union.
Instead, both backends' natural bounds are transformed onto the optimizer scale,
and the center uses the overlap:

```text
center_lower_t = max(invivo_lower_t, invitro_lower_t)
center_upper_t = min(invivo_upper_t, invitro_upper_t)
```

This guarantees that when `delta = 0`, both the in vivo and in vitro derived
values are equal to the center and are feasible for both backends. If no
transformed-scale overlap exists, joint context construction stops with an
error instead of silently using invalid bounds.

The initial center is the in vivo optimizer initial value clipped into the
center overlap:

```text
center_init_t = clip(invivo_init_t, center_lower_t, center_upper_t)
```

## How Delta Bounds Are Determined

Delta bounds are symmetric around zero and are initialized from the transformed
center span:

```text
center_span = center_upper_t - center_lower_t
delta_abs = joint_soft_coupling_delta_span_frac * center_span
delta_lower_t = -delta_abs
delta_upper_t =  delta_abs
```

With the default `joint_soft_coupling_delta_span_frac: 0.5`, the initial delta
range allows a diagnostic split up to half of the center range in either
direction.

During evaluation, the raw context values are clipped to each backend's own
transformed bounds:

```text
vivo_t  = clip(center + delta / 2, invivo_lower_t,  invivo_upper_t)
vitro_t = clip(center - delta / 2, invitro_lower_t, invitro_upper_t)
```

The clipped values, not the raw values, are decoded and passed to the in vivo
and in vitro objectives.

## Warm-Start and Start-Table Handling

To generate a labelled joint soft-coupling start table directly from separate
best-seed directories:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R \
  --invivo-seed-dir oxygen/results/fit_invivo_O2G_buffering_500seed/seed50 \
  --invitro-seed-dir oxygen/results/fit_invitro_O2G_buffering_500seed/seed350 \
  --seed-label invivo_seed50__invitro_seed350
```

The default output is:

```text
oxygen/data/O2G_supply_demand/joint_soft_coupling_parameters_table__invivo_seed50__invitro_seed350.csv
```

The CSV contains `param_name`, `value`, `scale`, `seed_label`,
`invivo_seed_label`, and `invitro_seed_label`. The joint backend reads the
`param_name`, `value`, and `scale` columns and ignores the label columns for
optimization; the labels are there to keep parallel warm-start submissions
traceable.

When `joint_warmup_enable: TRUE`, the joint backend reads selected in vivo and
in vitro best-seed parameter tables. For a soft-coupled parameter:

```text
center = (best_invivo_transformed + best_invitro_transformed) / 2
delta  = best_invivo_transformed - best_invitro_transformed
```

For hard-shared parameters, the warm-start value is the transformed-scale mean
of the in vivo and in vitro best values. In vivo-only and in vitro-only
parameters are initialized from their own best seed.

Warm-start bound behavior:

- center and other non-delta warm-start values are clipped into their current
  optimizer bounds;
- delta warm-start values are allowed to expand the delta lower/upper bound so
  the separate-fit difference can be represented;
- every applied warm-start value records its source and bound action in
  `joint_warmup_initial_values.tsv`.

After warm-start initialization, the optional labelled
`joint_soft_coupling_parameters_table__<seed_label>.csv` is applied as the final
override. It accepts `param_name`, `value`, and `scale` columns. Supported
scales are:

- `transformed`
- `log10`
- `identity`
- `natural`

For center parameters, `natural` values are transformed to the optimizer scale.
For delta parameters on log10 scales, a natural-scale difference is converted to
an optimizer-scale log difference using the current center value.

Start-table bound behavior:

- if the start-table value is below the current optimizer lower bound, the lower
  bound is expanded to that value;
- if it is above the current optimizer upper bound, the upper bound is expanded;
- expanded center bounds are copied back into the soft-coupling metadata;
- expanded delta bounds are copied back into the metadata delta bounds;
- every applied row records the before/after bounds and `bound_action` in
  `joint_soft_coupling_initial_values.tsv`.

The final start table currently includes explicit center and delta values for
the soft-coupled parameters listed above, plus non-soft-coupled shared starts
such as `log10_alpha_o2` and `gamma_growth`.

## Backend Bound Clipping During Objective Evaluation

Optimizer-bound expansion allows the optimizer to start from useful
separate-fit-derived values. It does not bypass backend validity constraints.

For every objective evaluation, the joint backend reports and applies:

```text
vivo_raw_transformed
vitro_raw_transformed
vivo_transformed
vitro_transformed
vivo_clipped
vitro_clipped
boundary_status_vivo
boundary_status_vitro
```

These fields are written to `joint_soft_coupling.tsv`. If a raw value exceeds a
backend-specific transformed range, the effective context value is clipped and
the boundary status is recorded as `clipped_lower` or `clipped_upper`. If it is
exactly on a boundary, the status is `at_lower` or `at_upper`; otherwise it is
`inside`.

## Joint Soft-Coupling Outputs

A joint seed directory can contain:

- `joint_soft_coupling.tsv`: center/delta values, context-specific transformed
  and natural values, ratios/fold changes, penalty paid, and clipping status;
- `joint_components.tsv`: objective components including
  `objective_soft_coupling`;
- `fit_summary.tsv`: summary rows such as `joint_soft_coupling_enabled`,
  `joint_soft_coupling_params`, `joint_soft_coupling_n_params`,
  `joint_soft_coupling_sigma_default`, and
  `joint_soft_coupling_delta_span_frac`;
- `joint_soft_coupling_initial_values.tsv`: start-table overrides and bound
  expansions;
- `joint_warmup_initial_values.tsv`: warm-start sources, values, and bound
  actions;
- `fit_config.rds` and `fit_result.rds`.

The report renderer and extra-results workflow read these files to display
soft-coupled parameter diagnostics separately from the old hard-shared
parameter table.

## Postprocessing and Reports

Extra-results aggregation:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/extra_results.R \
  --run_dir=oxygen/results/fit_joint_O2G_buffering_500seed
```

Extra-results HTML report:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/analysis/extra_results_report.R \
  --extra_results_dir=oxygen/results/fit_joint_O2G_buffering_500seed/extra_results
```

Per-seed fit report:

```bash
Rscript oxygen/code/O2G_supply_demand_MAP/report/render_fit_report.R \
  --fit_dir=oxygen/results/fit_joint_O2G_buffering_500seed/seed1
```

The extra-results workflow also writes combined soft-coupling tables and plots
when `joint_soft_coupling.tsv` is present in seed directories.
