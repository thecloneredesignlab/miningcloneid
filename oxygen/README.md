# O2 Supply-Demand MAP Fitting Workflow

## R Package Requirements

Install the R packages used by the O2 supply-demand fitting, visualization,
reporting, analysis, and tests:

```r
install.packages(c(
  "Matrix", "Rcpp", "DEoptim", "yaml",
  "dplyr", "tidyr", "ggplot2", "readxl", "readr", "data.table",
  "patchwork", "cowplot", "ggh4x", "ggalluvial", "ggrepel", "gridExtra",
  "rmarkdown", "magick", "base64enc",
  "testthat"
))
```

This directory contains the in vivo, in vitro, and joint O2 supply-demand MAP
fitting workflow. The main implementation lives under
`code/O2_supply_demand_MAP/`.

The current workflow supports:

- separate in vivo fitting against tumor-burden and terminal-ploidy data;
- separate in vitro fitting against passage growth, ploidy, karyotype, and flow
  data;
- joint in vivo/in vitro fitting;
- joint soft coupling for selected parameters that should be similar, but not
  forced to be identical, between in vivo and in vitro contexts.

## Main Paths

- Unified optimizer entry point:
  `code/O2_supply_demand_MAP/optimizer/fit_model_O2_supply_demand_MAP.R`
- Local shell runner:
  `code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh`
- Joint-only shell runner:
  `code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh`
- Unified HPC submitter:
  `code/O2_supply_demand_MAP/hpc/submit_o2_fit.sh`
- Default O2 config:
  `config/O2_supply_demand.yaml`
- Warm-start joint submission:
  use `config/O2_supply_demand.yaml` with the joint warm-start seed-dir and
  soft-coupling start-table options described below.
- Labelled soft-coupling start tables:
  `data/O2_supply_demand/joint_soft_coupling_parameters_table__<seed_label>.csv`

## Requirements

The workflow is R-based and uses the packages listed above plus the in-house
helper code in this repo.

Run commands from the repository root unless an absolute path is used:

```bash
cd /Users/4482173/Documents/GitHub/soft_coupling
```

## Separate In Vivo Fitting

The in vivo fit is selected with `--fit_invivo`.

For normal use, run through YAML-driven runner mode. A one-seed local run is
just runner mode with `--seeds_csv=1`:

```bash
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh \
  --fit_invivo \
  --mode=run \
  --config=oxygen/config/O2_supply_demand.yaml \
  --out_root=oxygen/results \
  --run_prefix=fit_invivo_example \
  --append_run_prefix_timestamp=FALSE \
  --seeds_csv=1 \
  --data_dir=data/InVivoData_Gemcitabine
```

The in vivo backend loads:

- the in vivo parameter table from the config or default O2 table;
- tumor burden workbook `dt_Gem_VT_20260209_v5.xlsx`;
- terminal ploidy data resolved from `data_dir`;
- config flags such as `O2_growth`, `fit_treatment`, `fit_tau_O2`,
  `truncate_at_treatment`, `ploidy_at_harvest`, and `paired_only`.

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
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh \
  --fit_invivo \
  --mode=run \
  --config=oxygen/config/O2_supply_demand.yaml \
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
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh \
  --fit_invitro \
  --seed=1 \
  --out_dir=oxygen/results/fit_invitro_example/seed1 \
  --parameter_table=oxygen/data/O2_supply_demand/parameter_table_invitro_buffering.csv \
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

The optimizer also runs DEoptim followed by optional `L-BFGS-B` local
refinement. Main outputs include:

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
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh \
  --mode=run \
  --config=oxygen/config/O2_supply_demand.yaml
```

The equivalent unified-runner call is:

```bash
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_model_O2_supply_demand_MAP.sh \
  --fit_joint \
  --mode=run \
  --config=oxygen/config/O2_supply_demand.yaml
```

For one joint seed:

```bash
bash oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh \
  --mode=fit_seed \
  --seed=1 \
  --out_dir=oxygen/results/fit_joint_example/seed1 \
  --data_dir=data/InVivoData_Gemcitabine \
  --config=oxygen/config/O2_supply_demand.yaml
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
bash oxygen/code/O2_supply_demand_MAP/hpc/submit_o2_fit.sh \
  --fitting_mode=invivo \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --invivo_run_prefix=fit_invivo_O2_buffering_500seed \
  --invivo_total_seeds=500 \
  --invivo_array_tasks=500 \
  --invivo_seeds_per_task=1
```

```bash
bash oxygen/code/O2_supply_demand_MAP/hpc/submit_o2_fit.sh \
  --fitting_mode=invitro \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --invitro_run_prefix=fit_invitro_O2_buffering_500seed \
  --invitro_total_seeds=500 \
  --invitro_array_tasks=500 \
  --invitro_seeds_per_task=1
```

For joint fitting:

```bash
bash oxygen/code/O2_supply_demand_MAP/hpc/submit_o2_fit.sh \
  --fitting_mode=joint \
  --joint_fitting_mode=DIRECT \
  --config_path=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/config/O2_supply_demand.yaml \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --joint_run_prefix=fit_joint_O2_buffering_500seed \
  --joint_job_name=o2_joint_B \
  --invivo_best_seed_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invivo_O2_buffering_500seed/seed50 \
  --invitro_best_seed_dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results/fit_invitro_O2_buffering_500seed/seed350 \
  --joint_warmup_seed_label=invivo_seed50__invitro_seed350 \
  --joint_soft_coupling_sigma_default=0.35 \
  --joint_soft_coupling_huber_k=1.5 \
  --joint_total_seeds=500 \
  --joint_array_tasks=500 \
  --joint_seeds_per_task=1
```

When both best-seed directories are provided, the submitter first runs
`make_joint_soft_coupling_parameters_table.R`, writes a labelled start table
under `data/O2_supply_demand/`, passes that exact CSV path to the joint array,
and appends the same label to `joint_run_prefix` unless it is already present.
If `--joint_warmup_seed_label` is omitted, the label is inferred from the two
seed directory basenames, for example `invivo_seed50__invitro_seed350`.
Use `--joint_soft_coupling_sigma_default=<value>` and
`--joint_soft_coupling_huber_k=<value>` to override the config's soft-coupling
penalty settings for that submission without creating a separate YAML file.
Include those values in `--joint_warmup_seed_label` when submitting multiple
penalty settings in parallel.

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
joint_soft_coupling_params: O2_crit,mu_hp,p_misseg,k_o_mis,buffer_smax,buffer_beta,buffer_n_exp,n_O,alpha_o2,gamma_growth,lam_max,p_mis_base,p_wgd,gamma_mu
joint_soft_coupling_sigma_default: 0.35
joint_soft_coupling_huber_k: 1.5
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
| `alpha_o2` | `log10_alpha_o2` |
| `gamma_growth` | `gamma_growth` |
| `lam_max` | `log10_lam_max` |
| `p_mis_base` | `log10_p_mis_base` |
| `p_wgd` | `log10_p_wgd` |
| `gamma_mu` | `gamma_mu` |

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

The soft-coupling penalty is Huber-style on the standardized transformed-scale
split:

```text
z = delta / sigma_delta

penalty = 0.5 * z^2                    if |z| <= huber_k
penalty = huber_k * (|z| - 0.5*huber_k) if |z| > huber_k
```

where `sigma_delta` comes from `joint_soft_coupling_sigma_<parameter>` when that
parameter-specific key exists, otherwise from
`joint_soft_coupling_sigma_default`. `joint_soft_coupling_huber_k` controls how
many sigma units are treated quadratically; beyond `huber_k * sigma_delta`,
additional separation is penalized linearly.

## How Joint Bounds Are Determined

Before soft coupling is applied, the joint backend computes merged natural
bounds for all shared parameters. For hard-shared parameters, it keeps the
existing merged behavior:

```text
joint_lower = min(invivo_lower, invitro_lower)
joint_upper = max(invivo_upper, invitro_upper)
```

For soft-coupled parameters, each backend's natural lower and upper bounds are
transformed onto the optimizer scale and then combined into a single joint union
bound:

```text
joint_union_lower_t = min(vivo_bound_lower_t, vitro_bound_lower_t)
joint_union_upper_t = max(vivo_bound_upper_t, vitro_bound_upper_t)
```

The per-backend transformed bounds are not stored in soft-coupling metadata or
reported in soft-coupling summaries. The joint union bound is the only active
admissibility rule during joint fitting. The center optimizer bound is the joint
union bound. The delta optimizer bound is the full transformed union span in
either direction:

```text
delta_abs = joint_union_upper_t - joint_union_lower_t
delta_lower_t = -delta_abs
delta_upper_t =  delta_abs
```

During evaluation, context-specific values are reconstructed directly:

```text
vivo_t  = center + delta / 2
vitro_t = center - delta / 2
```

No runtime clipping is applied. If either reconstructed value lies outside the
joint union transformed bound, the optimizer point is infeasible and receives a
large penalty before either likelihood component is evaluated.

## Warm-Start and Start-Table Handling

To generate a labelled joint soft-coupling start table directly from separate
best-seed directories:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/make_joint_soft_coupling_parameters_table.R \
  --invivo-seed-dir oxygen/results/fit_invivo_O2_buffering_500seed/seed50 \
  --invitro-seed-dir oxygen/results/fit_invitro_O2_buffering_500seed/seed350 \
  --seed-label invivo_seed50__invitro_seed350
```

The default output is:

```text
oxygen/data/O2_supply_demand/joint_soft_coupling_parameters_table__invivo_seed50__invitro_seed350.csv
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

- warm-start values must already fall inside the optimizer bounds;
- soft-coupled warm starts must reconstruct in vivo and in vitro transformed
  values inside the joint union bound;
- invalid warm starts stop with a clear error instead of clipping or expanding
  bounds;
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

- start-table values are init overrides only;
- start-table values must already fall inside the optimizer bounds;
- soft-coupled start-table values must reconstruct in vivo and in vitro
  transformed values inside the joint union bound;
- invalid start-table rows stop with a clear error instead of clipping or
  expanding bounds;
- every applied row records the unchanged before/after bounds and `bound_action`
  in `joint_soft_coupling_initial_values.tsv`.

The final start table currently includes explicit center and delta values for
the soft-coupled parameters listed above, plus any remaining non-soft-coupled
shared starts.

## Feasibility During Objective Evaluation

For every objective evaluation, the joint backend reconstructs and reports:

```text
vivo_transformed
vitro_transformed
joint_union_lower_transformed
joint_union_upper_transformed
feasible_at_point
```

These fields are written to `joint_soft_coupling.tsv` for accepted solutions.
If a trial point leaves the joint union bound, the point is rejected with a
large objective penalty before either backend likelihood is evaluated.

## Joint Soft-Coupling Outputs

A joint seed directory can contain:

- `joint_soft_coupling.tsv`: center/delta values, context-specific transformed
  and natural values, ratios/fold changes, penalty paid, joint union bounds, and
  feasibility status;
- `joint_components.tsv`: objective components including
  `objective_soft_coupling`;
- `fit_summary.tsv`: summary rows such as `joint_soft_coupling_enabled`,
  `joint_soft_coupling_params`, `joint_soft_coupling_n_params`,
  `joint_soft_coupling_sigma_default`, and `joint_soft_coupling_huber_k`;
- `joint_soft_coupling_initial_values.tsv`: start-table init overrides;
- `joint_warmup_initial_values.tsv`: warm-start sources, values, and bound
  actions;
- `fit_config.rds` and `fit_result.rds`.

The report renderer and extra-results workflow read these files to display
soft-coupled parameter diagnostics separately from the old hard-shared
parameter table.

## Postprocessing and Reports

Extra-results aggregation:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/extra_results.R \
  --run_dir=oxygen/results/fit_joint_O2_buffering_500seed
```

Extra-results HTML report:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/analysis/extra_results_report.R \
  --extra_results_dir=oxygen/results/fit_joint_O2_buffering_500seed/extra_results
```

Per-seed fit report:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R \
  --fit_dir=oxygen/results/fit_joint_O2_buffering_500seed/seed1
```

The extra-results workflow also writes combined soft-coupling tables and plots
when `joint_soft_coupling.tsv` is present in seed directories.
