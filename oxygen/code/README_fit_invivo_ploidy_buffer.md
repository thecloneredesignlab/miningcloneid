# fit_invivo_ploidy_buffer User Guide

Script path: `/Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R`

## 1. What this script does

This script estimates parameters of the `ploidy_buffer` ODE model using:

- in vivo tumor volume time series
- terminal ploidy distributions

By default, it runs in **two-stage fitting mode**:

- Stage 1 (growth): fit `log10_R`, `beta`, `log10_eta`
- Stage 2 (ploidy): fix Stage 1 values, then fit remaining parameters (for example `log10_pwgd`, `mr0`, `mr1`, `log10_pmis`)

## 2. Input data

Default input directory (override with `--data_dir`):

`/Users/4482173/Documents/GitHub/miningcloneid/data/InVivoData_Gemcitabine`

Expected files:

- `dt_Gem_VT_20260209_v5.xlsx`
- `all_ploidy.tsv`

## 3. Dependencies

Required R packages:

- `Matrix`
- `dplyr`
- `tidyr`
- `readxl`

Optional:

- `DEoptim` (if unavailable, the script falls back to multi-start `optim(L-BFGS-B)`)

## 4. Basic usage

Run from repo root:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R
```

Set output directory:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --out_dir=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/my_run
```

Append timestamp to a custom output directory:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --out_dir=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/my_run \
  --append_timestamp_out_dir=TRUE
```

## 5. Two-stage vs single-stage

Default: `--two_stage=TRUE`

- Stage 1 defaults: `stage1_w_burden=1`, `stage1_w_ploidy=0`
- Stage 2 defaults: `stage2_w_burden=0`, `stage2_w_ploidy=1`

Switch back to single-stage joint fitting:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --two_stage=FALSE
```

## 6. Key arguments (current script defaults)

### Data filtering

- `--dose_zero_only` (default `TRUE`): use only `Dose=0` rows
- `--truncate_at_treatment` (default `FALSE`): if `TRUE`, use only `Day <= Day of 1st treatment` for burden fitting
- `--ploidy_at_harvest` (default `TRUE`): align terminal ploidy at harvest day
- `--max_scenarios` (default `Inf`): cap number of scenarios

Backward-compatible alias:

- `--pretreat_only` (legacy alias for `--truncate_at_treatment`)

### Optimization

- `--itermax` (default `40`)
- `--NP` (default `80`)
- `--n_cores` (default `auto`): number of CPU cores for parallel optimization.  
  Default auto strategy is `max(1, detectCores() - 1)` (leave one core free).
- `--seed` (default `1`)
- `--n_starts` (fallback `optim` only, default `20`)
- `--optim_maxit` (fallback `optim` only, default `max(200, itermax*50)`)
- `--optim_trace` (default `TRUE`): print clean multi-start progress for optim backend
- `--optim_trace_every` (default `1`): print every N starts (for example, `5` prints 5/10/15/...)
- `--use_deoptim` (default `TRUE`): enable DEoptim backend
- `--deoptim_parallel` (default `FALSE`): allow DEoptim parallel workers when `n_cores > 1`
- `--init_params_tsv` (optional): warm-start file for optimizer initialization.  
  Supported formats:
  - `best_params.tsv` (`parameter`, `value`)
  - `fit_parameter_stages.tsv` (`transformed_parameter`, `transformed_value`)
  - one-row transformed parameter table
- `--append_timestamp_out_dir` (default `FALSE`): append timestamp suffix to `--out_dir`
- `--timestamp_format` (default `%Y%m%d_%H%M%S`): timestamp format used in output directory naming

Parallel behavior:

- Default behavior with `n_cores > 1`: use parallel multi-start `optim` backend.
- DEoptim is used by default for `n_cores = 1`.
- If you explicitly set `--deoptim_parallel=TRUE`, the script will try DEoptim parallel mode.
- If DEoptim is unavailable or fails, the script falls back to multi-start `optim`.

Progress behavior for `optim` backend:

- In serial mode, logs are printed as starts finish, e.g. `start 3/20 finished: val=..., best=...`
- In parallel mode, starts are evaluated concurrently, so logs are reported after collection.

### Objective weights

- Global (single-stage): `--w_burden`, `--w_ploidy`
- Two-stage: `--stage1_w_burden`, `--stage1_w_ploidy`, `--stage2_w_burden`, `--stage2_w_ploidy`

For single-stage mode, `--w_burden` and `--w_ploidy` can be scalars or comma-separated arrays:

- Scalar example: `--w_burden=1 --w_ploidy=1` (one optimization pass)
- Array example: `--w_burden=0.1,1 --w_ploidy=5,1` (two passes in one run)

Array behavior in single-stage mode:

- pass 1 uses first values
- pass 2 warm-starts from pass 1 and uses second values
- and so on for longer arrays

Rules:

- arrays must have same length, or one side can be length 1 (auto-recycled)
- weight arrays are not allowed when `--two_stage=TRUE`

Note: `objective` in `fit_summary.tsv` is reported using the final pass weights.  
Weight schedule details are also recorded in `weight_passes`, `w_burden_schedule`, and `w_ploidy_schedule`.

### Other model/numerical settings

- `--dt` (default `0.5`)
- `--O2` (default `1.0`)
- `--K` (default `1e12`)
- `--crowding` (`logistic` or `gompertz`, default `logistic`)
- `--init_total_size` (default `1e6`)
- `--N_UNIT` (default `22`)
- `--N_MIN` (default `22`)
- `--N_MAX` (default `154`)
- `--dose_ref` (default `30`)
- `--tx_mult_min` (default `0.05`)
- `--huber_k` (default `0.1`)
- `--eps_prob` (default `1e-12`)
- `--trace_obj` (default `FALSE`)
- `--fit_full_pmis` (default `FALSE`)
- `--fit_treatment` (default `FALSE`)

## 7. Output files

Default output location:

`/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_invivo_ploidy_buffer_YYYYmmdd_HHMMSS`

Main files:

- `best_params.tsv`: final fitted parameters (real scale)
- `stage1_best_params.tsv`: parameters at end of Stage 1 (two-stage only)
- `fit_parameter_stages.tsv`: which transformed parameter was optimized in which stage
- `fit_summary.tsv`: fit summary and objective decomposition
- `single_stage_pass_summary.tsv`: per-pass objective trace for single-stage weight schedule mode
- `burden_fit.tsv`: observed vs predicted burden trajectories
- `terminal_ploidy_fit.tsv`: observed vs predicted terminal ploidy distribution
- `deoptim_result.rds`: raw optimizer outputs
- `fit_config.rds`: run configuration used for this fit

## 8. Recommended starting command

Use pretreatment burden data while keeping all dose groups:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --dose_zero_only=FALSE \
  --truncate_at_treatment=TRUE \
  --n_cores=4 \
  --itermax=120 \
  --NP=140
```

## 9. Warm-start workflow (single-stage -> single-stage)

Run 1: ploidy-prioritized single-stage fit.

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --two_stage=FALSE \
  --w_burden=0.1 --w_ploidy=5 \
  --out_dir=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_joint_ploidy_bias
```

Run 2: equal-weight single-stage fit, warm-start from run 1.

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --two_stage=FALSE \
  --w_burden=1 --w_ploidy=1 \
  --init_params_tsv=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_joint_ploidy_bias/best_params.tsv \
  --out_dir=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_joint_equal_weight_warm
```

## 10. One-run two-pass single-stage optimization

This runs both passes in a single execution:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --two_stage=FALSE \
  --w_burden=0.1,1 \
  --w_ploidy=5,1 \
  --n_cores=4 \
  --out_dir=/Users/4482173/Documents/GitHub/miningcloneid/oxygen/results/fit_joint_schedule
```

In this example:

- pass 1 is ploidy-prioritized (`0.1, 5`)
- pass 2 warm-starts from pass 1 and uses equal weights (`1, 1`)
