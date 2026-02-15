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
- `--pretreat_only` (default `TRUE`): use only `Day <= Day of 1st treatment` for burden fitting
- `--ploidy_at_harvest` (default `TRUE`): align terminal ploidy at harvest day
- `--max_scenarios` (default `Inf`): cap number of scenarios

### Optimization

- `--itermax` (default `40`)
- `--NP` (default `80`)
- `--seed` (default `1`)
- `--n_starts` (fallback `optim` only, default `20`)
- `--optim_maxit` (fallback `optim` only, default `max(200, itermax*50)`)

### Objective weights

- Global: `--w_burden` (default `1`), `--w_ploidy` (default `1`)
- Two-stage: `--stage1_w_burden`, `--stage1_w_ploidy`, `--stage2_w_burden`, `--stage2_w_ploidy`

Note: `objective` in `fit_summary.tsv` is computed with global `w_burden/w_ploidy`; stage-specific objective values are reported as `stage1_*` and `stage2_*`.

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
- `burden_fit.tsv`: observed vs predicted burden trajectories
- `terminal_ploidy_fit.tsv`: observed vs predicted terminal ploidy distribution
- `deoptim_result.rds`: raw optimizer outputs
- `fit_config.rds`: run configuration used for this fit

## 8. Recommended starting command

Use pretreatment burden data while keeping all dose groups:

```bash
Rscript /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/fit_invivo_ploidy_buffer.R \
  --dose_zero_only=FALSE \
  --pretreat_only=TRUE \
  --itermax=120 \
  --NP=140
```
