# run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh

This script runs `fit_invivo_model_buffering_align_with_Richard.R` with an iterative stage-wise warm-start chain plus a final equal-weight callback pass.

Workflow:

1. Stage-wise chain passes (single-stage fitting, `two_stage=FALSE`)
2. Warm-start each pass from the previous pass output (`fit_parameter_stages.tsv`)
3. Final callback pass (default `(w_burden, w_ploidy) = (1,1)`), initialized from `step11` by default
4. Multi-seed execution
5. Optional auto-tuning of iteration / optimizer settings from input data complexity
6. Optional resume-from-pass and skip-existing support

Default oxygen baseline:

- If `--O2` is not provided, the fit script now uses `O2_fixed=5` (5%).

Burden observation model (current default in the fit script):

- burden data are fit in volume space (`mm^3`) via a ploidy-aware observation layer
- single-cell volume is modeled via estimated `rho_2N` (2N effective cell density; `cells/mm^3`) and `beta_size`
- `c_vol_2N_mm3` is retained only for legacy warm-start compatibility (`c_scale` -> `rho_2N` conversion)
- burden loss is `log-volume Huber` (not the older normalized cell-count shape loss)

Default chain in this runner:

- `w_burden_chain = 1,0.8,0.6,0.4,0.3,0.25,0.2,0.175,0.15,0.1,0.05,0`
- `w_ploidy_chain = 0,0.2,0.4,0.6,0.7,0.75,0.8,0.825,0.85,0.9,0.95,1`
- callback: `(1,1)`
- callback warm-start: `step11` by default (`callback_init_pass=11`)
- callback loss-scale mode: auto re-estimate by default (`callback_auto_rescale=TRUE`)

## Run

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh
```

## Pass Parameters Directly (`--key=value`)

You can pass parameters directly to the runner:

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh \
  --run_prefix=fit_stagewise_chain_new_deoptim_02212026_001 \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --seeds_csv=1,2,3 \
  --n_cores=24 \
  --fit_treatment=FALSE \
  --paired_only=TRUE \
  --auto_tune_iters=FALSE \
  --w_burden_chain=1,0.8,0.6,0.4,0.3,0.25,0.2,0.175,0.15,0.1,0.05,0 \
  --w_ploidy_chain=0,0.2,0.4,0.6,0.7,0.75,0.8,0.825,0.85,0.9,0.95,1 \
  --callback_init_pass=11 \
  --callback_auto_rescale=TRUE \
  --callback_w_burden=1 \
  --callback_w_ploidy=1 \
  --loss_rescale=TRUE \
  --use_deoptim=TRUE \
  --deoptim_parallel=TRUE \
  --pass_itermax=220 \
  --callback_itermax=480 \
  --np=260 \
  --pass_n_starts=80 \
  --callback_n_starts=140 \
  --pass_optim_maxit=15000 \
  --callback_optim_maxit=28000 \
  --rho_2N_min=3.2e4 \
  --rho_2N_max=5.6e4 \
  --burden_log_eps=1e-12 \
  --huber_k_burden_log=0.1
```

Note:

- For the current `log-volume Huber` burden loss, old manually tuned values such as `loss_scale_burden=0.003` / `loss_scale_ploidy=17.6477` may no longer be appropriate.
- Recommended first pass: keep `loss_rescale=TRUE` and omit manual `loss_scale_*` so the fit script can estimate scales for the new loss definition.

## Supported Keys

Supported `--key=value` options (same as the runner usage):

- `out_root`, `run_prefix`, `data_dir`, `seeds_csv`, `k`, `n_cores`, `max_scenarios`
- `pass_itermax`, `callback_itermax`, `np`
- `pass_n_starts`, `callback_n_starts`
- `pass_optim_maxit`, `callback_optim_maxit`
- `use_deoptim`, `deoptim_parallel`
- `fit_treatment`, `dose_zero_only`, `paired_only`, `truncate_at_treatment`, `ploidy_at_harvest`
- `loss_rescale`, `loss_scale_burden`, `loss_scale_ploidy`, `loss_scale_eps`
- `use_soft_prior`, `lambda_prior`
- `scenario_agg`, `scenario_agg_burden`, `scenario_agg_ploidy`, `scenario_agg_huber_k`
- `scenario_weight_burden`, `scenario_weight_ploidy`
- `prior_center_log10_k_o`, `prior_sd_log10_k_o`
- `prior_center_log10_K_O2`, `prior_sd_log10_K_O2`
- `prior_center_beta_size`, `prior_sd_beta_size`
- `prior_center_log10_n_exp`, `prior_sd_log10_n_exp`
- `prior_center_log10_rho_2N`, `prior_sd_log10_rho_2N`
- `rho_2N_min`, `rho_2N_max`, `c_vol_2N_mm3` (legacy warm-start compatibility), `burden_log_eps`, `huber_k_burden_log`, `burden_exclude_day0`
- `w_burden_chain`, `w_ploidy_chain`, `callback_w_burden`, `callback_w_ploidy`
- `callback_init_pass`, `callback_auto_rescale`
- `auto_tune_iters`
- `resume_from_pass`, `resume_init_tsv_template`, `resume_skip_existing`
- `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`

## Resume Options

This runner supports restarting from a later chain pass.

- `resume_from_pass`: 1-based pass index to start from (`1` means start from the first chain pass)
- `resume_init_tsv_template`: optional warm-start file template path; supports `{seed}` or `__SEED__`
- `resume_skip_existing`: if `TRUE`, existing step/callback outputs are reused and skipped

Examples:

```bash
# Resume from pass 6 using outputs already in the same run directory
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh \
  --run_prefix=my_run \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --seeds_csv=1 \
  --resume_from_pass=6 \
  --resume_skip_existing=TRUE
```

```bash
# Resume from another run root by providing a seed-specific warm-start template
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh \
  --run_prefix=my_new_run \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --seeds_csv=1,2,3 \
  --resume_from_pass=4 \
  --resume_init_tsv_template=/share/.../previous_run_seed{seed}/step03_wb0p6_wp0p4/fit_parameter_stages.tsv
```

## Environment Variables (Alternative)

The same settings can also be provided as environment variables:

```bash
export OUT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results"
export RUN_PREFIX="fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain"
export DATA_DIR="/share/lab_crd/lab_crd/taoli/Project/miningcloneid/data/InVivoData_Gemcitabine"
export SEEDS_CSV="1,2,3"
export K="1e12"
export N_CORES=24

export FIT_TREATMENT=FALSE
export DOSE_ZERO_ONLY=TRUE
export PAIRED_ONLY=TRUE
export TRUNCATE_AT_TREATMENT=FALSE
export PLOIDY_AT_HARVEST=TRUE

export LOSS_RESCALE=TRUE
export LOSS_SCALE_EPS=1e-8
export USE_SOFT_PRIOR=TRUE
export LAMBDA_PRIOR=0.1
export SCENARIO_AGG_BURDEN=huber
export SCENARIO_AGG_PLOIDY=huber
export SCENARIO_AGG_HUBER_K=1.5
export SCENARIO_WEIGHT_BURDEN=TRUE
export SCENARIO_WEIGHT_PLOIDY=TRUE

# Burden observation model (volume-space loss; optional overrides)
export RHO_2N_MIN=3.2e4
export RHO_2N_MAX=5.6e4
# Legacy warm-start conversion baseline (only needed when reading old c_scale-based warm starts)
# export C_VOL_2N_MM3=4.19e-06
export BURDEN_LOG_EPS=1e-12
export HUBER_K_BURDEN_LOG=0.1
# Exclude day=0 burden points from burden loss (recommended TRUE)
export BURDEN_EXCLUDE_DAY0=TRUE

export W_BURDEN_CHAIN="1,0.8,0.6,0.4,0.3,0.25,0.2,0.175,0.15,0.1,0.05,0"
export W_PLOIDY_CHAIN="0,0.2,0.4,0.6,0.7,0.75,0.8,0.825,0.85,0.9,0.95,1"
export CALLBACK_W_BURDEN=1
export CALLBACK_W_PLOIDY=1
export CALLBACK_INIT_PASS=11
export CALLBACK_AUTO_RESCALE=TRUE

export AUTO_TUNE_ITERS=TRUE
export USE_DEOPTIM=FALSE
export DEOPTIM_PARALLEL=FALSE

# Optional manual overrides (used if AUTO_TUNE_ITERS=FALSE or to override auto-tuned values)
export PASS_ITERMAX=220
export CALLBACK_ITERMAX=480
export NP=260
export PASS_N_STARTS=80
export CALLBACK_N_STARTS=140
export PASS_OPTIM_MAXIT=15000
export CALLBACK_OPTIM_MAXIT=28000

# Prevent BLAS/OpenMP oversubscription on shared nodes (runner defaults to 1)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_mono/run_fit_invivo_model_buffering_align_with_Richard_it_stagewise_chain_callback.sh
```

## Auto-Tune Runtime Defaults

If `auto_tune_iters=TRUE` (default), the runner inspects the input data and estimates reasonable defaults for:

- `pass_itermax`, `callback_itermax`, `np`
- `pass_n_starts`, `callback_n_starts`
- `pass_optim_maxit`, `callback_optim_maxit`

If estimation fails, the runner falls back to static defaults.

Static fallback defaults:

- `pass_itermax=150`
- `callback_itermax=280`
- `np=180`
- `pass_n_starts=45`
- `callback_n_starts=90`
- `pass_optim_maxit=9000`
- `callback_optim_maxit=15000`

## Output Layout

For each seed:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step01_wb..._wp.../`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step02_wb..._wp.../`
- `...`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/callback_equal/`

Per-step and callback directories contain the standard fit outputs from `fit_invivo_model_buffering_align_with_Richard.R`, including:

- `best_params.tsv`
- `fit_parameter_stages.tsv` (used as warm-start input for the next pass)
- `fit_summary.tsv`
- `burden_fit.tsv`
- `terminal_ploidy_fit.tsv`
- `fit_config.rds`
- `deoptim_result.rds` (written even when DEoptim is not used; contains optimizer result object)

With the current burden observation model, `burden_fit.tsv` includes volume-space prediction columns such as:

- `pred_burden_volume_mm3`
- `obs_log_burden`
- `pred_log_burden`

Per-step and callback logs are written as sibling files:

- `<step_dir>.log`
- `<callback_dir>.log`

Example:

- `${OUT_ROOT}/${RUN_PREFIX}_seed1/step03_wb0p6_wp0p4.log`
- `${OUT_ROOT}/${RUN_PREFIX}_seed1/callback_equal.log`

Global callback summary across seeds:

- `${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv`

Columns:

- `seed`
- `objective`
- `objective_burden`
- `objective_ploidy`
- `rmse_4N_burden`
- `mean_nll_4N_ploidy`

## How to Choose the Final Seed

Use `${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv` to compare seeds after the callback `(1,1)` pass.

In practice:

- prioritize `objective` (or inspect `objective_burden` / `objective_ploidy` separately if tradeoff matters)
- verify `rmse_4N_burden`
- verify `mean_nll_4N_ploidy`

## Notes

- This runner always calls the fit script with `--two_stage=FALSE`.
- Warm starts are passed via `--init_params_tsv=<previous_pass>/fit_parameter_stages.tsv`.
- Callback defaults to warm-start from `step11` (falls back to nearest available earlier step if missing).
- With `callback_auto_rescale=TRUE` (default), callback run omits manual `loss_scale_*` and lets fit script re-estimate scales from callback init parameters under current loss.
- `paired_only=TRUE` is supported and commonly used so burden/ploidy losses are computed on the same subset of scenarios.
- The fit script now estimates burden observation-model parameters (`rho_2N`, `beta_size`) in addition to the original model parameters.
- If `use_deoptim=TRUE` and `DEoptim` is unavailable in the R environment, the fit script will error.
