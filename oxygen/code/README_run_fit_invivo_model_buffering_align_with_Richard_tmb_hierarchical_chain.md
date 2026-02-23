# run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh

This script runs `fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.R` for the TMB hierarchical chain fit of `model_buffering_align_with_Richard.R`.

It supports:

1. Weight-chain fitting across multiple `(w_burden, w_ploidy)` pairs (single run produces all steps)
2. Multi-seed execution (default `1`)
3. TMB hierarchical alternating optimization (`E-step local` + `TMB pooling` + `global shared`)
4. Optional DEoptim for both local and global optimizers
5. Optional loss rescaling for burden/ploidy (`loss_rescale=TRUE` by default)
6. Realtime logs per seed (`run_realtime.log`)

Burden observation model (current default in the fit script):

- burden is fit in volume space (`mm^3`) via a ploidy-aware observation layer
- 2N baseline single-cell volume is controlled by `c_vol_2N_mm3` (runner passthrough optional)
- the fit estimates shared observation-model parameters `c_scale` and `beta_size`
- burden loss is `log-volume Huber`

Default weight chain in this runner:

- `w_burden_chain = 1,0.8,0.6,0.4,0.2,0.175,0.15,0.1,0.05,0`
- `w_ploidy_chain = 0,0.2,0.4,0.6,0.8,0.825,0.85,0.9,0.95,1`

## Run

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh
```

## Pass Parameters Directly (`--key=value`)

You can pass parameters directly without `export`:

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh \
  --run_prefix=fit_stagewise_chain_new_TMB_02212026_001 \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --seeds_csv=1,2,3 \
  --n_cores=32 \
  --fit_treatment=FALSE \
  --paired_only=TRUE \
  --w_burden_chain=1,0.8,0.6,0.4,0.3,0.25,0.2,0.175,0.15,0.1,0.05,0 \
  --w_ploidy_chain=0,0.2,0.4,0.6,0.7,0.75,0.8,0.825,0.85,0.9,0.95,1 \
  --loss_rescale=TRUE \
  --n_alt_iter=3 \
  --use_deoptim_local=TRUE \
  --use_deoptim_global=TRUE \
  --deoptim_itermax_local=220 \
  --deoptim_np_local=260 \
  --deoptim_itermax_global=480 \
  --deoptim_np_global=260 \
  --deoptim_trace=TRUE \
  --c_vol_2N_mm3=4.19e-06 \
  --burden_log_eps=1e-12 \
  --huber_k_burden_log=0.1 \
  --select_rule=min_objective_data
```

Note:

- Old manually tuned `loss_scale_burden` / `loss_scale_ploidy` values from the previous burden-loss definition may not transfer to the current `log-volume Huber` burden loss.
- Recommended first pass: keep `loss_rescale=TRUE` and omit manual `loss_scale_*`, then inspect the logged auto-estimated scales.

## Supported Keys

The runner accepts these `--key=value` options (same names as in the shell script usage):

- `out_root`, `run_prefix`, `data_dir`, `seeds_csv`, `n_cores`, `max_scenarios`
- `fit_treatment`, `dose_zero_only`, `paired_only`, `truncate_at_treatment`, `ploidy_at_harvest`
- `w_burden_chain`, `w_ploidy_chain`
- `loss_rescale`, `loss_scale_burden`, `loss_scale_ploidy`, `loss_scale_eps`
- `c_vol_2N_mm3`, `burden_log_eps`, `huber_k_burden_log`
- `n_alt_iter`, `n_starts_local`, `n_starts_global`, `maxit_local`, `maxit_global`
- `use_deoptim_local`, `use_deoptim_global`
- `deoptim_itermax_local`, `deoptim_np_local`
- `deoptim_itermax_global`, `deoptim_np_global`, `deoptim_trace`
- `lambda_shrink`, `tau_floor`, `tmb_tau_min`, `tmb_log_tau_prior_sd`, `tmb_maxit`, `tmb_rebuild`
- `select_rule`
- `N_MIN`, `N_MAX`, `N_UNIT`, `dt`, `O2`, `K`

## Environment Variables (Alternative)

The same settings can also be provided via environment variables:

```bash
export OUT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results"
export RUN_PREFIX="fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain"
export DATA_DIR="/share/lab_crd/lab_crd/taoli/Project/miningcloneid/data/InVivoData_Gemcitabine"
export SEEDS_CSV="1,2,3"
export N_CORES=32

export FIT_TREATMENT=FALSE
export DOSE_ZERO_ONLY=TRUE
export PAIRED_ONLY=TRUE
export TRUNCATE_AT_TREATMENT=FALSE
export PLOIDY_AT_HARVEST=TRUE

export W_BURDEN_CHAIN="1,0.8,0.6,0.4,0.2,0.175,0.15,0.1,0.05,0"
export W_PLOIDY_CHAIN="0,0.2,0.4,0.6,0.8,0.825,0.85,0.9,0.95,1"

export LOSS_RESCALE=TRUE
export LOSS_SCALE_EPS=1e-8

# Burden observation model (volume-space loss; optional overrides)
export C_VOL_2N_MM3=4.19e-06
export BURDEN_LOG_EPS=1e-12
export HUBER_K_BURDEN_LOG=0.1

export N_ALT_ITER=3
export N_STARTS_LOCAL=6
export N_STARTS_GLOBAL=12
export MAXIT_LOCAL=2500
export MAXIT_GLOBAL=6000

export USE_DEOPTIM_LOCAL=TRUE
export USE_DEOPTIM_GLOBAL=TRUE
export DEOPTIM_ITERMAX_LOCAL=220
export DEOPTIM_NP_LOCAL=260
export DEOPTIM_ITERMAX_GLOBAL=480
export DEOPTIM_NP_GLOBAL=260
export DEOPTIM_TRACE=TRUE

export LAMBDA_SHRINK=1.0
export TAU_FLOOR=0.05
export TMB_TAU_MIN=1e-3
export TMB_LOG_TAU_PRIOR_SD=2.0
export TMB_MAXIT=200
export TMB_REBUILD=FALSE
export SELECT_RULE="min_objective_data"

bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_model_buffering_align_with_Richard_tmb_hierarchical_chain.sh
```

## Output Layout

For each seed, the runner creates:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/run_realtime.log`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/fit_config.rds`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/chain_global_summary.tsv`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/chain_global_summary_ranked.tsv`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/selected_best_step.tsv`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/run_config_extra.tsv`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step01_wb..._wp.../`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step02_wb..._wp.../`
- `...`

Each `stepXX_*` directory contains the main step-level outputs:

- `alt_iter_summary.tsv` (alternating iterations summary)
- `per_sample_theta_i.tsv` (local parameters, natural scale)
- `per_sample_theta_i_transformed.tsv` (local parameters, transformed scale)
- `theta0_robust.tsv` (robustly aggregated parameter center)
- `global_best_params.tsv` (shared/global parameters for that step)
- `per_sample_loss.tsv` (per-scenario losses at final step solution)
- `global_burden_fit.tsv`
- `global_terminal_ploidy_fit.tsv`
- `global_fit_summary.tsv`

With the current burden observation model, `global_burden_fit.tsv` includes volume-space columns such as:

- `pred_burden_volume_mm3`
- `obs_log_burden`
- `pred_log_burden`

## How to Choose the Final Step

This pipeline writes one row per weight step to:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/chain_global_summary.tsv`

and the automatically selected step to:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/selected_best_step.tsv`

By default:

- `select_rule = min_objective_data`

This means the selected step is the one with the minimum data objective (burden + ploidy, after the configured weights and optional loss rescaling).

## Notes on Parallelism / Logs

- The runner writes realtime logs to `run_realtime.log` (via `tee`; `stdbuf` is used when available).
- `n_cores` is passed to the R fit script and affects local/global DEoptim parallel behavior.
- The fit script estimates shared burden observation-model parameters (`c_scale`, `beta_size`) as part of the global/shared parameter vector.
- If DEoptim parallel is requested but unavailable in the R environment (e.g., package missing), the fit script will fail.
- If strict parallel mode is enabled in the fit script and workers cannot be started, the fit will stop instead of silently falling back.
