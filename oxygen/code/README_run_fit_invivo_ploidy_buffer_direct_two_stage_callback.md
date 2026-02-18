# run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh

This script runs `fit_invivo_ploidy_buffer.R` with a stage-wise warm-start chain:

1. Single-stage chained passes (warm-started from previous pass):  
   Default chain: `(1,0) -> (0.8,0.2) -> (0.6,0.4) -> (0.4,0.6) -> (0.2,0.8) -> (0,1)`
2. Final equal-weight callback pass:  
   Default: `(1,1)`
3. Multi-seed execution (default `1,2,3`)
4. Optional loss rescaling to align burden/ploidy magnitudes before weighting  
   Default in this script: `loss_rescale=TRUE`

## Run

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh
```

## Pass Parameters Directly (`--key=value`)

You can pass parameters directly when running the script, without `export`:

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh \
  --seeds_csv=1,2,3 \
  --n_cores=34 \
  --pass_itermax=120 \
  --callback_itermax=220 \
  --pass_n_starts=40 \
  --callback_n_starts=80 \
  --pass_optim_maxit=8000 \
  --callback_optim_maxit=14000 \
  --loss_rescale=TRUE \
  --loss_scale_eps=1e-8 \
  --w_burden_chain=1,0.8,0.6,0.4,0.2,0 \
  --w_ploidy_chain=0,0.2,0.4,0.6,0.8,1 \
  --callback_w_burden=1 \
  --callback_w_ploidy=1 \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --run_prefix=fit_stagewise_chain_3seed
```

Supported keys:

- `out_root`, `run_prefix`, `seeds_csv`, `k`, `n_cores`, `max_scenarios`
- `pass_itermax`, `callback_itermax`, `np`
- `pass_n_starts`, `callback_n_starts`
- `pass_optim_maxit`, `callback_optim_maxit`
- `use_deoptim`, `deoptim_parallel`
- `dose_zero_only`, `truncate_at_treatment`, `ploidy_at_harvest`
- `loss_rescale`, `loss_scale_burden`, `loss_scale_ploidy`, `loss_scale_eps`
- `w_burden_chain`, `w_ploidy_chain`, `callback_w_burden`, `callback_w_ploidy`

## Environment Variables (Alternative)

The same settings can also be provided via environment variables:

```bash
export SEEDS_CSV="1,2,3"
export N_CORES=34
export USE_DEOPTIM=FALSE
export DEOPTIM_PARALLEL=FALSE
export PASS_ITERMAX=120
export CALLBACK_ITERMAX=220
export PASS_N_STARTS=40
export CALLBACK_N_STARTS=80
export PASS_OPTIM_MAXIT=8000
export CALLBACK_OPTIM_MAXIT=14000
export LOSS_RESCALE=TRUE
export LOSS_SCALE_EPS=1e-8
export OUT_ROOT="/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results"
export RUN_PREFIX="fit_stagewise_chain_3seed"

bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/run_fit_invivo_ploidy_buffer_direct_two_stage_callback.sh
```

## Output Layout

For each seed:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step01_...`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step02_...`
- `...`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/callback_equal`

A summary file is generated at:

- `${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv`

Summary columns:

- `objective`
- `objective_burden`
- `objective_ploidy`
- `rmse_4N_burden`
- `mean_nll_4N_ploidy`

Use this table to select the final seed. In practice, prioritize callback `objective_burden` and `objective_ploidy`, and explicitly check `rmse_4N_burden`.
