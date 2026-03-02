# run_fit_invivo_model_O2_invivo_it_stagewise_chain_callback.sh

This runner executes `fit_invivo_model_O2_invivo.R` with iterative stage-wise warm-start chain plus final callback pass, analogous to the Richard runner workflow.

Workflow:

1. Single-stage chain passes (`two_stage=FALSE`), with warm-start from previous pass.
2. Final callback pass (default `(w_burden, w_ploidy) = (1,1)`), initialized from `step11` by default.
3. Multi-seed execution.
4. Optional auto-tuned optimizer iterations from data complexity.
5. Optional resume-from-pass / skip-existing.
6. Callback warm-start fallback: if `step11` is missing, runner falls back to nearest earlier available step.
7. Callback loss-scale mode: `callback_auto_rescale=TRUE` (default) omits manual `loss_scale_*` so fit re-estimates scales from callback init.

Default chain in this runner:

- `w_burden_chain = 1,0.8,0.6,0.4,0.3,0.25,0.2,0.175,0.15,0.1,0.05,0`
- `w_ploidy_chain = 0,0.2,0.4,0.6,0.7,0.75,0.8,0.825,0.85,0.9,0.95,1`
- callback: `(1,1)`
- callback warm-start default: `callback_init_pass=11`
- callback loss-scale mode default: `callback_auto_rescale=TRUE`

## Model Difference (O2 Only)

Karyotype dynamics are kept aligned with the current Richard model.  
Only the burden-to-O2 relation is changed to a windowed angiogenesis shape:

\[
O2_{eff}(N)=\mathrm{clip}\left(
o2_{min}+\frac{O2_{base}-o2_{min}}{1+(N/K_{down})^{h_{down}}}
+A_{ang}\cdot W(\log_{10}(N+\epsilon)),
0,100
\right)
\]

\[
W(x)=\sigma\left(\frac{x-m_{on}}{s_{on}}\right)\cdot\left(1-\sigma\left(\frac{x-m_{off}}{s_{off}}\right)\right), \quad
m_{off}=m_{on}+\Delta m
\]

Interpretation: early decline, mid-burden recovery, late decline.

## O2 Parameters

Fitted O2-window parameters in this model:

- `K_down` (baseline decline scale)
- `h_down` (baseline decline steepness; optimized as `log10_h_down`)
- `A_ang` (angiogenesis recovery amplitude)
- `m_on` (window on center on `log10(N)`)
- `delta_m` (window width proxy; `m_off = m_on + delta_m`)
- `s_on` (window on slope)
- `s_off` (window off slope)

Configured (not fitted by default) O2 constants:

- `O2` (base O2 level in %, passed as `--O2`, default `5.0`)
- `o2_min` (floor, default `0`)
- `h_down_init` (initial value for fitted `h_down`, default `1`)
- `o2_logn_eps` (epsilon for `log10(N+eps)`, default `1.0`)
- `h_down_min`, `h_down_max` (fit bounds for `h_down`; defaults `0.2`, `5.0`)

Optional initialization defaults for O2-window warm start:

- `o2_a_ang_default` (default `25.0`, percent points)
- `o2_m_on_default` (default `9.0`)
- `o2_delta_m_default` (default `1.0`)
- `o2_s_on_default` (default `0.3`)
- `o2_s_off_default` (default `0.3`)

## Fitted Parameter List (Cited)

- With `--fit_treatment=FALSE`, this model fits **18** parameters (`decode_params` non-treatment branch). [D1]
- With `--fit_treatment=TRUE`, it fits **20** parameters (adds `log10_alpha`, `gamma`). [D2]

| Optimizer parameter (transformed scale) | Natural-scale parameter | Fitted when `fit_treatment=FALSE` | Fitted when `fit_treatment=TRUE` | Meaning in the model | Key citation |
| --- | --- | --- | --- | --- | --- |
| `log10_lam_min` | `lam_min` | Yes | Yes | Lower asymptote of oxygen-modulated growth rate | [D1], [M1] |
| `delta_lam` | `lam_max` | Yes | Yes | Growth-rate gap parameter with `lam_max = lam_min + exp(delta_lam)` (thus `lam_max > lam_min`) | [D1], [M1] |
| `log10_k_o` | `k_o` | Yes | Yes | O2 half-saturation scale for growth modulation | [D1], [M1] |
| `log10_p_misseg` | `p_misseg` | Yes | Yes | Baseline missegregation intensity at low oxygen | [D1], [M2] |
| `log10_k_o_mis` | `k_o_mis` | Yes | Yes | O2 scale controlling how fast missegregation decreases with O2 | [D1], [M2] |
| `beta_buffer` | `beta_buffer` | Yes | Yes | Strength of ploidy-dependent buffering in missegregation transition weights | [D1], [F1] |
| `log10_n_exp` | `n_exp` | Yes | Yes | Exponent controlling the ploidy-buffering shape | [D1], [F1] |
| `log10_smax` | `smax` | Yes | Yes | Maximum survival/retention factor used in missegregation transitions | [D1], [F1] |
| `log10_p_wgd` | `p_wgd` | Yes | Yes | Whole-genome-doubling branch probability | [D1], [F1], [M4] |
| `log10_K_down` | `K_down` | Yes | Yes | Tumor burden scale for baseline O2 decline | [D1], [M3] |
| `log10_h_down` | `h_down` | Yes | Yes | Steepness of baseline O2 decline term | [D1], [M3] |
| `A_ang` | `A_ang` | Yes | Yes | Amplitude of angiogenesis recovery window in O2(N) | [D1], [M3] |
| `m_on` | `m_on` | Yes | Yes | Center of window onset in `log10(N)` | [D1], [M3] |
| `log10_delta_m` | `delta_m` | Yes | Yes | Window width proxy (`m_off = m_on + delta_m`) | [D1], [M3] |
| `log10_s_on` | `s_on` | Yes | Yes | Onset slope of angiogenesis window | [D1], [M3] |
| `log10_s_off` | `s_off` | Yes | Yes | Offset slope of angiogenesis window | [D1], [M3] |
| `log10_rho_2N` | `rho_2N` | Yes | Yes | 2N cell density used in burden-volume observation model | [D1], [V1] |
| `beta_size` | `beta_size` | Yes | Yes | Ploidy-size exponent in cell-volume scaling `(P/2)^beta_size` | [D1], [V1] |
| `log10_alpha` | `alpha` | No | Yes | Treatment effect strength in `exp(-alpha * dose_scaled^gamma)` | [D2], [F2] |
| `gamma` | `gamma` | No | Yes | Treatment dose-response exponent | [D2], [F2] |

### Citation Index

- [D1] Non-treatment parameter vector and mapping: `decode_params(..., fit_treatment=FALSE)` in [`fit_invivo_model_O2_invivo.R`](./fit_invivo_model_O2_invivo.R#L291).
- [D2] Treatment-enabled parameter vector and mapping: `decode_params(..., fit_treatment=TRUE)` in [`fit_invivo_model_O2_invivo.R`](./fit_invivo_model_O2_invivo.R#L240).
- [M1] Growth-rate O2 modulation (`growth_lambda`): [`model_O2_invivo.R`](./model_O2_invivo.R#L212).
- [M2] O2-dependent missegregation formula (`.pmisseg_of_O2`): [`model_O2_invivo.R`](./model_O2_invivo.R#L226).
- [M3] O2 window function parameters (`.o2_window_supply_from_burden`): [`model_O2_invivo.R`](./model_O2_invivo.R#L162).
- [M4] WGD probability usage in generator construction (`wgd_prob_vec`): [`model_O2_invivo.R`](./model_O2_invivo.R#L321).
- [V1] Burden observation model volume scaling with `rho_2N` and `beta_size`: [`fit_invivo_model_O2_invivo.R`](./fit_invivo_model_O2_invivo.R#L88).
- [F1] `beta_buffer`, `n_exp`, `smax`, `p_wgd` passed to generator build: [`fit_invivo_model_O2_invivo.R`](./fit_invivo_model_O2_invivo.R#L757).
- [F2] Treatment multiplier definition (`alpha`, `gamma`): [`fit_invivo_model_O2_invivo.R`](./fit_invivo_model_O2_invivo.R#L867).

## Run

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_dynamic/run_fit_invivo_model_O2_invivo_it_stagewise_chain_callback.sh
```

## Example (Direct Args)

```bash
bash /Users/4482173/Documents/GitHub/miningcloneid/oxygen/code/O2_dynamic/run_fit_invivo_model_O2_invivo_it_stagewise_chain_callback.sh \
  --run_prefix=fit_stagewise_chain_O2_invivo \
  --out_root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid/oxygen/results \
  --seeds_csv=1,2,3 \
  --n_cores=63 \
  --use_deoptim=TRUE \
  --deoptim_parallel=TRUE \
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
  --loss_scale_burden=0.73625 \
  --loss_scale_ploidy=23.6334 \
  --pass_itermax=220 \
  --callback_itermax=480 \
  --np=756 \
  --use_soft_prior=TRUE \
  --lambda_prior=0.1 \
  --scenario_agg_burden=huber \
  --scenario_agg_ploidy=huber \
  --scenario_agg_huber_k=1.5 \
  --scenario_weight_burden=TRUE \
  --scenario_weight_ploidy=TRUE \
  --O2=5 \
  --o2_burden_feedback=TRUE \
  --o2_min=0 \
  --h_down_init=1 \
  --h_down_min=0.2 \
  --h_down_max=5 \
  --o2_logn_eps=1 \
  --o2_a_ang_default=25 \
  --o2_m_on_default=9.0 \
  --o2_delta_m_default=1.0 \
  --o2_s_on_default=0.3 \
  --o2_s_off_default=0.3
```

## Supported Keys (Runner)

- `out_root`, `run_prefix`, `data_dir`, `seeds_csv`, `k`, `n_cores`, `max_scenarios`
- `O2`, `o2_burden_feedback`, `o2_min`, `h_down_init`, `h_down_min`, `h_down_max`, `o2_logn_eps`
- `o2_a_ang_default`, `o2_m_on_default`, `o2_delta_m_default`, `o2_s_on_default`, `o2_s_off_default`
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
- `prior_center_log10_K_down`, `prior_sd_log10_K_down`
- `prior_center_log10_h_down`, `prior_sd_log10_h_down`
- `prior_center_beta_size`, `prior_sd_beta_size`
- `prior_center_log10_n_exp`, `prior_sd_log10_n_exp`
- `prior_center_log10_rho_2N`, `prior_sd_log10_rho_2N`
- `burden_exclude_day0`
- `rho_2N_min`, `rho_2N_max`, `burden_log_eps`, `huber_k_burden_log`
- `w_burden_chain`, `w_ploidy_chain`, `callback_w_burden`, `callback_w_ploidy`
- `callback_init_pass`, `callback_auto_rescale`
- `auto_tune_iters`
- `resume_from_pass`, `resume_init_tsv_template`, `resume_skip_existing`

## Outputs

Per-seed directory:

- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/step01_wb..._wp.../`
- `...`
- `${OUT_ROOT}/${RUN_PREFIX}_seed<seed>/callback_equal/`

Important files per step/callback:

- `best_params.tsv`
- `fit_parameter_stages.tsv` (used as warm-start for next pass)
- `fit_summary.tsv`
- `burden_fit.tsv`
- `terminal_ploidy_fit.tsv`
- `fit_config.rds`
- `deoptim_result.rds`

Cross-seed summary:

- `${OUT_ROOT}/${RUN_PREFIX}_callback_metrics.tsv`

## Notes

- This runner always calls fit with `--two_stage=FALSE`.
- Callback warm-start defaults to `step11`; if missing, it falls back to the nearest earlier available step.
- With `callback_auto_rescale=TRUE` (default), callback run omits manual `loss_scale_*` and lets fit re-estimate scales under callback loss settings.
