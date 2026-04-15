# Doxo Fit Implementation Notes

This note records the corrected implementation map for the doxorubicin fitting work after resolving the runtime/source mismatch between the older ploidy-buffer stack and the current oxygen supply-demand workflow.

## Current conclusion

The active invivo fitting implementation in this repository is the `oxygen/code/O2_supply_demand_MAP/` workflow, not `oxygen/code/model_functions.R` and not the older `oxygen/code/fit_invivo_ploidy_buffer.R` stack.

That conclusion is based on:

- current runner and config wiring
- direct runtime tracing of sourced R files and compiled C++ files
- recent git history showing active April 2026 development in `O2_supply_demand_MAP`
- older and much less recent history for `model_functions.R` and `scr/model_functions_ploidy_buffer.R`

There is no checked-in `oxygen/results/` directory in this workspace, so I could not verify the latest successful run from saved output artifacts. The active-path conclusion therefore rests on the code path that the current runner and config would execute now, plus commit history.

## Active runtime code path

### Entry points

Primary shell entrypoint:

- `oxygen/code/O2_supply_demand_MAP/runner/run_fit_invivo_model_O2_supply_demand_MAP.sh`

Primary fitter entrypoint:

- `oxygen/code/O2_supply_demand_MAP/optimizer/fit_invivo_model_O2_supply_demand_MAP.R`

The shell wrapper simply forwards to the fitter in `--mode=run`.

### Config used by the active workflow

Primary config:

- `oxygen/config/O2_supply_demand.yaml`

Relevant defaults currently present in that config include:

- `run_prefix: fit_invivo_model_O2_supply_demand_MAP`
- `fit_treatment: FALSE`
- `dose_zero_only: TRUE`
- `paired_only: TRUE`
- `o2_burden_feedback: TRUE`
- `O2_growth: TRUE`
- `ploidy_O2_death: ploidy_related`
- `start_with: chr_number`

The `start_with: chr_number` setting matters because parts of the code and some historical filenames still use ploidy-oriented terminology even though the active endpoint mode can now work on chromosome-number-like observations.

### Sourced R files at runtime

The active fitter sources:

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_shared.R`
- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_common_semantics.R`
- `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R`

The key runtime source call is in:

- `fit_invivo_model_O2_supply_demand_MAP.R`, where `model_path <- file.path(workflow_root, "model", "model_O2_supply_demand_MAP.R")` and then `source(model_path)`

### Compiled C++ backend used at runtime

`model_O2_supply_demand_MAP.R` compiles and loads:

- `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp`

through `Rcpp::sourceCpp(...)`, with caching in the dedicated workflow cache directory.

The required exported symbols include:

- `cpp_o2simps_pr_delta_vec`
- `cpp_o2simps_loss_survival_nullisomy`
- `cpp_o2simps_build_B_total_triplet`
- `cpp_o2simps_build_B_WGD_triplet`
- `cpp_o2simps_o2_window_supply`
- `cpp_o2simps_build_G_for_o2_triplet`
- `cpp_o2simps_simulate_one`
- `cpp_o2simps_objective_components_map`

## Active function map

### Parameter decoding and transforms

Main decoding function:

- `decode_params(...)` in `fit_invivo_model_O2_supply_demand_MAP.R`

This function decodes transformed optimization parameters into natural-scale model parameters, including:

- growth parameters `lam_min`, `lam_max`, `k_o`
- baseline and amplified missegregation terms `p_mis_base`, `p_misseg`, `k_o_mis`
- nullisomy/buffering term `gamma_loss`
- WGD term `p_wgd`
- oxygen supply-demand parameters `o2_S0`, `kappa_O`, `eta_o2`, `tau_O2`
- oxygen-growth parameters `alpha_o2`, `gamma_growth`
- hypoxia-death parameters `mu_hp`, `gamma_mu`, `O2_crit`, `n_O`
- dead-compartment clearance `k_clear`
- burden scale `sigma_burden`
- optional treatment parameters `alpha`, `gamma`

Shared config normalization and run-parameter normalization happen in:

- `o2_supply_demand_map_common_semantics.R`

This file also normalizes booleans and canonicalizes fields like `p_mis_base`, `ploidy_O2_death`, `start_with`, `tau_O2`, and treatment settings.

### Simulation

Main R simulation wrapper:

- `simulate_one(...)` in `fit_invivo_model_O2_supply_demand_MAP.R`

This function directly calls the compiled backend:

- `cpp_o2simps_simulate_one(...)` in `model_O2_supply_demand_MAP.cpp`

The simulation returns live and dead observables, including:

- `Ntot_live_obs`
- `Ntot_dead_hypoxia_obs`
- `Ntot_dead_buffer_obs`
- `Ntot_dead_total_obs`
- `Ntot_total_obs`
- the same burden outputs in `Vmm3_*`
- `frac_N_live`

### Objective evaluation

Raw objective decomposition:

- `evaluate_objective_components_raw(...)` in `fit_invivo_model_O2_supply_demand_MAP.R`

This function decodes parameters, prepares C++ scenario structures, and calls:

- `cpp_o2simps_objective_components_map(...)` in `model_O2_supply_demand_MAP.cpp`

Final optimizer-facing objective:

- `evaluate_objective_components(...)`
- `evaluate_objective(...)`

The active scalar objective is:

- `L = L_b + L_p + L_prior`

where:

- `L_b` is the burden likelihood contribution
- `L_p` is the endpoint likelihood contribution
- `L_prior = lambda_prior * L_prior_raw` when the soft prior is enabled

`run_optimizer(...)` is the active optimization driver.

## Active model behavior

### Direct hypoxia-linked death `mu_eff`

The active implementation does include direct hypoxia-linked death.

R mirror:

- `.mu_eff_of_O2(...)` in `model_O2_supply_demand_MAP.R`

C++ core:

- `mu_eff_soft_cpp(...)` in `model_O2_supply_demand_MAP.cpp`

The hypoxia weighting is a Hill-type function of oxygen using `O2_crit` and `n_O`. The death mode depends on `ploidy_O2_death`, which currently supports:

- `uniform`
- `diploid_NULL`
- `ploidy_related`

The current config default is `ploidy_related`.

### Current implementation of `p_mis`

In the active workflow, `p_mis` is not computed by direct interpolation in oxygen. It is linked to `mu_eff`.

R mirror:

- `.pmisseg_of_O2(...)` in `model_O2_supply_demand_MAP.R`

The active formula is:

- `p_mis = clamp01(p_mis_base + p_misseg * mu_eff / (mu_eff + k_o_mis))`

Interpretation:

- `p_mis_base` is the baseline per-chromosome-copy missegregation probability
- `p_misseg` is a death-linked missegregation amplification scale
- `k_o_mis` is a half-saturation scale on `mu_eff`

So the currently active model is death-linked missegregation, not direct oxygen-linked interpolation.

### Missegregation kernel

The active kernel lives in the compiled backend, with R accessors in `model_O2_supply_demand_MAP.R`.

Relevant pieces include:

- `cpp_o2simps_pr_delta_vec`
- `cpp_o2simps_build_B_total_triplet`
- `cpp_o2simps_build_B_WGD_triplet`
- `cpp_o2simps_build_G_for_o2_triplet`

The model remains a per-chromosome-copy missegregation model over chromosome-count states.

### Nullisomy-risk buffering

The active workflow already contains explicit nullisomy-risk-based buffering using hidden balanced chromosome-class copy vectors.

Key C++ functions include:

- `representative_balanced_copy_vector(...)`
- `cached_nullisomy_risk_curve(...)`
- `nullisomy_risk_for_loss(...)`
- `asymmetric_loss_survival_modifier(...)`

R wrapper/helper:

- `.loss_survival_nullisomy(...)` in `model_O2_supply_demand_MAP.R`

This is more specific than the older generic `mr_lethality` mechanism. The implementation explicitly constructs a hidden balanced chromosome-class representation and uses it to estimate nullisomy risk under loss events.

### Meaning of dropped mass and `m_drop`

The active code distinguishes two different kinds of dropped mass.

1. Biologically nonviable daughter mass from nullisomy-based loss survival failure.
2. Boundary-drop mass when a live offspring state falls outside the modeled chromosome grid under `boundary = "drop"`.

Inside the kernel, `mass_dropped` corresponds to true nonviable offspring mass associated with failed survival of loss outcomes, not merely numerical tail truncation.

In the active implementation, both the nullisomy nonviable component and the boundary-drop component feed into the dead-buffer bookkeeping.

This is therefore not equivalent to the older ploidy-buffer notes that interpreted dropped mass mainly as numerical truncation.

### Dead compartments

The active simulator maintains explicit dead compartments:

- `v_dead_hypoxia`
- `v_dead_buffer`

These are updated at every simulation step in `cpp_o2simps_simulate_one(...)`.

Interpretation:

- hypoxia death contributes to the hypoxia-dead compartment
- nullisomy/buffer-drop/nonviable daughter mass contributes to the buffer-dead compartment
- both dead compartments clear at rate `k_clear`

Important nuance:

- hypoxia death is a direct death flow from the live state
- buffer-dead inflow is division-linked and therefore scales with the division generator

### Treatment entry point in the active model

Treatment currently enters as a multiplicative suppression of division-linked dynamics.

In `cpp_o2simps_simulate_one(...)`:

- `tx_mult = exp(-alpha * (dose / dose_ref)^gamma)` after treatment start
- `tx_mult` is clipped to `[tx_mult_min, 1]`
- the division-linked update uses `scalar = DT * crowd_mult * tx_mult`

This multiplier affects the generator-driven live dynamics and the division-linked dead-buffer inflow.

It does not replace or absorb direct hypoxia death, which remains a separate flow through `mu_eff`.

So, in the active code, treatment is neither a pure direct-death term nor a pure `p_mis` driver. It is a generator multiplier layered on top of oxygen-linked growth/death and death-linked missegregation.

### Objective computation

The active objective is a MAP likelihood.

In the compiled backend:

- burden is modeled with a log-normal NLL on predicted total burden
- endpoint data are modeled with a continuous mixture likelihood
- endpoint contributions from 2N and 4N cohorts are balanced when both are present

In the R optimizer wrapper:

- the burden and endpoint contributions are combined as `L_b + L_p`
- an optional soft prior penalty is added

This is the current objective that the optimizer minimizes.

## Legacy and parallel implementations

### Older ploidy-buffer stack

Still present and still callable, but not the current primary workflow:

- `oxygen/code/fit_invivo_ploidy_buffer.R`
- `oxygen/code/scr/model_functions_ploidy_buffer.R`
- `oxygen/code/run_single_sim_v2.R`
- `oxygen/code/run_fit_invivo_ploidy_buffer.sh`

These files are part of an older fitting stack. `fit_invivo_ploidy_buffer.R` explicitly sources `scr/model_functions_ploidy_buffer.R`, and `run_single_sim_v2.R` also defaults to that file.

This stack is not the current active target for new doxorubicin work unless backward compatibility with the older ploidy-buffer workflow is also required.

### Older `model_functions.R` stack

Older still:

- `oxygen/code/model_functions.R`
- `oxygen/code/run_optim.R`
- `oxygen/code/viz_results.R`

This layer appears to be the earliest R implementation in the repository and is no longer the right primary modification target for the current invivo fitter.

### Python exploration layer

Separate from the active R/C++ fitter:

- `code/Missegregation_Model.py`
- `code/FakeMCG4.py`

These are still useful reference code for forward simulation ideas, but they are not the active fitting stack.

## Git-history evidence

### `O2_supply_demand_MAP` is current

Recent commits touching the active workflow include:

- `e61103c` on 2026-04-03: `refactor(o2_supply_demand_map): reorganize workflow code and add profiling/calibration analysis suite`
- `6f95d92` on 2026-04-03: `feat(o2_supply_demand_map): add chr_number endpoint mode with backward-compatible outputs`
- `5b7871b` on 2026-04-04: `Support chromosome-number burden scaling in chr_number mode`
- `a274b12` on 2026-04-06: `Remove beta_size and make burden volume independent of ploidy`
- `a955603` on 2026-04-15: `feat: add O2G supply-demand workflow with dynamic glucose and automated reports`

The April 2026 density and recency of edits strongly indicates that this is the current maintained workflow.

### `scr/model_functions_ploidy_buffer.R` is older

Recent history for that stack is concentrated in February to early March 2026:

- `2f04e77` on 2026-02-10
- `1141988` on 2026-02-15
- `4607066` on 2026-02-16
- `4925a07` on 2026-02-23
- `f4f4b44` on 2026-03-02: reorganization of the older O2 pipelines

This indicates an older parallel implementation, not the most recent active one.

### `model_functions.R` is oldest

History for `oxygen/code/model_functions.R` is sparse:

- `b24826a` on 2025-11-17: initial addition
- `98fd66f` on 2026-02-10: minor debug-era modification

That file is therefore the least current of the relevant modeling layers.

## Discrepancies between earlier notes and the active code

The earlier orientation notes were incomplete in several important ways.

1. They treated `oxygen/code/model_functions.R` as the most relevant inspected model file. That is no longer the correct implementation target.
2. They described `p_mis` mainly as direct oxygen interpolation. The active code uses death-linked missegregation through `mu_eff`.
3. They described the absence of explicit `mu_eff`. The active code does include an explicit hypoxia-linked death function.
4. They interpreted dropped mass mainly as numerical truncation. The active code uses explicit biologically nonviable mass from nullisomy buffering, plus separate boundary-drop bookkeeping.
5. They did not reflect the explicit dead-hypoxia and dead-buffer compartments in the current simulator.
6. They understated the role of the compiled C++ backend. In the active workflow, the core simulation and raw objective evaluation both run through the C++ layer.

## Correct modification target for doxorubicin work

If the doxorubicin-nullisomy mode is implemented for the current invivo fitter, the correct primary modification targets are:

- `oxygen/code/O2_supply_demand_MAP/optimizer/fit_invivo_model_O2_supply_demand_MAP.R`
- `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R`
- `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp`
- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_common_semantics.R`
- plus any new data-loading or diagnostics code needed for the drug-response JSON workflow

The older ploidy-buffer and `model_functions.R` stacks should be treated as legacy or reference implementations unless there is an explicit requirement to support them too.
