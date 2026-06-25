# Soft Joint Parameter Coupling Implementation Plan

## Goal

Add an opt-in diagnostic mode for the joint fitter where selected currently shared parameters become softly coupled rather than hard shared.

Target parameterization for each selected shared transformed parameter `c`:

- center parameter: existing transformed shared parameter
- delta parameter: new transformed parameter
- in vivo transformed value: `c + delta / 2`
- in vitro transformed value: `c - delta / 2`

This keeps the implementation on the current optimizer scale. For log-scale parameters this is a symmetric ratio model on the natural scale. For identity-scale parameters it is a symmetric additive split on the optimizer scale.

## Minimal Scope Recommendation

First-stage split set:

- `O2_crit`
- `alpha_o2`
- `mu_hp`
- `p_misseg`

Keep hard-shared for stage 1:

- `lam_max`
- `p_mis_base`
- `n_O`
- `gamma_growth`
- `gamma_mu`
- `k_o_mis`
- `buffer_smax`
- `buffer_beta`
- `buffer_n_exp`
- `p_wgd`

## A. Affected Files And Functions

### Core implementation

`oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R`

Functions to modify or extend:

- `shared_invitro_param_names()`
- `joint_shared_natural_param_names()`
- `merge_joint_shared_optimizer_bounds()`
- `build_joint_context()`
- `build_invitro_transformed_from_joint()`
- `joint_objective_components()`
- `write_joint_outputs()`
- `split_joint_natural_parameter_tables()`

New helper functions to add in this file:

- `joint_soft_split_natural_param_names()`
- `joint_soft_split_transformed_name_map()`
- `joint_soft_split_delta_name()`
- `joint_build_context_specific_transformed_vectors()`
- `joint_soft_coupling_penalty_components()`
- `joint_soft_coupling_summary_table()`

### Reporting

`oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R`

Add:

- loading and display of `joint_soft_coupling.tsv`
- display of the new joint objective components rows

`oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results.R`

Add:

- collection of delta and penalty values into extra-results summaries

`oxygen/code/O2_supply_demand_MAP/analysis/fit_results/extra_results_report.R`

Add plots:

- delta magnitude ranking
- penalty ranking
- vivo vs vitro paired parameter plot
- hard-shared vs soft-coupled objective comparison if baseline is available

### Config and runners

`oxygen/config/O2_supply_demand.yaml`

Add diagnostic config keys.

Potentially:

- `oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_model_O2_supply_demand_MAP.sh`
- `oxygen/code/O2_supply_demand_MAP/runner/run_fit_joint_buffering_local.sh`

Only if CLI pass-through needs to be exposed explicitly.

### Tests

Add a new test file:

- `oxygen/tests/testthat/test-joint-soft-coupling.R`

## B. Naming Convention

### Internal optimizer names

For minimal disruption:

- keep the current shared transformed parameter name as the center parameter
- add a new delta parameter name

Examples:

- center: `log10_O2_crit`
- delta: `delta__log10_O2_crit`

- center: `log10_alpha_o2`
- delta: `delta__log10_alpha_o2`

- center: `log10_mu_hp`
- delta: `delta__log10_mu_hp`

- center: `log10_p_misseg`
- delta: `delta__log10_p_misseg`

This avoids changing `decode_params()` in the in vivo backend.

### Reported natural-scale names

For outputs and report tables:

- `O2_crit_center`
- `O2_crit_vivo`
- `O2_crit_vitro`
- `O2_crit_delta`

and similarly for the rest.

## C. Bounds Handling

### Center bounds

For parameters that remain hard-shared, keep the existing merged shared bounds exactly as today.

For parameters that are soft-coupled, use the overlap of the in vivo and in vitro transformed backend bounds as the center bounds.

Reason:

- when `delta = 0`, both derived context values equal the center,
- the center is feasible for both backends,
- the current hard-shared objective is preserved unless a nonzero delta itself pushes one context to a backend bound.

### Delta bounds

Add bounded deltas per split parameter.

Recommended diagnostic default:

- define delta bounds from the current transformed span:
- `delta_lower = -delta_span_frac * (upper_center - lower_center)`
- `delta_upper = +delta_span_frac * (upper_center - lower_center)`

Recommended default:

- `delta_span_frac = 0.5`

Then clip derived context-specific transformed values back to the corresponding backend bounds before decode/evaluation.

This clipping must be reported, because a clipped vivo/vitro value means the effective split is hitting a bound.

### Why not unconstrained delta

Unconstrained delta will produce:

- boundary pileups,
- wasted optimizer exploration,
- harder interpretation of the penalty.

## D. Default Sigma Choice

Add config keys such as:

- `joint_soft_coupling_enable: TRUE`
- `joint_soft_coupling_params: O2_crit,alpha_o2,mu_hp,p_misseg`
- `joint_soft_coupling_sigma_default: 0.65`
- `joint_soft_coupling_welsch_c: 0.4`

Optional parameter-specific overrides:

- `joint_soft_coupling_sigma_O2_crit`
- `joint_soft_coupling_sigma_alpha_o2`
- `joint_soft_coupling_sigma_mu_hp`
- `joint_soft_coupling_sigma_p_misseg`

Recommended first diagnostic defaults:

- log-scale parameters: `sigma_delta = 0.65`
- identity-scale shape parameters, when eventually used: scale by transformed bound width, for example `0.15 * (upper - lower)`

Interpretation for log-scale parameters:

- `sigma = 0.65` on log10 scale allows roughly a factor of `10^0.65 ~= 4.47` at 1 SD between contexts.

## E. Objective Penalty

Add a new helper in the joint backend:

- `joint_soft_coupling_penalty_components(par_t, ctx)`

It should:

- extract each enabled delta parameter,
- look up its sigma,
- compute the standardized Welsch penalty:

```text
z = delta / sigma
penalty = 0.5 * c^2 * (1 - exp(-(|z| / c)^2))
```

- return per-parameter terms plus total.

Then modify `joint_objective_components()` so that:

- `objective_unpenalized` / `objective_without_soft_or_constraints` is the current weighted in vivo plus in vitro objective before soft-coupling and biological gate penalties
- `objective_soft_coupling` is added separately
- `objective_constraints` remains the current biological gate penalty
- `objective_total` is the sum of all three

Important: the in vivo component here remains the existing `invivo_comp$L`, so it still includes the existing in vivo prior term when soft priors are enabled. Do not silently change it to `invivo_comp$L_data`.

## F. Reporting Changes

### Files to write

From `write_joint_outputs()` add:

- `joint_soft_coupling.tsv`
- additional rows in `joint_components.tsv`

### Suggested `joint_soft_coupling.tsv` columns

- `parameter`
- `split_enabled`
- `center_transformed`
- `delta_transformed`
- `center_natural`
- `vivo_natural`
- `vitro_natural`
- `delta_interpretable`
- `ratio_vivo_to_vitro`
- `regularization_sigma`
- `penalty_paid`
- `center_lower_bound`
- `center_upper_bound`
- `vivo_clipped`
- `vitro_clipped`
- `boundary_status_vivo`
- `boundary_status_vitro`

### Report sections

Add to `render_fit_report.R`:

- a “Soft Coupling Diagnostics” table
- the soft-coupling contribution inside the objective summary

### Extra-results plots

Add to `extra_results_report.R`:

- ranked `|delta|`
- ranked penalty
- paired vivo/vitro values
- objective comparison against hard-shared baseline when present

## G. Unit Tests And Regression Tests

### New unit tests

1. Parameterization round-trip

- with `delta = 0`, derived vivo and vitro transformed values equal center exactly.

2. Symmetry test

- `vivo - vitro = delta` on the transformed scale.

3. Penalty formula

- total penalty equals the sum of the per-parameter standardized Welsch terms.

4. Hard-shared equivalence

- when all enabled deltas are fixed to zero and no context-specific clipping occurs, the joint objective equals the current hard-shared objective to numerical tolerance.

5. Bound clipping behavior

- if a delta pushes vivo or vitro beyond bounds, the derived transformed values are clipped and the summary table flags clipping.

6. Reporting smoke test

- `write_joint_outputs()` produces `joint_soft_coupling.tsv` and new rows in `joint_components.tsv`.

### Regression tests

1. Existing joint fit unchanged when `joint_soft_coupling_enable=FALSE`.

2. Existing report rendering unchanged when no soft-coupling file is present.

3. Existing in vitro and in vivo standalone fits unchanged.

## Proposed Implementation Sequence

1. Add config parsing in `build_joint_context()` for:

- enable flag
- split parameter list
- sigma defaults and overrides
- delta span fraction

2. Add helper functions for:

- mapping natural shared names to transformed center names
- naming delta parameters
- deriving context-specific transformed vectors
- computing penalty terms

3. Extend `build_joint_context()` to:

- append delta parameters to `init`, `lower`, `upper`
- store split metadata in `ctx`

4. Modify `joint_objective_components()` to:

- build separate transformed vectors for in vivo and in vitro
- evaluate existing in vivo/in vitro objectives with those derived vectors
- add the soft-coupling penalty

The implementation must not route soft-coupled values through only `invivo_run_params` and then back-transform in vitro values from that object; doing so would reintroduce hard sharing. Build and apply context-specific transformed vectors explicitly.

5. Extend `write_joint_outputs()` to emit:

- `joint_soft_coupling.tsv`
- expanded `joint_components.tsv`

6. Extend report/analysis code to display the new diagnostics.

Report compatibility note:

- old reports move hard-shared `O2_crit` into `O2_crit_vivo` / `O2_crit_vitro`;
- soft-coupled parameters should be displayed from `joint_soft_coupling.tsv` as their own table, not duplicated through the old hard-shared parameter section.

7. Add unit tests and verify that the default path remains bitwise or near-bitwise unchanged.

## Warnings

### Do not split all candidate parameters at once

The following groups are likely to produce non-identifiability if split simultaneously:

- `O2_crit`, `n_O`, `alpha_o2`, `gamma_growth`
- `mu_hp`, `gamma_mu`
- `p_misseg`, `k_o_mis`

### Do not migrate probability parameters to logit in the first diagnostic implementation

Scientifically, logit is preferable for:

- `p_mis_base`
- `p_misseg`
- `p_wgd`

But migrating base transforms in both backends is a larger refactor. For the first diagnostic version, use the current transformed scale directly and document that this is an approximation.

## Recommended First Deliverable

Implement an opt-in diagnostic mode with only:

- `O2_crit`
- `alpha_o2`
- `mu_hp`
- `p_misseg`

This is the smallest useful experiment that can answer the scientific question:

- which mechanisms really need context-specific values,
- without making the optimizer and interpretation problems much worse than they already are.
