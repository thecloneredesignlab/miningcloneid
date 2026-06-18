# Softly Coupled Joint Fitting Report

Branch snapshot: `feat-O2-G-resource` at `dc17ea8`

## Executive Assessment

This change is feasible in the current codebase without rewriting the in vivo or in vitro fit backends. The current joint fitter already has a clear seam where shared parameters are:

1. carried in the in vivo transformed optimizer vector,
2. projected into the in vitro transformed vector, and
3. combined into a single joint objective.

That seam lives in `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R`.

The lowest-risk diagnostic implementation is:

- keep the current shared transformed parameter as the center parameter;
- add one extra transformed delta parameter per selected shared parameter;
- derive `theta_vivo` and `theta_vitro` from `(center, delta)` inside the joint backend only;
- add a separate soft-coupling penalty to the joint objective;
- report center, context-specific values, delta, and penalty contribution explicitly.

This is a moderate change, not a small one, but it is localized. It does not require changing the C++ runtime or the standalone in vivo/in vitro optimizers.

## Current Hard-Shared Parameters

The current joint fitter hard-shares the following natural-scale parameters between in vivo and in vitro:

- `lam_max`
- `p_mis_base`
- `p_misseg`
- `k_o_mis`
- `buffer_smax`
- `buffer_beta`
- `buffer_n_exp`
- `alpha_o2`
- `gamma_growth`
- `mu_hp`
- `gamma_mu`
- `O2_crit`
- `n_O`
- `p_wgd`

Source:

- `shared_invitro_param_names()`
- `joint_shared_natural_param_names()`
- `build_invitro_transformed_from_joint()`

## Current Optimizer Scale By Shared Parameter

| parameter | current in vivo transformed name | in vivo transform | current in vitro transformed name | in vitro transform | recommendation for diagnostic split |
| --- | --- | --- | --- | --- | --- |
| `lam_max` | `log10_lam_max` | log10 | `log10_lam_max` | log10 | use current transformed scale |
| `p_mis_base` | `log10_p_mis_base` | log10 | `log10_p_mis_base` | log10 | use current transformed scale first; scientific ideal is logit |
| `p_misseg` | `log10_p_misseg` | log10 | `log10_p_misseg` | log10 | use current transformed scale first; scientific ideal is logit |
| `k_o_mis` | `log10_k_o_mis` | log10 | `log10_k_o_mis` | log10 | use current transformed scale |
| `buffer_smax` | `buffer_smax` | identity | `buffer_smax` | identity | use current transformed scale |
| `buffer_beta` | `log10_buffer_beta` | log10 | `log10_buffer_beta` | log10 | use current transformed scale |
| `buffer_n_exp` | `log10_buffer_n_exp` | log10 | `log10_buffer_n_exp` | log10 | use current transformed scale |
| `alpha_o2` | `log10_alpha_o2` | log10 | `log10_alpha_o2` | log10 | use current transformed scale |
| `gamma_growth` | `gamma_growth` | identity | `gamma_growth` | identity | use current transformed scale |
| `mu_hp` | `log10_mu_hp` | log10 | `log10_mu_hp` | log10 | use current transformed scale |
| `gamma_mu` | `gamma_mu` | identity | `gamma_mu` | identity | use current transformed scale |
| `O2_crit` | `log10_O2_crit` | `log10_nonnegative` | `log10_O2_crit` | log10 | use current transformed scale |
| `n_O` | `n_O` | identity | `n_O` | identity | use current transformed scale |
| `p_wgd` | `log10_p_wgd` | log10 | `log10_p_wgd` | log10 | use current transformed scale first; scientific ideal is logit |

## Parameter Classes

### Positive parameters that are naturally compatible with log-ratio style divergence

These are the cleanest candidates because the current transformed scale is already log-like:

- `lam_max`
- `p_mis_base`
- `p_misseg`
- `k_o_mis`
- `buffer_beta`
- `buffer_n_exp`
- `alpha_o2`
- `mu_hp`
- `O2_crit`
- `p_wgd`

For the diagnostic implementation, the simplest version is to define the center directly on the current transformed scale and use an additive delta on that same scale. For log10-transformed parameters this is equivalent to a ratio penalty on the natural scale up to a constant factor.

### Probabilities that are scientifically better handled on a logit divergence scale

- `p_mis_base`
- `p_misseg`
- `p_wgd`

These are all probabilities bounded in `(0,1)`. A logit difference is the cleanest scientific parameterization, but the current optimizer does not use logit for them. It uses log10. Because this is intended as a diagnostic version, the recommended first pass is:

- keep the existing log10 center parameterization for minimal code change;
- add deltas on the same transformed scale;
- defer a full migration of these probability parameters to `logit01` unless the diagnostic version proves useful.

### Shared parameters with bounded or special handling that should use the current optimizer scale directly

- `buffer_smax`
- `gamma_growth`
- `gamma_mu`
- `n_O`
- `O2_crit` because the in vivo side uses `log10_nonnegative`

Notes:

- `buffer_smax` is unit-interval bounded and currently identity-scaled.
- `gamma_growth`, `gamma_mu`, and `n_O` are positive but currently identity-scaled in both backends.

## Candidate Parameter Recommendations

### Recommended first-stage soft-split subset

Recommended first stage:

- `O2_crit`
- `alpha_o2`
- `mu_hp`
- `p_misseg`

Reasoning:

- these span the main growth, death, and missegregation mechanisms;
- each is easier to interpret biologically than the corresponding shape parameter in the same family;
- this limits the added dimensionality and reduces within-family confounding.

### Parameter-by-parameter recommendation

| parameter | recommendation | transformed scale for center/delta | rationale |
| --- | --- | --- | --- |
| `O2_crit` | split with soft regularization | current `log10_O2_crit` scale | plausible context shift; interpretable; both arms inform it |
| `n_O` | postpone | current `n_O` identity scale | strongly confounded with `O2_crit`, `alpha_o2`, and `gamma_growth` |
| `alpha_o2` | split with soft regularization | current `log10_alpha_o2` scale | plausible context shift in growth penalty strength; useful first-stage candidate |
| `gamma_growth` | keep hard-shared initially | current `gamma_growth` identity scale | shape parameter; likely weakly identified once `O2_crit` and `alpha_o2` are allowed to move |
| `mu_hp` | split with soft regularization | current `log10_mu_hp` scale | likely cleaner than `gamma_mu`; direct context shift in death strength |
| `gamma_mu` | keep hard-shared initially | current `gamma_mu` identity scale | shape parameter confounded with `mu_hp` and ploidy selection |
| `p_misseg` | split with soft regularization | current `log10_p_misseg` scale first | plausible context shift; easier than splitting `k_o_mis` first |
| `k_o_mis` | keep hard-shared initially | current `log10_k_o_mis` scale | strongly confounded with `p_misseg` as amplitude vs oxygen-sensitivity pair |
| `buffer_smax` | postpone | current `buffer_smax` identity scale | bounded near [0,1]; likely boundary-sensitive |
| `buffer_beta` | postpone | current `log10_buffer_beta` scale | buffer family likely internally confounded; postpone as a block |
| `buffer_n_exp` | postpone | current `log10_buffer_n_exp` scale | same reason as `buffer_beta` |
| `p_wgd` | postpone | current `log10_p_wgd` scale first | likely weakly identified and confounded with growth and missegregation effects |

### Parameters that should remain hard-shared initially

- `lam_max`
- `p_mis_base`
- `gamma_growth`
- `gamma_mu`
- `k_o_mis`

Reasons:

- `lam_max`: global timescale parameter; likely to absorb mismatch too easily
- `p_mis_base`: baseline floor parameter; often weakly separated from `p_misseg`
- `gamma_growth`, `gamma_mu`, `k_o_mis`: shape parameters that will be hard to interpret if split before their paired amplitude/threshold parameters

## Key Identifiability Risks

### Growth/resource family

Likely confounding group:

- `O2_crit`
- `n_O`
- `alpha_o2`
- `gamma_growth`

Interpretation:

- `O2_crit` and `n_O` both shape the hypoxia response curve.
- `alpha_o2` and `gamma_growth` both modulate the severity of growth suppression once stress is present.

Recommendation:

- do not split all of these at once;
- start with `O2_crit` and `alpha_o2`;
- keep `n_O` and `gamma_growth` fixed in the hard-shared layer initially.

### Death family

Likely confounding group:

- `mu_hp`
- `gamma_mu`

Interpretation:

- `mu_hp` sets death scale;
- `gamma_mu` controls ploidy dependence of that death pressure.

Recommendation:

- split `mu_hp` first;
- keep `gamma_mu` hard-shared initially.

### Missegregation family

Likely confounding group:

- `p_misseg`
- `k_o_mis`

Interpretation:

- `p_misseg` is the amplitude;
- `k_o_mis` controls how missegregation responds to oxygen stress.

Recommendation:

- split `p_misseg` first;
- keep `k_o_mis` hard-shared initially.

## Where To Introduce The Parameterization

The center/delta parameterization should be introduced in the joint backend only.

Primary file:

- `oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R`

Primary insertion points:

- `shared_invitro_param_names()`
- `joint_shared_natural_param_names()`
- `build_joint_context()`
- `build_invitro_transformed_from_joint()`
- `joint_objective_components()`
- `write_joint_outputs()`
- `split_joint_natural_parameter_tables()`

Recommended design:

- keep the existing shared transformed parameter names as the center variables for backward compatibility;
- add new transformed parameters only for the split set, for example `delta__log10_O2_crit`;
- reconstruct separate in vivo and in vitro transformed vectors before calling the existing decoders/evaluators.

This is preferable to renaming the whole optimizer vector because the existing in vivo decoder already expects the current names.

## Where To Add The Regularization Penalty

Add it in `joint_objective_components()` after:

- the in vivo objective components are evaluated, and
- the in vitro objective components are evaluated,

but before the final `objective` value is returned.

The penalty should be:

- computed per split parameter,
- summed into a `joint_soft_coupling_penalty`,
- kept separate from the current biological hard-constraint penalty `joint_constraint_penalty_total`.

Recommended decomposition:

- `objective_joint_unpenalized_data = w_vivo * L_vivo + w_vitro * L_vitro`
- `objective_joint_soft_coupling = sum_i welsch(delta_i / sigma_i; c)`
- `objective_joint_constraints = joint_constraint_penalty_total`
- `objective_joint_total = objective_joint_unpenalized_data + objective_joint_soft_coupling + objective_joint_constraints`

The Welsch soft-coupling penalty is approximately quadratic near zero and
smoothly saturates for large standardized splits:

```text
z_i = delta_i / sigma_i

welsch(z_i; c) = 0.5 * c^2 * (1 - exp(-(|z_i| / c)^2))
```

## Reporting Additions

### Tables

Add a new TSV and report table such as `joint_soft_coupling.tsv` with columns:

- `parameter`
- `split_enabled`
- `center_transformed`
- `delta_transformed`
- `theta_center`
- `theta_vivo`
- `theta_vitro`
- `delta_interpretable`
- `ratio_vivo_to_vitro`
- `regularization_sigma`
- `penalty_type`
- `welsch_c`
- `welsch_transition_delta`
- `abs_delta_over_sigma`
- `welsch_saturation_fraction`
- `penalty_region`
- `penalty_paid`
- `vivo_lower_bound`
- `vivo_upper_bound`
- `vitro_lower_bound`
- `vitro_upper_bound`
- `boundary_status_vivo`
- `boundary_status_vitro`

For probability-like parameters, replace `ratio_vivo_to_vitro` with:

- `odds_ratio_vivo_to_vitro`

### Joint components table

Extend `joint_components.tsv` to include:

- `objective_joint_soft_coupling`
- `objective_joint_constraints`
- `objective_joint_unpenalized_data`

and ideally one row per split parameter:

- `soft_coupling_penalty_O2_crit`
- `soft_coupling_penalty_alpha_o2`
- etc.

### Plots

Recommended new plots:

- ranked bar plot of `abs(delta_transformed)`
- ranked bar plot of `penalty_paid`
- paired point or dumbbell plot of `theta_vivo` vs `theta_vitro`
- comparison plot of hard-shared vs soft-coupled best objective
- optional per-parameter profile plot of objective improvement against allowed sigma

## Feasibility Conclusion

Feasible: yes.

Recommended scope for the first diagnostic version:

- split only `O2_crit`, `alpha_o2`, `mu_hp`, and `p_misseg`;
- keep all deltas on the current transformed optimizer scale;
- keep the standalone in vivo and in vitro fitters unchanged;
- add the penalty only in the joint backend;
- make the whole feature opt-in via config.

Main warning:

The change will improve flexibility, but it will also increase non-identifiability if too many within-family parameters are split at once. This should be introduced in stages, not as an all-parameter switch.
