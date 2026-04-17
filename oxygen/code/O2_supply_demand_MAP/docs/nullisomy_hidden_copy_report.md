# Nullisomy Hidden-Copy Approximation Report

## 1. What Was Implemented

An optional hidden-copy approximation was added to the active `O2_supply_demand_MAP` nullisomy/buffering layer. The model now supports:

- `balanced`
- `dirichlet_multinomial`

The default remains `balanced`, so existing workflows are unchanged unless the new mode is explicitly enabled.

The new Dirichlet mode does not change the coarse state space. It only changes how nullisomy risk is computed inside the loss-branch buffering helper before those risks are cached and reused by the existing missegregation kernel.

## 2. Files Changed

- [model_O2_supply_demand_MAP.cpp](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp)
- [model_O2_supply_demand_MAP.R](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R)
- [fit_invivo_model_O2_supply_demand_MAP.R](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/optimizer/fit_invivo_model_O2_supply_demand_MAP.R)
- [O2_supply_demand.yaml](/Users/4470246/Repositories/miningcloneid/oxygen/config/O2_supply_demand.yaml)
- [test_nullisomy_hidden_copy.R](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/analysis/test_nullisomy_hidden_copy.R)
- [plot_nullisomy_hidden_copy_diagnostics.R](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/analysis/plot_nullisomy_hidden_copy_diagnostics.R)

Key C++ entry points:

- balanced baseline helper: [representative_balanced_copy_vector](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:483)
- exact copy-vector risk helper: [nullisomy_risk_curve_from_copy_vector](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:535)
- Dirichlet-multinomial sampler: [sample_hidden_copy_vector_dirichlet_multinomial](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:587)
- mode-aware cached lookup: [cached_nullisomy_risk_curve](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:639)
- loss-branch lookup: [nullisomy_risk_for_loss](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:718)
- loss survival modifier: [asymmetric_loss_survival_modifier](/Users/4470246/Repositories/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:757)

## 3. Definition of the Dirichlet-Multinomial Approximation

For total chromosome count `N` and `M = N_unit` chromosome classes, the current Dirichlet mode uses a one-copy floor:

- `E = N - M`
- sample excess copies `(e_1,...,e_M) ~ DirichletMultinomial(E; alpha,...,alpha)`
- set `c_j = 1 + e_j`

For each sampled hidden copy vector `c`, the code computes the exact nullisomy risk curve for losing `m` copies using the same subset-count logic already used by the balanced model:

- safe-loss coefficient comes from truncating each chromosome-class polynomial at `c_j - 1`
- nullisomy risk is `1 - safe_prob`

The cached Dirichlet risk is the Monte Carlo average of those exact copy-vector risks over `nullisomy_dirichlet_mc_samples` draws.

## 4. Why the One-Copy Floor Was Used

The one-copy floor was used because viable modeled states are intended to represent cells with at least one copy of each chromosome class before the new missegregation event. That keeps the hidden-copy approximation aligned with the existing buffer interpretation, where nullisomy is created by the current loss event rather than being already present in the incoming hidden state.

The implemented edge-case behavior is:

- if `N < N_unit`, positive-loss nullisomy risk is forced to `1`
- if `N = N_unit`, any positive loss also gives risk `1`

This avoids silently creating invalid negative excess-copy counts.

## 5. Unit Test Results

Unit-test outputs are in:

- [unit_test_summary.tsv](/Users/4470246/Repositories/miningcloneid/code/nullisomy_hidden_copy_tests/unit_test_summary.tsv)
- [unit_test_report.md](/Users/4470246/Repositories/miningcloneid/code/nullisomy_hidden_copy_tests/unit_test_report.md)

Summary:

- balanced regression: `PASS`
- one-copy-floor edge cases: `PASS`
- Dirichlet reproducibility: `PASS`
- cache-key separation: `PASS`
- risk bounds: `PASS`
- monotonicity in loss size: `PASS`
- low-alpha heterogeneity: `PASS`
- Monte Carlo stability: `PASS`
- integration check: `PASS`
- high-alpha approximation to balanced: `WARN`

Important quantitative results:

- balanced mode matched the old exact balanced-copy calculation to numerical precision, with max absolute differences on the order of `1e-15`
- seed-to-seed Monte Carlo variation at `10000` samples was small:
  - `alpha = 1`: worst max absolute difference `0.001891`
  - `alpha = 100`: worst max absolute difference `0.001936`

## 6. Diagnostic Plot Interpretation

Diagnostic outputs are in:

- [nullisomy_risk_curves_hidden_copy.pdf](/Users/4470246/Repositories/miningcloneid/code/nullisomy_hidden_copy_diagnostics/nullisomy_risk_curves_hidden_copy.pdf)
- [copy_variance_hidden_copy.pdf](/Users/4470246/Repositories/miningcloneid/code/nullisomy_hidden_copy_diagnostics/copy_variance_hidden_copy.pdf)
- [diagnostic_note.md](/Users/4470246/Repositories/miningcloneid/code/nullisomy_hidden_copy_diagnostics/diagnostic_note.md)

The main interpretation is:

- `balanced / maximal buffering` is a separate deterministic baseline
- low-`alpha` Dirichlet curves increase nullisomy risk substantially, especially at high `N`
- even very large `alpha` Dirichlet curves remain well above the balanced baseline
- the copy-variance diagnostic explains why: the balanced model has minimal within-class variance, whereas Dirichlet-multinomial hidden states remain heterogeneous even for `alpha = 1000`

Example:

- at `N = 90`, balanced hidden-copy variance is about `0.0866`
- at `N = 90`, Dirichlet `alpha = 1000` mean within-class variance is about `3.11`

That variance gap drives a persistent nullisomy-risk gap.

## 7. Does `alpha = 100` Behave Like Balanced?

No.

Under the current one-copy-floor Dirichlet-multinomial definition, `alpha = 100` does not approximate the balanced deterministic vector. The same is true for `alpha = 1000`.

This is not a coding bug. It is a property of the chosen prior.

Under the implemented prior:

- `E = N - M`
- `(e_1,...,e_M) ~ DirichletMultinomial(E; alpha,...,alpha)`

as `alpha -> infinity`, the excess-copy distribution approaches:

- `Multinomial(E; 1/M,...,1/M)`

That is an equal-probability multinomial-like stochastic limit, not the deterministic balanced allocation. So the balanced model should be treated as a separate maximal-buffering baseline, not as the `alpha -> infinity` limit of this Dirichlet mode.

## 8. Recommended Alpha Values for Sensitivity Analysis

Reasonable sensitivity values for the current definition are:

- `alpha = 1000`
- `alpha = 100`
- `alpha = 10`
- `alpha = 1`
- `alpha = 0.5`

Interpretation:

- `alpha = 1000` and `alpha = 100` probe the equal-probability multinomial-like high-concentration regime
- `alpha = 10` is intermediate
- `alpha = 1` and `alpha = 0.5` emphasize stronger hidden-copy unevenness

These should be interpreted as a family of stochastic hidden-copy priors, not as a path interpolating back to the balanced baseline.

## 9. Limitations

- The Dirichlet mode uses Monte Carlo averaging, so it is approximate rather than closed-form.
- The one-copy-floor prior is only one possible hidden-karyotype mean-field assumption.
- The current Dirichlet family does not interpolate to the deterministic balanced vector.
- No full chromosome-by-chromosome state tracking was added.
- `alpha` is not yet a fitted parameter; it is only a sensitivity setting.
- The current diagnostics establish behavioral differences in nullisomy-risk buffering, but they do not by themselves determine which hidden-copy prior is biologically best.
