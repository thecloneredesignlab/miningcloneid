# Implementation plan: add Hill-gate support to the replicate/legacy live-dead fit path

## Goal

Bring the replicate-level / legacy live-dead fitting path to feature parity with the current primary partial-pooling model for the dFdCTP dose-gating layer.

Specifically, the replicate/legacy path should support the same effective intracellular driver options as the primary path:

1. **beta-only**
   - `use_hill_dose_gate = False`
   - `fit_beta_dose = True` or `False`

2. **Hill-only**
   - `use_hill_dose_gate = True`
   - `fit_beta_dose = False`
   - `fixed_beta_dose = 1.0`

3. **beta + Hill**
   - `use_hill_dose_gate = True`
   - `fit_beta_dose = True`

The key requirement is compatibility, not a redesign of the replicate path.

---

## Current problem

The primary partial-pooling path already supports:

- `beta_dose`
- optional Hill gating
- fixed or fitted Hill parameters
- normalized-objective optimization

But the legacy / replicate-level path still goes through `fit_live_dead_model(...)`, which currently raises:

```python
NotImplementedError("Hill dose gate is currently implemented for fit_joint_partial_pooling_model only.")
```

That means:

- turning replicate fits back on with current defaults fails immediately, because
  - `use_hill_dose_gate = True`
  - `fit_hill_dose_gate = True`

So the current blocker is not the ODE or the likelihood. It is that the legacy fitter has not been wired to pass and, when requested, estimate the Hill-gate parameters.

---

## Design constraints

This implementation should **not** change:

- the dFdCTP PK surface construction
- the ODE state vector
- the Model 1 / 2 / 3 biological structure
- the negative-binomial observation model
- the normalized-objective fix in the primary partial-pooling path
- output file names already used by the main script

The change should be narrowly scoped to:

- replicate/legacy parameter plumbing
- replicate/legacy optimizer setup
- replicate/legacy summaries / diagnostics

---

## Current desired model behavior

The effective dFdCTP signal used by the ODE should remain:

\[
C_{\mathrm{eff}}(t,d)
=
C_{\mathrm{dFdCTP}}(t,d)
\times
\left(\frac{d}{d_{\mathrm{ref}}}\right)^{\beta_{\mathrm{dose}} - 1}
\times
G(d; EC50, h, d_{\mathrm{ref}})
\]

where:

- the beta correction is ploidy-specific when fitted
- the Hill gate uses:
  - `dose_gate_ec50_uM`
  - `dose_gate_hill`
- the Hill gate is shared/global in the primary path today

For the legacy replicate path, we should preserve the same interpretation.

---

## Implementation strategy

### Preferred compatibility target

Make the replicate/legacy path accept the same `JointFitConfig` semantics as the primary partial-pooling path, while keeping its optimizer style intact.

That means:

- if Hill is disabled:
  - legacy path behaves as it does now
- if Hill is enabled but fixed:
  - legacy path uses fixed `dose_gate_ec50_uM` and `dose_gate_hill`
- if Hill is enabled and fitted:
  - legacy path includes Hill parameters in the fitted parameter vector

### Scope boundary

Do **not** attempt to make the replicate path share Hill parameters across all replicates globally in one pooled optimization. That is the primary path’s job.

Instead, for the replicate/legacy fit:

- fixed-Hill mode is straightforward and should be supported
- fitted-Hill mode should fit Hill parameters within that replicate-level optimization, using the same parameter meaning and bounds

This may be statistically weaker than the partial-pooling implementation, but it restores feature compatibility and avoids breaking the current default workflow.

---

## Step 1: audit and isolate the legacy fit entry points

Identify every legacy code path that still calls `fit_live_dead_model(...)` or otherwise uses the old treatment-parameter packing.

At minimum inspect:

- `fit_live_dead_model(...)`
- `fit_joint_one_replicate(...)`
- row-replicate loops in `__main__` / `main(...)`
- any plotting helpers that consume the legacy fit result structure

Expected result:

- one authoritative list of legacy fit entry points
- confirmation of which functions assume the Hill gate is unavailable

---

## Step 2: make treatment-parameter naming fully config-aware in the legacy path

The primary path already uses:

- `get_parameter_names_for_config(...)`
- `get_treatment_parameter_names_for_config(...)`
- `unpack_treatment_params_for_config(...)`

The legacy path should use the same config-aware helpers everywhere possible.

Required changes:

- remove any remaining hard-coded treatment parameter lists in the legacy path
- ensure the legacy path gets:
  - `beta_dose` only if `fit_beta_dose=True`
  - `k_cyto` only for Model 2 / 3
  - `mu_base_death` always

This should eliminate parameter-length mismatches before Hill support is added.

---

## Step 3: add Hill parameters to the legacy optimizer vector when requested

If:

```python
fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate
```

then the legacy fit vector should include:

- `dose_gate_ec50_uM`
- `dose_gate_hill`

These should be parameterized in log space, exactly like the primary path.

Bounds:

- `dose_gate_ec50_uM`: `0.001` to `0.1`
- `dose_gate_hill`: `0.5` to `8.0`

If Hill is enabled but not fitted:

- do not include those coordinates in the optimizer vector
- instead pass:
  - `fit_config.fixed_dose_gate_ec50_uM`
  - `fit_config.fixed_dose_gate_hill`

If Hill is disabled:

- do neither

---

## Step 4: thread Hill parameters through legacy simulation calls

The ODE and simulation layer already support:

- `use_hill_dose_gate`
- `dose_gate_ec50_uM`
- `dose_gate_hill`

The legacy fit path should simply pass them through.

Required updates:

- `fit_live_dead_model(...)` must build the correct simulation arguments
- legacy residual / NLL evaluation must call:
  - `simulate_joint_dfdctp(...)`
  - or `simulate_joint_dfdctp_safe(...)`
  with the same Hill configuration used by the primary path

Remove the current `NotImplementedError` only after this wiring is complete.

---

## Step 5: apply the same fixed-vs-fitted beta behavior in the legacy path

The current primary path supports:

- fitted ploidy-specific beta
- fixed beta via config

The legacy path should mirror that behavior:

- if `fit_beta_dose=True`:
  - include `beta_dose` in the fitted treatment vector
- if `fit_beta_dose=False`:
  - omit it from the fitted vector
  - inject `fit_config.fixed_beta_dose` into simulation

This matters because Hill-only mode depends on `beta_dose` being fixed at `1.0`.

---

## Step 6: update legacy starts conservatively

Do not create a wide new legacy start grid.

Use the same current legacy starts, but:

- add `beta_dose = 1.0` where needed
- add Hill starts only if Hill is fitted

Conservative Hill starts:

- `dose_gate_ec50_uM = 0.0125`
- `dose_gate_hill = 2.0`

If multiple legacy starts are already used, keep them, but append the same default Hill values to each.

The goal is compatibility, not a large search expansion.

---

## Step 7: update legacy summaries and plots

Legacy fit summaries should report the same effective-dose configuration fields as the primary path whenever they are relevant.

Add or ensure reporting of:

- `use_hill_dose_gate`
- `fit_hill_dose_gate`
- `dose_gate_ec50_uM`
- `dose_gate_hill`
- `fit_beta_dose`
- `beta_dose`

If a legacy fit title or annotation already prints fitted parameters, include Hill parameters there too.

If a parameter is fixed rather than fitted, make that explicit in the summary text or fit record.

---

## Step 8: preserve backward compatibility

The safest compatibility contract is:

- old calls with Hill disabled should still work without code changes
- new calls with Hill enabled should no longer raise
- result dictionaries should remain backward-compatible where possible

Suggested approach:

- keep the old keys already returned by legacy fit functions
- add Hill-related keys rather than renaming existing ones

---

## Step 9: verification tests

Add targeted tests before re-enabling replicate fits by default.

### Required tests

1. **Legacy path still works with Hill disabled**
   - existing behavior smoke test

2. **Legacy path works with fixed Hill**
   - `use_hill_dose_gate=True`
   - `fit_hill_dose_gate=False`
   - verify no `NotImplementedError`
   - verify fit result contains fixed Hill values

3. **Legacy path works with fitted Hill**
   - `use_hill_dose_gate=True`
   - `fit_hill_dose_gate=True`
   - verify returned parameter vector/result contains finite Hill values

4. **Hill-only mode works**
   - `fit_beta_dose=False`
   - `fixed_beta_dose=1.0`
   - Hill enabled

5. **Replicate wrapper compatibility**
   - `fit_joint_one_replicate(...)` runs with the current default config
   - no `NotImplementedError`

6. **Main path compatibility**
   - turning row/replicate fitting back on under current defaults does not immediately fail on configuration mismatch

### Optional but useful

7. **Legacy and primary simulation parity**
   - for the same explicit parameter values and same doses, both paths call the same simulation with the same effective-dose settings

---

## Step 10: runtime smoke test

After tests pass, run a small real-data smoke test with:

- current default config
  - Hill enabled
  - Hill fitted
  - beta fitted
  - Model 2
- row/replicate fitting enabled
- one or two replicate rows only if a narrow debug entry point exists

Verify:

- optimizer starts are finite
- no `NotImplementedError`
- no Hill-parameter unpacking mismatch
- no simulation failure caused by Hill parameter passing
- finite fit summaries are produced

This does not need to be a full exhaustive replicate run in the first pass.

---

## Step 11: only then decide whether replicate fits should be re-enabled by default

Do not flip any defaults until the following are true:

1. legacy/replicate Hill support is implemented
2. tests pass
3. at least one real-data replicate smoke fit completes with the current default config

If runtime is still too high, leave replicate fitting opt-in, but it should at least be **compatible** with the current defaults.

---

## Expected deliverable

Implementation is complete when:

1. The legacy replicate fit path no longer raises on Hill-gated configs.
2. The replicate path supports:
   - beta-only
   - Hill-only
   - beta + Hill
3. The same `JointFitConfig` fields control both primary and legacy paths.
4. Legacy fit summaries include Hill settings/values.
5. `fit_joint_one_replicate(...)` works under the current default config.
6. A small real-data replicate smoke test completes without Hill-specific failures.

---

## Non-goals for this pass

Do not do these in the same implementation unless clearly necessary:

- redesign the replicate optimizer
- add replicate-level partial pooling
- share Hill parameters across replicates globally in the legacy path
- change the PK driver or dose surface
- reinterpret the biology of the Hill gate
- change the main production objective or optimizer

Those are separate decisions.
