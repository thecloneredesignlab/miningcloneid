# Joint Soft-Coupling Bound Simplification Plan

## Goal

Simplify the joint soft-coupling implementation so that bounds are defined once and used consistently.

Target behavior:

- `in vivo` bounds: one explicit transformed bound interval per parameter.
- `in vitro` bounds: one explicit transformed bound interval per parameter.
- joint bounds: one explicit transformed bound interval per optimizer variable.
- for shared/soft-coupled parameters, the joint bound should be the union of the `in vivo` and `in vitro` transformed bounds.
- during joint fitting, only the joint union bound is active as an admissibility rule for reconstructed context-specific values.
- no additional clipping of reconstructed context-specific values during objective evaluation.
- no hidden bound reinterpretation through overlap ranges, start-table mutations, or separate clip maps.

The old implementation should be removed, not retained behind compatibility branches or deprecation paths.

## Current Problems

The current code has multiple overlapping meanings of “bounds”:

1. merged shared bounds in `merge_joint_shared_optimizer_bounds()`
2. overlap-based center bounds in `joint_soft_coupling_metadata()`
3. delta bounds derived from `joint_soft_coupling_delta_span_frac`
4. separate `invitro_clip_lower` / `invitro_clip_upper` maps
5. runtime clipping in `joint_build_context_specific_transformed_vectors()`
6. start-table logic that mutates center and delta bounds after the fact

This creates three kinds of ambiguity:

- the optimizer bounds are not the same as the effective bounds seen by the likelihood
- the warm start may not represent the intended separate-fit values because reconstructed context values can be clipped later
- reporting fields such as `vivo_clipped`, `vitro_clipped`, `boundary_status_*`, and center overlap bounds describe implementation artifacts rather than model assumptions

## Desired Model

For each soft-coupled parameter `Q`, define:

- provenance transformed bounds:
  - `z_Q^vivo in [l_Q^vivo, u_Q^vivo]`
  - `z_Q^vitro in [l_Q^vitro, u_Q^vitro]`
- optimizer variables:
  - `c_Q`
  - `delta_Q`
- reconstructed transformed values:
  - `z_Q^vivo = c_Q + delta_Q / 2`
  - `z_Q^vitro = c_Q - delta_Q / 2`

The optimizer should be allowed to search only in a region where both reconstructed values are valid under the joint union bounds, without any runtime clipping.

That means the optimizer bounds must be chosen to make clipping unnecessary.

There are two clean options:

1. Reparameterize directly with `z_Q^vivo` and `z_Q^vitro`.
2. Keep `center + delta`, but derive admissible `center` and `delta` bounds from the joint union bounds and reject any point whose reconstruction leaves those bounds.

Recommendation:

- keep `center + delta`
- remove all runtime clipping
- treat points outside the valid joint union range as infeasible and return a large objective penalty

This preserves the current scientific parameterization while simplifying the implementation.

## Recommended Simplification

### 1. Keep only one source of truth for transformed bounds

Introduce one bound object for shared/soft-coupled parameters that stores:

- `invivo_lower_t`
- `invivo_upper_t`
- `invitro_lower_t`
- `invitro_upper_t`
- `joint_union_lower_t`
- `joint_union_upper_t`

The joint union bounds should be:

- `joint_union_lower_t = min(invivo_lower_t, invitro_lower_t)`
- `joint_union_upper_t = max(invivo_upper_t, invitro_upper_t)`

This object should be computed once and then reused everywhere.

Important interpretation:

- `invivo_*` and `invitro_*` bounds are provenance inputs
- `joint_union_*` bounds are the only active admissibility bounds during joint fitting

### 2. Remove overlap-based center bounds

Delete the logic that sets:

- `center_lower_t = max(invivo_lower_t, invitro_lower_t)`
- `center_upper_t = min(invivo_upper_t, invitro_upper_t)`

This overlap construction is the main reason separate-fit warm starts become non-representable even when each context-specific value is individually valid under the intended joint union bound.

### 3. Remove `invitro_clip_lower` and `invitro_clip_upper`

These maps should not exist.

The in vitro transformed vector should be built directly from the optimizer state and should already be valid if the optimizer point is valid.

### 4. Remove runtime clipping of reconstructed context-specific values

Delete the clipping in:

- `joint_build_context_specific_transformed_vectors()`
- `build_invitro_transformed_from_joint()`

Specifically remove:

- `vivo <- clip(...)`
- `vitro <- clip(...)`
- `pmin(pmax(...))` for shared/soft-coupled transformed in vitro values

The optimizer should never rely on post hoc projection to produce admissible likelihood inputs.

### 5. Replace clipping with explicit feasibility checks

After reconstructing `z_Q^vivo` and `z_Q^vitro`, check:

- `joint_union_lower_t <= z_Q^vivo <= joint_union_upper_t`
- `joint_union_lower_t <= z_Q^vitro <= joint_union_upper_t`

If any soft-coupled parameter violates the joint union bound:

- mark the point infeasible
- return a large objective penalty

This should happen before evaluating either likelihood component.

This makes the effective search region explicit and removes hidden projection behavior.

Do not reject a reconstructed `z_Q^vivo` merely because it lies outside the original `in vivo` standalone bound, and do not reject a reconstructed `z_Q^vitro` merely because it lies outside the original `in vitro` standalone bound. Those original per-context bounds are used only to define the union.

### 6. Stop mutating bounds from the start table

The start table should provide:

- warm-start init values only

It should not widen:

- center bounds
- delta bounds
- any optimizer bound

If a proposed warm-start value is outside the joint union bounds, the code should stop with a clear error.

This makes representability failures explicit and keeps bounds interpretable.

### 7. Remove `joint_soft_coupling_delta_span_frac`

This setting encodes a second bound system on top of the parameter tables.

After the simplification, it should be removed from:

- config
- context construction
- reporting
- tests

If the project wants narrower soft-coupling search regions later, that should be implemented as explicit optimizer bounds or priors, not hidden span heuristics.

## Required Code Changes

### A. `oxygen/code/O2G_supply_demand_MAP/util/o2g_supply_demand_map_fit_joint_backend.R`

This file will carry most of the change.

#### Remove or rewrite these functions

- `joint_soft_coupling_metadata()`
  - remove overlap-based center bounds
  - remove `delta_span_frac`-based delta bounds
  - replace with simple transformed bound metadata plus optimizer variable names

- `merge_joint_shared_optimizer_bounds()`
  - simplify to compute transformed context bounds and joint union bounds once
  - remove `invitro_clip` output

- `joint_apply_soft_coupling_start_table()`
  - remove bound mutation behavior
  - allow init assignment only
  - fail if requested init is out of bounds

- `joint_build_context_specific_transformed_vectors()`
  - remove clipping
  - reconstruct raw values only
  - add explicit feasibility checks
  - return an infeasibility flag and offending parameter names if violated

- `joint_apply_soft_coupling_to_invitro()`
  - keep only direct transformed assignment
  - no clipping logic

- `build_invitro_transformed_from_joint()`
  - remove `clip_lower` / `clip_upper` arguments
  - remove final `pmin(pmax(...))`

#### Remove these context fields

- `invitro_clip_lower`
- `invitro_clip_upper`
- `delta_span_frac`

#### Update objective evaluation

In `joint_objective_components()`:

- detect infeasible reconstructed context-specific values before calling either backend
- return a large penalty immediately for infeasible points
- do not allow clipped surrogate values to enter the likelihood

### B. `oxygen/config/O2G_supply_demand.yaml`

Remove:

- `joint_soft_coupling_delta_span_frac`

Possibly add, only if needed:

- nothing

The goal is simplification, so avoid replacing it with another hidden bound knob.

### C. Reporting outputs

Current reporting fields tied to clipping and overlap should be removed or rewritten.

Likely removals:

- `vivo_clipped`
- `vitro_clipped`
- `boundary_status_vivo`
- `boundary_status_vitro`
- `center_lower_bound`
- `center_upper_bound`
- `center_lower_transformed`
- `center_upper_transformed`

Replace with simpler fields:

- `invivo_lower_transformed`
- `invivo_upper_transformed`
- `invitro_lower_transformed`
- `invitro_upper_transformed`
- `joint_union_lower_transformed`
- `joint_union_upper_transformed`
- `feasible_at_solution`

Files to update:

- `oxygen/code/O2G_supply_demand_MAP/util/o2g_supply_demand_map_fit_joint_backend.R`
- `oxygen/code/O2G_supply_demand_MAP/report/render_fit_report.R`
- any downstream scripts expecting the old columns

### D. Tests

Update or replace tests in:

- `oxygen/tests/testthat/test-joint-soft-coupling.R`

Remove tests that validate:

- clipping behavior
- delta bound widening via start table
- overlap-based center bounds

Add tests for:

1. union bound construction
2. exact warm-start reconstruction when separate-fit values are within the joint union bounds
3. explicit failure when warm-start values are not representable
4. infeasible optimizer points receiving a large penalty rather than being clipped
5. no hidden bound mutation from start tables

## Implementation Sequence

### Step 1. Simplify metadata and bound construction

- rewrite `merge_joint_shared_optimizer_bounds()`
- rewrite `joint_soft_coupling_metadata()`
- remove `invitro_clip_*` fields from `build_joint_context()`

### Step 2. Remove clipping from reconstruction and in vitro vector assembly

- rewrite `joint_build_context_specific_transformed_vectors()`
- rewrite `build_invitro_transformed_from_joint()`
- simplify `joint_apply_soft_coupling_to_invitro()`

### Step 3. Make infeasibility explicit in the objective

- change `joint_objective_components()` so invalid reconstructed context-specific values produce an immediate large penalty

### Step 4. Simplify warm-start logic

- rewrite `joint_apply_soft_coupling_start_table()`
- remove all bound widening behavior
- fail loudly on non-representable warm starts

### Step 5. Remove old config and reporting concepts

- remove `joint_soft_coupling_delta_span_frac`
- remove clipping/overlap fields from TSV/report outputs

### Step 6. Rewrite tests to match the new model

- delete tests that codify clipping or bound mutation
- add tests for explicit representability and infeasibility handling

## Deletions That Should Happen

These behaviors should disappear completely:

- overlap-based center bounds
- `delta_span_frac`
- `invitro_clip_lower`
- `invitro_clip_upper`
- runtime clipping of reconstructed soft-coupled values
- start-table-driven bound expansion
- reporting of clipping status as if clipping were part of the intended model

## Expected Benefits

- one unambiguous meaning of bounds
- warm starts either represent the intended separate-fit values under the joint union bounds or fail explicitly
- no hidden projection between optimizer space and likelihood inputs
- simpler reporting
- easier manuscript description
- easier debugging of representability problems

## Open Design Decision

The one remaining design choice is whether optimizer bounds for `center` and `delta` should be:

1. broad union-based bounds, with infeasible reconstructed points rejected during objective evaluation
2. analytically tightened center/delta bounds that guarantee feasibility by construction under the joint union bounds

Recommendation:

- implement option 1 first

Reason:

- it is simpler
- it avoids deriving and maintaining another layer of coupled bounds
- it removes clipping immediately
- it keeps the code easier to reason about

If later performance shows too much optimizer time spent on infeasible points, option 2 can be considered separately. It should not be mixed into this cleanup unless strictly necessary.

## Summary

The core simplification is:

- define transformed `in vivo` bounds once
- define transformed `in vitro` bounds once
- define joint union bounds once
- use only the joint union bounds as active admissibility bounds during joint fitting
- reconstruct context-specific values without clipping
- reject infeasible points explicitly only when they leave the joint union bounds
- remove all hidden bound mutation and overlap logic

That produces a simpler, cleaner, and manuscript-aligned implementation than the current clipping-based system.
