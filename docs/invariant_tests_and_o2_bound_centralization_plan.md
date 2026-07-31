# Invariant Tests And O2 Bound Centralization Plan

## Goal

Add explicit biological-accounting tests for the O2 supply-demand model and centralize the `o2_S0_upper_bound` contract across R and C++ so that:

- event semantics for division, missegregation, WGD, and death are testable and unambiguous;
- simulation-level compartment growth obeys intended invariants in simple limiting regimes;
- oxygen upper-bound handling is enforced consistently even when the C++ entry points are called directly.

This plan is intentionally conservative:

- extend the current test suite instead of creating a parallel testing framework;
- keep the current model semantics unless a test demonstrates a mismatch;
- make the R-to-C++ oxygen bound contract explicit rather than relying on wrapper-side pre-clamping.

## Current State

### Existing useful coverage

The current suite already covers several low-level buffering and drop-routing invariants in:

- [oxygen/tests/testthat/test-buffer-missegregation.R](/Users/4470246/Downloads/miningcloneid/oxygen/tests/testthat/test-buffer-missegregation.R:1)
- [oxygen/tests/testthat/helper-oracle.R](/Users/4470246/Downloads/miningcloneid/oxygen/tests/testthat/helper-oracle.R:1)

Specifically, it already checks:

- symmetric daughter survival under buffering;
- per-division offspring accounting conservation for the missegregation kernel;
- the `p_misseg = 0` limit at the daughter-routing level;
- `boundary = "drop"` versus `boundary = "absorb_minmax"` routing differences;
- dead-buffer routing for nonviable daughters.

### Important current gaps

The following are not yet covered explicitly enough:

1. Full generator-level no-death/no-WGD/no-missegregation invariant.
2. Explicit WGD accounting and semantic interpretation.
3. Simulation-level dead-compartment invariants when `mu_hp = 0`.
4. A centralized oxygen upper-bound contract shared by R and C++.

## Workstream 1: Generator-Level No-Death / No-WGD / No-Missegregation Invariant

### Question to lock down

When:

- `p_misseg = 0`
- `p_wgd = 0`
- `mu_hp = 0`
- buffering is irrelevant
- boundary effects are irrelevant

what is the expected live-mass change per `dt` implied by the generator?

The answer must be made explicit in the tests, not inferred informally.

### Likely intended invariant

For a single-state setup with no death and no alternative daughter routing:

- each live state should contribute growth according to the exact meaning of `lam_max` / `lambda_eff`;
- the diagonal/off-diagonal structure of the generator should produce a predictable net live-mass change per time step;
- that expected change should be derived once and asserted in code.

The test should not assume the answer. It should first document the intended semantics:

- whether `lambda` is a parent event rate;
- whether each event produces two daughters and removes one parent;
- therefore whether the net live gain per event is `+1` cell, `+2` cells, or some other encoded quantity.

### Files to inspect and anchor

- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:1)
- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R:1)
- [oxygen/tests/testthat/test-buffer-missegregation.R](/Users/4470246/Downloads/miningcloneid/oxygen/tests/testthat/test-buffer-missegregation.R:1)

### Implementation steps

1. Add a short comment block in the new test section that states the intended generator semantics in plain language.
2. Build a minimal one-state or narrow-grid test case using `cpp_o2simps_build_G_for_o2_triplet`.
3. Set:
   - `p_mis_base = 0`
   - `p_misseg = 0`
   - `p_wgd = 0`
   - `mu_hp = 0`
   - `O2_growth = FALSE` or choose parameters so `lambda_eff` is easy to compute analytically
4. Convert the returned triplet to a sparse matrix.
5. Compute:
   - per-state live outflow;
   - per-state live inflow;
   - net live mass change implied by the generator.
6. Compare the result against an explicit oracle formula derived from the intended division semantics.

### Test shape

Add one or more tests in `test-buffer-missegregation.R`, unless the file is becoming too overloaded. If it grows too much, split into a new file such as:

- `oxygen/tests/testthat/test-event-accounting.R`

Recommended tests:

- `test_that("no-missegregation no-WGD no-death generator has expected live mass gain", { ... })`
- `test_that("no-missegregation no-WGD generator keeps all live growth on self state", { ... })`

### Acceptance criteria

- The expected net live mass change is written explicitly in the test comments and assertions.
- The test passes for both the C++ generator and any R-side helper path that constructs the same object.
- Future changes to generator semantics will fail this test immediately.

## Workstream 2: Dedicated WGD Biological Accounting Test

### Question to lock down

When `p_wgd > 0`, what does WGD mean mathematically in the implementation?

Possible semantics include:

- WGD is an alternative daughter outcome during division;
- WGD is a parent-state transition;
- WGD is an additional branch layered on top of ordinary daughter accounting;
- WGD is a one-off event with separate mass bookkeeping.

The code needs a test that makes the chosen interpretation explicit.

### Why this matters

Without this test, fitted values of:

- `lam_max`
- `p_wgd`
- `p_misseg`
- `mu_hp`

can produce plausible curves while still having ambiguous biological meaning.

### Files to inspect and anchor

- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:1)
- [oxygen/tests/testthat/helper-oracle.R](/Users/4470246/Downloads/miningcloneid/oxygen/tests/testthat/helper-oracle.R:1)

### Implementation steps

1. Read the generator path in C++ and identify exactly where `p_wgd` affects:
   - live daughter routing;
   - state transitions;
   - any dead-buffer or dropped-mass channels.
2. Write a small oracle helper in `helper-oracle.R` describing the intended WGD event accounting for one parent division.
3. Add a narrow-grid test where:
   - `p_mis_base = 0`
   - `p_misseg = 0`
   - `mu_hp = 0`
   - `p_wgd > 0`
   - oxygen growth effects are disabled or simplified
4. Assert:
   - where the WGD offspring mass lands;
   - whether the source state loses parent mass in the expected way;
   - whether total live offspring mass per division matches the intended biology.

### Recommended test cases

1. `p_wgd > 0`, interior state, no boundary interaction.
2. `p_wgd > 0`, high state near the upper boundary, to verify what happens when WGD would leave the grid.
3. `p_wgd = 1` or a very high synthetic value in a controlled generator-only test, if the code path permits this safely, to make the routing visually obvious.

### Recommended tests

- `test_that("WGD routing matches intended biological event accounting", { ... })`
- `test_that("WGD near upper boundary uses the documented routing rule", { ... })`
- `test_that("WGD total live offspring mass matches oracle semantics", { ... })`

### Acceptance criteria

- The test states plainly whether WGD is an alternative daughter fate, a parent transition, or something else.
- The generator behavior and the documented semantics agree.
- If `p_wgd` semantics drift in the future, the test fails.

## Workstream 3: Simulation-Level Dead-Compartment Invariant For `mu_hp = 0`

### Question to lock down

When `mu_hp = 0`, dead-hypoxia mass should not grow. The only dead-mass growth allowed should come from explicitly modeled nonviable daughter routing, such as:

- dropped out-of-grid daughters under `boundary = "drop"`;
- buffered nonviable daughters.

This needs a simulation-level test, not only a generator-level one.

### Why current coverage is insufficient

The existing tests show this indirectly in low-level generator pieces, but they do not prove that the full simulator:

- keeps dead-hypoxia at zero when hypoxia death is turned off;
- only increases dead-buffer through explicit daughter-routing mechanisms.

### Files to inspect and anchor

- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:1478)
- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R:1)
- [oxygen/tests/testthat/test-buffer-missegregation.R](/Users/4470246/Downloads/miningcloneid/oxygen/tests/testthat/test-buffer-missegregation.R:1)

### Implementation steps

1. Build a minimal simulation argument list for `cpp_o2simps_simulate_one`.
2. Use a simple initial state concentrated on one or two ploidy states.
3. Set:
   - `mu_hp = 0`
   - `p_wgd = 0` for the first test
   - `p_mis_base = 0`, `p_misseg = 0` for the strictest null case
4. Assert over all observation steps:
   - `Ntot_dead_hypoxia_obs == 0`
   - `Ntot_dead_total_obs == 0`
   - dead state matrices remain zero
5. Add a second test with:
   - `mu_hp = 0`
   - `p_misseg > 0`
   - `boundary = "drop"`
   - optional buffering parameters
6. Assert:
   - `Ntot_dead_hypoxia_obs == 0`
   - dead-buffer may grow;
   - dead-total equals dead-buffer, because hypoxia death is disabled.

### Recommended tests

- `test_that("mu_hp zero keeps dead-hypoxia compartment identically zero", { ... })`
- `test_that("with mu_hp zero dead growth occurs only through nonviable daughter routing", { ... })`

### Acceptance criteria

- In the strict null case, all dead compartments remain zero.
- In the drop/buffering case, dead-hypoxia remains zero while dead-buffer behaves as expected.
- The assertions use full simulation outputs, not only one-step generator summaries.

## Workstream 4: Centralize `o2_S0_upper_bound` Between R And C++

### Problem statement

The current implementation uses:

- R-side clamping to `[0, o2_S0_upper_bound]`
- C++-side clamping to `[0, 100]`

This is safe enough for current main call paths, but it leaves a contract mismatch:

- direct C++ calls can accept oxygen values above the active model cap;
- wrapper-side normalization is doing semantic work that the low-level implementation does not enforce.

### Desired design

Make `o2_S0_upper_bound` an explicit low-level simulation argument and enforce it in C++ wherever `o2_S0` and `o2_min` are accepted or derived.

That means:

- R continues to normalize inputs early for readability and convenience;
- C++ becomes the final authoritative guardrail for the same bound;
- both layers enforce the same admissible range.

### Files to modify

- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:1)
- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R:1)
- [oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R:1695)
- [oxygen/code/O2_supply_demand_MAP/vis/viz_invivo_model_O2_supply_demand_MAP_results.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/vis/viz_invivo_model_O2_supply_demand_MAP_results.R:238)
- any other R call sites that construct `sim_args` for `cpp_o2simps_simulate_one`

### C++ changes

1. Add `o2_S0_upper_bound` to the exported C++ argument lists for:
   - `cpp_o2simps_o2_window_supply`
   - `cpp_o2simps_simulate_one`
   - any related scalar/vector helpers exposed through Rcpp where `o2_S0` or `o2_min` are interpreted
2. Introduce a helper like:

   - `clamp_o2_pct_to_upper(x, o2_upper)`

   with behavior:

   - finite positive upper bound required;
   - oxygen values clamped to `[0, min(o2_upper, 100)]`.

3. Replace direct `clamp_o2_pct(...)` uses in supply-target handling with the new bound-aware helper where model semantics require the active upper cap.
4. Add explicit assertions or guarded fallback for invalid `o2_S0_upper_bound`.

### R changes

1. Pass `o2_S0_upper_bound` explicitly into every call to:
   - `cpp_o2simps_o2_window_supply`
   - `cpp_o2simps_simulate_one`
2. Keep current R-side normalization, but document that it is now duplicated intentionally as an early check rather than a semantic necessity.
3. Remove any misleading implication that `[0, 100]` alone is the true low-level admissible region for the active model.

### New tests for the oxygen bound contract

Add tests in a new file such as:

- `oxygen/tests/testthat/test-o2-upper-bound-contract.R`

Recommended tests:

1. `test_that("cpp_o2simps_o2_window_supply respects o2_S0_upper_bound", { ... })`
   - pass `o2_S0` above the configured cap;
   - expect outputs to stay at or below the explicit cap, not merely below `100`.

2. `test_that("cpp_o2simps_simulate_one respects o2_S0_upper_bound for O2_state and targets", { ... })`
   - use exaggerated oxygen inputs in a minimal simulation;
   - verify observed/returned oxygen trajectories never exceed the passed cap.

3. `test_that("R and C++ oxygen supply helpers agree on the active upper bound", { ... })`
   - compare `.o2_supply_demand_from_burden(...)` against `cpp_o2simps_o2_window_supply(...)` under a non-default upper cap.

### Acceptance criteria

- Every simulation/objective path passes `o2_S0_upper_bound` explicitly.
- No meaningful C++ oxygen path relies only on the hard-coded `[0, 100]` clamp when an active model cap is smaller.
- R and C++ agree numerically under a non-default cap.

## Recommended Test File Organization

### Keep in existing file

Add to:

- `oxygen/tests/testthat/test-buffer-missegregation.R`

if the new assertions are short extensions of existing generator-level logic.

### Create new files if needed

Prefer creating:

- `oxygen/tests/testthat/test-event-accounting.R`
- `oxygen/tests/testthat/test-o2-upper-bound-contract.R`

if `test-buffer-missegregation.R` becomes too broad.

This split is preferable if:

- generator semantics and simulation semantics need different fixtures;
- oxygen-bound contract tests are logically separate from buffering semantics.

## Suggested Implementation Order

1. Add generator-level no-death/no-WGD/no-missegregation invariant tests.
2. Add WGD accounting tests.
3. Add simulation-level `mu_hp = 0` dead-compartment tests.
4. Introduce explicit `o2_S0_upper_bound` into C++ entry points.
5. Update all R call sites.
6. Add oxygen upper-bound contract tests.
7. Run the targeted test suite and verify no behavior regressions.

## Minimum Validation Commands

After implementation, run at least:

```bash
Rscript oxygen/tests/run_unit_tests.R test-buffer-missegregation
Rscript oxygen/tests/run_unit_tests.R test-o2-upper-bound-contract
Rscript oxygen/tests/run_unit_tests.R test-joint-soft-coupling
Rscript oxygen/tests/run_unit_tests.R test-invitro-defaults
```

If generator or simulation interfaces change, also run any additional local smoke tests that exercise:

- one in vivo seed fit;
- one joint fit report render;
- one visualization path that calls `cpp_o2simps_simulate_one`.

## Deliverable Definition

This work is complete when:

- the biological meaning of division, missegregation, WGD, and hypoxia death is pinned down by executable invariant tests;
- the oxygen upper bound is an explicit contract shared by R and C++;
- there are no remaining important call paths where semantic validity depends on wrapper-only clamping.

## Workstream 5: Increase Useful Documentation Density

### Problem statement

The codebase contains a large amount of documentation text, but much of it is low-value boilerplate. In particular, many comments look like:

- `Function-specific input argument`
- `Object used by downstream model fitting/simulation steps`

These comments increase volume without clarifying:

- the semantic contract of the function;
- the scale of inputs and outputs;
- which layer owns normalization;
- whether a function returns raw model quantities, objective components, or report-ready summaries.

For a fitting codebase with multiple layers:

- natural-scale parameters;
- transformed optimizer parameters;
- scenario rows;
- normalized configs;
- generator outputs;
- simulation summaries;
- objective component bundles;

this lack of precise local contracts is a real maintenance risk.

### Documentation goal

Replace broad low-information comments with short contract-style comments that answer the questions a developer actually needs:

- what are the inputs and what scale are they on?
- what assumptions have already been enforced by the caller?
- what exactly is returned?
- is the function low-level math, simulation orchestration, objective assembly, or report formatting?

The target style is:

```text
Inputs:
- run_params: natural-scale parameters after transform_params()
- scenario: one row/list from scenario table
- cfg: normalized config from normalize_cfg()

Returns:
- burden: predicted burden at observed days
- ploidy: predicted distribution over chromosome-number states
- components: objective terms before weighting
```

### Scope

This workstream is not about documenting every helper uniformly. It is about increasing useful documentation density in the files that matter most for correctness and auditability.

Priority files:

- [oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invivo_backend.R:1)
- [oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invitro_backend.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_invitro_backend.R:1)
- [oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_fit_joint_backend.R:1)
- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R:1)
- [oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp:1)
- [oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R](/Users/4470246/Downloads/miningcloneid/oxygen/code/O2_supply_demand_MAP/report/render_fit_report.R:1)

### Prioritization rule

Only rewrite comments where the function meets at least one of these criteria:

1. The function is a public or semi-public entry point used by multiple files.
2. The function translates between parameter scales or parameter representations.
3. The function assembles objective terms or weighting.
4. The function constructs simulation inputs or interprets simulation outputs.
5. The function encodes semantics that are easy to misunderstand:
   - live vs dead compartments;
   - generator vs simulator outputs;
   - natural scale vs transformed scale;
   - chromosome number vs ploidy mode;
   - joint union bounds vs provenance bounds.

### Documentation style guide for this cleanup

#### Keep comments short

Do not replace boilerplate with long prose blocks. Prefer:

- 4 to 10 lines of concrete contract notes
- no generic “Purpose” text unless it adds real information

#### State parameter scale explicitly

For any function touching model parameters, say whether inputs are:

- natural scale
- transformed optimizer scale
- center/delta joint representation
- decoded run parameters

This is one of the most valuable pieces of local documentation in this codebase.

#### Distinguish simulation layer from objective layer

For functions that return tables or lists, say clearly whether outputs are:

- raw simulator trajectories
- observation-day summaries
- weighted objective terms
- unweighted objective components
- report-ready display tables

#### Document boundary and constraint assumptions

If a function assumes:

- normalized `cfg`
- already-clamped `run_params`
- joint union bound semantics
- no missing columns

state that directly rather than describing the inputs abstractly.

#### Prefer contract comments over “Purpose/Parameters/Returns” templates when the template adds no value

Bad:

```text
Purpose: Internal helper used by the model fitting and simulation pipeline.
Parameters:
- x: Function-specific input argument.
Returns:
- Object used by downstream model fitting/simulation steps.
```

Better:

```text
Inputs:
- par_t: transformed optimizer vector on joint center/delta scale
- ctx: normalized joint context with union-bound metadata

Returns:
- invivo_par: transformed in vivo vector after center/delta reconstruction
- soft_derived: per-parameter reconstruction table used by in vitro mapping and reporting
- feasible: TRUE only if all reconstructed values stay within joint union bounds
```

### Concrete target functions to rewrite first

#### Joint backend

Highest-priority comment rewrites in:

- `build_joint_context()`
- `joint_build_context_specific_transformed_vectors()`
- `joint_objective_components()`
- `joint_apply_warmup_initial_values()`
- `joint_apply_soft_coupling_start_table()`
- `merge_joint_shared_optimizer_bounds()`

Each of these should explicitly state:

- input scale and ownership;
- whether provenance bounds or joint bounds are being used;
- what “feasible” means;
- what outputs are for.

#### In vivo backend

Highest-priority comment rewrites in:

- `simulate_one()`
- `evaluate_objective_components()`
- any parameter decode/encode helpers
- config normalization / parameter-table ingestion helpers

These should state:

- whether `run_params` are natural scale;
- whether `scenario` is one scenario row/list;
- whether outputs are full trajectories or observation-day summaries;
- whether returned objectives are pre- or post-weighting.

#### In vitro backend

Highest-priority comment rewrites in:

- the in vitro objective assembly function(s)
- flow-density handling functions
- chromosome-number / ploidy observation interpretation helpers
- helper functions delegated into `in-vitro-utils`

These should explicitly state:

- whether the function operates in chromosome-number mode, ploidy mode, or both;
- where observation weighting is applied;
- whether likelihood terms are returned separately or combined.

#### Model layer

Highest-priority rewrites in:

- oxygen supply helpers
- `lambda_eff` / `mu_eff` helpers
- C++ exported simulation entry points

These should state:

- acceptable oxygen scale;
- whether the active oxygen upper bound is enforced here;
- whether outputs are rates, probabilities, targets, or state trajectories.

### C++ documentation cleanup strategy

The C++ file currently contains many repeated structured comment blocks. For this cleanup:

1. Keep comments on exported or semantically important functions.
2. Remove or compress comments on trivial helpers where the code is already obvious.
3. For retained comments, replace generic template entries with:
   - exact input interpretation
   - exact output interpretation
   - any clamping/bound assumptions

For example, for an oxygen helper:

```text
Inputs:
- o2_S0: natural-scale oxygen baseline in percent O2
- o2_min: lower oxygen floor in percent O2
- o2_S0_upper_bound: active model cap for oxygen values

Returns:
- oxygen target after burden-dependent reduction, clamped to [0, min(o2_S0_upper_bound, 100)]
```

### R documentation cleanup strategy

For R functions, use short comment blocks immediately above definitions. Prefer:

- `Inputs:`
- `Assumes:`
- `Returns:`

only when needed.

Add `Assumes:` specifically for functions that rely on prior normalization, for example:

```text
Assumes:
- cfg has already been normalized by normalize_cfg()
- run_params are on natural scale
- o2_S0_upper_bound has already been validated as positive
```

### Relationship to the testing workstreams

Documentation cleanup should follow or accompany the invariant tests, not precede them blindly.

Reason:

- the tests force the code semantics to be made explicit;
- the best comments can then describe those tested semantics accurately;
- this reduces the chance of documenting ambiguous behavior incorrectly.

Recommended sequencing:

1. Write or update invariant tests.
2. Confirm intended semantics from passing tests.
3. Rewrite comments for the affected functions to match the tested contract.

### Deliverables for this documentation workstream

1. Replace low-value boilerplate comments in the priority files.
2. Ensure every major fitting/simulation/objective entry point has a useful contract comment.
3. Ensure every major parameter-translation function states the parameter scale explicitly.
4. Ensure oxygen-bound semantics are documented consistently in both R and C++ after centralization.

### Acceptance criteria

This workstream is complete when:

- major public or semi-public functions in the fitting stack have concise contract comments;
- low-information boilerplate has been removed from the most important files;
- comments explicitly distinguish natural scale, transformed scale, and joint center/delta scale where relevant;
- comments on simulation/objective functions say what is returned in operational terms;
- a reviewer can audit the core semantics of fitting and simulation without reverse-engineering every call chain from code alone.
