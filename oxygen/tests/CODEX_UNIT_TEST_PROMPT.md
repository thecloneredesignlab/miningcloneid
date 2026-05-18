# Codex task: add and run unit tests for `model_O2_supply_demand_MAP`

## Files in scope

- `model_O2_supply_demand_MAP.cpp`
- `model_O2_supply_demand_MAP.R`
- `fit_invivo_model_O2_supply_demand_MAP.R`
- `viz_invivo_model_O2_supply_demand_MAP_results.R`
- `run_fit_invivo_model_O2_supply_demand_MAP.sh`

---

## Primary objective

Create a robust unit-test suite for the missegregation / buffer-layer logic in `model_O2_supply_demand_MAP`, then run the tests and report which tests pass or fail.

The immediate goal is to catch bugs such as:

- signed daughter transitions being mirrored incorrectly,
- self-transition / zero-shift mass being double-counted,
- dropped daughter mass not being conserved correctly,
- boundary handling behaving inconsistently,
- buffer asymmetry being lost during matrix construction.

---

## Test philosophy

Use **behavior-level tests** derived from the model specification.

Do **not** write tests like “function X internally called helper Y”.

Prefer tests that would still be valid after refactoring, such as:

- monotonicity,
- conservation,
- symmetry vs asymmetry,
- no-missegregation limits,
- boundary-mode consistency,
- compartment separation.

Where possible, implement a **small independent oracle** in the test file instead of calling the production helper being tested.

---

## Test framework requirements

Use R `testthat`.

### Expected test layout

Create something along these lines:

- `oxygen/test/test-buffer-missegregation.R`
- optional helper file:
  - `oxygen/test/helper-oracle.R`
- optional lightweight runner:
  - `run_unit_tests.R`

If there is no package structure, set up the tests so they can still run from the repo root with a plain `Rscript` command.

---

## Environment / loading requirements

The code is not a formal R package, so make the tests self-contained.

The tests should:

1. source `model_O2_supply_demand_MAP.R`
2. allow the R file to compile / load `model_O2_supply_demand_MAP.cpp` via `Rcpp::sourceCpp`
3. fail clearly if required packages are unavailable

Assume the test runner is executed from the repo root.

---

## Functions / interfaces that are especially relevant

These exported Rcpp interfaces are the main test targets:

- `cpp_o2simps_pr_delta_vec(N, p, eps_tail = 1e-8, buffer_smax = 1.0, buffer_beta = 0.0, buffer_n_exp = 1.0, N_unit = 22)`
- `cpp_o2simps_build_B_total_triplet(Nmin, Nmax, p_vec, boundary = "drop", eps_tail = 1e-8, buffer_smax = 1.0, buffer_beta = 0.0, buffer_n_exp = 1.0, N_unit = 22)`
- `cpp_o2simps_build_B_WGD_triplet(...)`
- if needed for one-step checks: `cpp_o2simps_build_G_for_o2_triplet(...)` and/or `cpp_o2simps_simulate_one(...)`

Prefer the narrowest interface that is sufficient for each behavioral check.

---

## Required tests

Implement at least the following tests.

### 1. Independent oracle for buffering survival

Write a tiny standalone oracle that does **not** call production survival code.

For a given total chromosome count `N` and missegregation size `m`, compute
`sN = buffer_smax * exp(-buffer_beta * ((2 * N_unit) / N)^buffer_n_exp)` and
survival `sN^m`, clamped to `[0, 1]`.
- for fixed `N`, survival is nonincreasing in `m`
- increasing `buffer_beta` decreases survival for fixed `N` and `m > 0`

### 2. Signed-daughter asymmetry test

This is the highest-priority test.

For a state where loss survival is strictly below 1, the negative branch must be smaller than the positive branch.

Use a small nonzero missegregation probability, no WGD, and compare `+n` vs `-n` for a case like:

- `N = 33`
- `n = 1`
- small nonzero `p`

Use the oracle to compute the expected ratio:

- positive branch weight = base `d_N(+n)`
- negative branch weight = base `d_N(-n) * Sloss(N, n)`

Because the base kernel is symmetric before the loss penalty, the expected ratio should be exactly `Sloss(N, n)`.

Required assertions:

- weight to `N - n` is strictly smaller than weight to `N + n`
- ratio `(negative / positive)` matches `Sloss(N, n)` within tolerance

Run this through the generator / output level, not only the raw survival function.

### 3. Per-division offspring accounting / conservation

For one mother state with fixed `N`, fixed `p`, no WGD:

- compute expected live offspring mass from the oracle,
- compute expected dropped fraction,
- assert:

`live surviving daughters + dropped daughters = 2` per mitosis

Also compare the dropped mass used by the model to the oracle’s dropped mass.

This test should fail if live offspring mass is doubled or if dropped mass bookkeeping is inconsistent.

### 4. Zero-missegregation sanity test

At `p = 0`:

- all offspring remain at the same state,
- self-transition mass equals exactly 2 daughters per mitosis,
- all off-diagonal missegregation transitions are zero,
- dropped mass is zero.

This is meant to catch accidental duplication of the no-shift branch.

### 5. Boundary-rule consistency test

Compare at least two boundary modes:

- `drop`
- `absorb_minmax`

Pick a state near the upper boundary (for example `N = Nmax - 1`) with nonzero chance of leaving the grid.

Required assertions:

- within-grid weights are unchanged between boundary modes,
- only out-of-grid handling differs,
- under `absorb_minmax`, overflow mass is redirected to the boundary state,
- under `drop`, overflow mass is not inserted into the live grid.

If you can test dropped-mass routing as well, do it.

### 6. Loss-only penalty invariance mini-panel

Keep a small fast CI-friendly panel.

Examples:

- `N = 44, n = 1`: equality expected if survival is 1
- `N = 33, n = 1`: strict inequality expected
- `N = 88, n = 3`: equality expected
- `N = 88, n = 4`: likely strict inequality expected once loss becomes possible

For any `(N, n)`:

- if `Sloss(N, n) == 1`, allow equality of `+n` and `-n`
- if `Sloss(N, n) < 1`, require negative branch < positive branch

### 7. Optional but recommended: one-step compartment test

If practical with the current interfaces, run a one-step simulation from a single mother state with:

- one live cell at state `N`,
- no pre-existing dead cells,
- small `dt`,
- fixed oxygen,
- no crowding,
- no clearance,
- no WGD.

Then verify that:

- live compartment change matches surviving offspring allocation minus mother removal,
- hypoxia-dead compartment only reflects hypoxia death,
- buffer-dead compartment only reflects dropped daughter mass.

If this is too heavy for a first pass, document it as a follow-up and still deliver tests 1–6.

---

## Practical implementation details

### Build sparse matrices from triplets

`cpp_o2simps_build_B_total_triplet` returns triplets (`i`, `j`, `x`, `nrow`, `ncol`).

Write a tiny test helper to convert that output into an ordinary dense matrix for small ranges, or a sparse Matrix object using `Matrix::sparseMatrix`.

### Suggested small grids for tests

Use tiny grids whenever possible, for example:

- `Nmin = 30`, `Nmax = 36`
- or `Nmin = 40`, `Nmax = 50`

This keeps tests fast and makes expectations easy to inspect.

### Tolerances

Use clear tolerances, e.g.:

- `1e-12` for exact structural expectations,
- `1e-8` or `1e-6` for floating-point comparisons.

---

## Deliverables

1. Add the test files.
2. Add any minimal helper / runner files needed.
3. Run the tests.
4. Report:
   - number of tests passed,
   - number failed,
   - exact failing assertions,
   - whether the failures support the suspected bugs.
5. If failures occur, do **not** silently fix production code unless explicitly requested.
   - First show the failing tests and explain what they imply.

---

## Execution commands

Use commands of this general form from the repo root:

```bash
Rscript -e "testthat::test_dir('oxygen/test', reporter = 'summary')"
```

If needed, create a dedicated runner:

```bash
Rscript run_unit_tests.R
```

---

## Completion checklist

Do not stop at writing test files. Finish only after all of the following are done:

- [ ] tests created
- [ ] tests run
- [ ] output captured
- [ ] pass/fail summary written
- [ ] likely bug implications summarized

---

## Important constraints

- Keep tests behavior-level and specification-driven.
- Avoid assertions tied only to the current implementation layout.
- Prefer independent-oracle checks whenever feasible.
- Keep the first version lean, readable, and fast.
- Do not rewrite the model yet; this task is only to create and execute tests.
