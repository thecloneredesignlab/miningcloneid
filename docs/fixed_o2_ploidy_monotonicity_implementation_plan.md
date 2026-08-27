# Fixed-O2 Ploidy Monotonicity Implementation Plan

## Goal

Add an analytical workflow to determine, for each fitted in vivo parameter set, whether the fixed-O2 dominant ploidy attractor changes monotonically or non-monotonically as oxygen varies.

The motivating question is whether a given parameter set `theta` implies a simple directional relationship between oxygen and long-term ploidy composition, or whether the model predicts a more complex O2-dependent attractor landscape.

This should use the existing fixed-O2 analytical eigenmode machinery, not stochastic simulation.

## Conceptual Definition

For each fitted parameter set `theta` and fixed oxygen value `O2`, the current code can compute:

```text
dominant_mean_ploidy(theta, O2)
```

by building the fixed-O2 linear generator and extracting the dominant eigenvector. The relevant existing implementation is:

- `oxygen/code/O2_supply_demand_MAP/analysis/fixed_o2/FixO2_invivo.R`
- `fo2_dominant_attractor_one()`
- `generate_fixo2_attractor_mode_table()`

The proposed monotonicity analysis evaluates the function:

```text
O2 -> dominant_mean_ploidy(theta, O2)
```

over a dense oxygen grid and classifies the curve shape.

This is numerical evaluation over O2, but each O2 point is calculated analytically from the model generator. It should therefore be described as an analytical grid evaluation, not empirical simulation.

## Why A Grid-Based Analytical Evaluation Is Preferred

A closed-form symbolic proof of monotonicity is not practical for this model. The dominant ploidy is derived from the dominant eigenvector of an O2-dependent matrix:

```text
dx/dt = M(theta, O2) x
```

where:

```text
M(theta, O2) = G(theta, O2) - diag(mu(theta, O2))
```

The dominant eigenvector can change nonlinearly with O2. Apparent direction changes can also occur near eigenvalue crossings or near-ties, where the dominant eigenvalue has a small spectral gap. A dense O2 grid paired with spectral-gap diagnostics is therefore more robust and easier to audit than symbolic differentiation.

## Scope

Primary scope:

- in vivo best-fit parameter sets
- fixed-O2 attractor dominant mean ploidy
- O2 grid covering the current biologically relevant range, at minimum `0, 0.1, 0.5, 1, 2, 5`
- dense O2 grid for monotonicity classification
- per-seed curve classification
- cross-seed summary of curve classes

Out of scope for the first implementation:

- refitting any model
- stochastic fixed-O2 simulations
- changing how fixed-O2 modes are assigned
- proving monotonicity algebraically

## Recommended O2 Grid

Use two grids:

1. A dense classification grid, for example:

```text
seq(0, 5, by = 0.025)
```

or a log-enriched grid with more resolution near low O2:

```text
c(0, 10^seq(log10(0.01), log10(5), length.out = 250))
```

2. A reporting grid anchored to familiar values:

```text
0, 0.1, 0.5, 1, 2, 5
```

The dense grid should drive classification. The reporting grid should be used for tables and figure labels.

## Computational Feasibility

This analysis is feasible for all 500 best-fit in vivo parameter sets on a laptop-scale machine.

The current fixed-O2 analytical attractor path uses `N_MIN = 22` and `N_MAX = 154`, so each seed/O2 calculation requires an eigen decomposition of an approximately `133 x 133` matrix. A local benchmark of the existing analytical kernel gave an elapsed time of about `0.011` seconds per seed/O2 point after model load/compile.

Approximate runtime estimates:

| Seeds | O2 grid | Evaluations | Single core | 8 workers |
|---:|---:|---:|---:|---:|
| 500 | 101 points, `0.05` step | 50,500 | ~9.5 min | ~1.2 min |
| 500 | 201 points, `0.025` step | 100,500 | ~18.8 min | ~2.4 min |
| 500 | 501 points, `0.01` step | 250,500 | ~46.9 min | ~5.9 min |
| 5,000 | 201 points | 1,005,000 | ~3.1 h | ~23.5 min |
| 10,000 | 201 points | 2,010,000 | ~6.3 h | ~47 min |
| 128,000 | 201 points | 25,728,000 | ~80 h | ~10 h |

Recommended first pass:

- run all 500 best-fit seeds
- use `seq(0, 5, by = 0.025)` as the default dense grid
- enable seed-level parallelization with `n_workers`
- write intermediate curve tables so the run can be resumed or audited

Recommended refinement pass:

- identify seeds with sign changes, small spectral gaps, or ambiguous classifications
- rerun only those seeds on a denser local grid, such as `0.01` spacing or adaptive refinement around sign changes

Do not run this analysis over all 128,000 DE initial parameter vectors as a default. That is technically possible on a cluster, but it is not necessary for the first biological question. If initial-population behavior becomes important later, use a coarse grid or stratified sample first.

## Curve Features To Compute

For each seed:

- `dominant_mean_ploidy` at every dense-grid O2 value
- dominant growth rate at every O2 value
- spectral gap at every O2 value
- finite differences in dominant mean ploidy
- number of sign changes in the finite-difference curve
- maximum ploidy value and O2 at maximum
- minimum ploidy value and O2 at minimum
- ploidy range across O2
- dominant mode label at each reporting-grid O2
- O2 values where dominant mean ploidy crosses key thresholds, such as 1.5 and 2.0

Suggested output columns:

```text
seed_id
O2_pct
dominant_mean_ploidy
dominant_growth_rate
spectral_gap
mode_label
finite_difference_next
local_slope_sign
```

## Monotonicity Classification

Classify each seed using the dense-grid curve.

Recommended classes:

- `monotone_increasing`
- `monotone_decreasing`
- `approximately_flat`
- `u_shaped`
- `inverted_u_shaped`
- `single_transition_increase_then_plateau`
- `single_transition_decrease_then_plateau`
- `complex_nonmonotone`
- `unreliable_small_spectral_gap`

Use a tolerance so small numerical noise does not create false direction changes.

Recommended starting tolerance:

```text
slope_epsilon = max(0.01, 0.02 * ploidy_range)
```

Operational rule:

1. Smooth only if necessary, and always retain the unsmoothed curve in output.
2. Convert finite differences to signs: `+1`, `0`, `-1`.
3. Collapse adjacent repeated signs and remove zeros.
4. Classify based on the sign sequence.

Examples:

- all `+1`: `monotone_increasing`
- all `-1`: `monotone_decreasing`
- no nonzero signs and small range: `approximately_flat`
- `-1, +1`: `u_shaped`
- `+1, -1`: `inverted_u_shaped`
- more than one direction change: `complex_nonmonotone`

If spectral gap is below a threshold over a meaningful part of the O2 range, append or assign an uncertainty flag.

Recommended spectral-gap flags:

- `min_spectral_gap`
- `fraction_o2_gap_below_0p005`
- `fraction_o2_gap_below_0p01`
- `monotonicity_reliability = reliable | caution_small_gap | unreliable_small_gap`

## Spectral-Gap Interpretation

The spectral gap is essential for interpreting apparent non-monotonicity.

If the top two eigenvalues are close, the dominant eigenvector may change rapidly with O2. A sign change in the ploidy curve near a small spectral gap may indicate a near-tie or mode switch rather than a smooth biological transition.

The report should distinguish:

- robust non-monotonicity with adequate spectral gap
- non-monotonicity localized to small-gap regions
- monotonic curves with weak spectral separation

## Optional Local Derivative Analysis

Do not make eigenvector sensitivity derivatives part of the first implementation.

They may be useful later, but they are more fragile because they depend on differentiating an eigenvector of an O2-dependent matrix and handling near-degenerate eigenvalues. Dense analytical grid evaluation is simpler, easier to validate, and sufficient for the current biological question.

## Biological Interpretation Outputs

The report should help answer:

- Does each seed predict increasing, decreasing, or non-monotone ploidy as O2 changes?
- Are Mode1 and Mode2 seeds associated with different curve classes?
- Are low-O2 and high-O2 mode labels part of the same monotone curve or opposite sides of a non-monotone response?
- Do ambiguous cases occur near small spectral gaps?
- Are parameter families that fit equally well biologically different in their O2-ploidy response curve?

Useful summaries:

- count and fraction of seeds in each monotonicity class
- class distribution by reference-O2 mode label
- class distribution by objective rank/family
- median dominant-ploidy curve per class
- representative seed per class
- parameter differences across curve classes

## Figures

Add figures under the fixed-O2 or parameter-landscape result directory.

Recommended figures:

1. All-seed dominant ploidy versus O2 curves, colored by monotonicity class.
2. Median and IQR dominant ploidy curve per monotonicity class.
3. Heatmap of `dominant_mean_ploidy(seed, O2)`, ordered by class and objective.
4. Heatmap of fixed-O2 mode labels across O2.
5. Spectral-gap heatmap across seeds and O2.
6. Scatter of minimum spectral gap versus ploidy range, colored by monotonicity class.
7. Representative curves for each class, including dominant ploidy and spectral gap panels.

For paper-facing figures, avoid overemphasizing individual noisy curves. Use class-level summaries and representative examples.

## Suggested Implementation Location

Add a new analysis script:

```text
oxygen/code/O2_supply_demand_MAP/analysis/fixed_o2/fixo2_ploidy_monotonicity.R
```

or, if the intent is to keep all parameter-landscape downstream analyses together:

```text
oxygen/code/O2_supply_demand_MAP/analysis/parameter_landscape_clustering/fixo2_ploidy_monotonicity.R
```

Preferred first location: `analysis/fixed_o2`, because the calculation is conceptually part of the fixed-O2 attractor machinery.

## Implementation Tasks

### A. Add Curve Table Generation

Read in vivo best-fit parameter sets and compute fixed-O2 attractor outputs over the dense O2 grid using the same helper path as `generate_fixo2_attractor_mode_table()`.

Output:

```text
fixed_o2_ploidy_monotonicity_curves.tsv
```

Each row should be one seed and one O2 value.

### B. Add Curve Classification

Implement a deterministic classifier for each seed's O2-ploidy curve.

Output:

```text
fixed_o2_ploidy_monotonicity_by_seed.tsv
```

Include:

- curve class
- number of sign changes
- ploidy range
- O2 at minimum ploidy
- O2 at maximum ploidy
- spectral-gap summary
- reliability flag
- mode labels at reporting-grid O2 values

### C. Add Cross-Reference Mode Comparison

Join the monotonicity classifications to existing fixed-O2 mode labels at reporting-grid O2 values.

Output:

```text
fixed_o2_ploidy_monotonicity_mode_crosswalk.tsv
```

This should make it easy to ask whether a seed that is Mode1 at one O2 and Mode2 at another O2 is part of a monotone transition or a non-monotone response.

### D. Add Figures

Generate the recommended curve, heatmap, and spectral-gap figures.

All figures should be reproducible from the TSV outputs.

### E. Add HTML Or Markdown Summary

Add a compact report summarizing:

- run arguments
- dense O2 grid
- classification thresholds
- class counts
- class-by-mode tables
- key figures
- reliability caveats

The report should explicitly state that each O2 point is analytically evaluated from the fixed-O2 generator, while the monotonicity class is inferred numerically from the dense O2 grid.

## Validation Plan

Add lightweight validation:

1. Parse the new R script.
2. Run on a small fixture with 3-5 seeds and a short O2 grid.
3. Verify that constant synthetic curves classify as `approximately_flat`.
4. Verify that strictly increasing synthetic curves classify as `monotone_increasing`.
5. Verify that strictly decreasing synthetic curves classify as `monotone_decreasing`.
6. Verify that one-turn synthetic curves classify as `u_shaped` or `inverted_u_shaped`.
7. Verify that low spectral gap triggers the reliability flag.
8. Verify that the reporting-grid O2 values are present in the dense-grid output or are explicitly interpolated/added.
9. Verify that the existing fixed-O2 mode labels at `0, 0.1, 0.5, 1, 2, 5` agree with the corresponding dense-grid curve rows.

## Acceptance Criteria

The implementation is complete when:

- every best-fit seed has a dense fixed-O2 dominant-ploidy curve
- every seed has a monotonicity class and reliability flag
- monotonicity outputs include spectral-gap diagnostics
- fixed-O2 mode labels across reporting-grid O2 values are joined to the curve classification
- summary figures show curve classes and spectral-gap reliability
- validation covers monotone, flat, non-monotone, and small-gap cases

## Interpretation Guidance

Use monotonicity classes to refine fixed-O2 mode interpretation.

A seed can be Mode1 at high O2 and Mode2 at low O2 because it follows a monotone O2 response curve, or because it has a non-monotone attractor landscape. Those are biologically different cases and should not be collapsed into the same interpretation.

If non-monotonicity is robust to spectral-gap filtering, it suggests that different O2 ranges select for different ploidy regimes through different balances of growth, death, missegregation, and buffering. If non-monotonicity appears mainly where the spectral gap is small, interpret it cautiously as a possible near-tie between attractor modes.

The key result should not be only whether the ploidy curve increases or decreases with O2. The more biologically useful output is whether fitted parameter families imply stable, fragile, or multi-regime O2-ploidy responses.
