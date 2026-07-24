# Figure 5 drafting notes

## Delivered scope

Figure 5 contains the required five panels in the required order:

- A: direct in-sample fit adequacy for in-vivo burden, in-vivo terminal
  chromosome number, in-vitro passage growth rate, and in-vitro mean
  chromosome number.
- B: all 14 in-vivo/in-vitro soft-coupled parameter ratios for all six approved
  July joint-pair winners.
- C: fitted proliferation functions for the same six winners.
- D: fitted stress-linked missegregation functions for the same six winners.
- E: fitted post-missegregation survival functions for the same six winners.

The single generator is `scripts/make_figure5.R`. It reads existing tabular
exports only and was run through `scripts/agentRrunner.sh`. It does not source
or evaluate the model, read an RDS, optimize, simulate, refit, or write a
derived data table.

## Directive coverage

- FD06: Figure 5 leads with four direct observed-versus-predicted summaries
  spanning both fitted contexts and both burden/growth and karyotype/ploidy
  observations.
- FD10: Every panel uses the six `joint_pair_best` records in
  `selected_results.tsv`. Panels B-E are programmatically checked to contain
  six solutions at every parameter or displayed function coordinate.
- FD11: The order is fit adequacy, full parameter ratios, proliferation,
  stress-linked missegregation, and post-missegregation survival.
- FD12: Panel B retains the full 14-parameter view. It does not substitute a
  compact robustness score.
- FD17: The main fit block is compact and direct. Objective distributions and
  optimizer diagnostics were not added to the main figure.
- FD18: Work is restricted to Figure 5 drafting from approved saved results.
  The 0-5% oxygen grid in C-D is display interpolation within the complete
  saved range, not a new model-evaluation or parameter-analysis grid.

## Selected-winner universe

All six winners use in-vitro anchor seed 10:

- C01/Sc01: in-vivo seed 366, selected joint seed 472, objective 18.8523.
- C01/Sc02: in-vivo seed 290, selected joint seed 155, objective 19.7913.
- C02/Sc01: in-vivo seed 25, selected joint seed 497, objective 18.9705.
- C02/Sc02: in-vivo seed 322, selected joint seed 54, objective 18.8901.
- C03/Sc01: in-vivo seed 138, selected joint seed 122, objective 19.9782.
- C03/Sc02: in-vivo seed 311, selected joint seed 18, objective 19.4145.

The objective range is 18.8523-19.9782. It is not used to select one of these
six as uniquely representative.

## Immutable inputs and transformations

The active generator reads byte-identical copies under
`source_tables/frozen_inputs/F5/`. The upstream selection source is
`oxygen/results/fitting_output_bundle_20260722/selected_results.tsv`.
`source_tables/frozen_input_manifest.csv` records each upstream-to-frozen
mapping, byte count, and SHA-256. Per winner, the generator reads:

- `invivo_burden_fit.tsv` — 120 rows. Panel A retains the 112 finite, post-day-0
  rows and plots the saved within-series normalized observed and predicted
  burden.
- `invivo_terminal_ploidy_fit.tsv` — 1,064 rows. Panel A reduces each of eight
  harvest histograms to observed and predicted count-weighted mean chromosome
  number.
- `invitro_growth_loglik.tsv` — 142 finite observed/predicted passage-growth
  rows used in Panel A.
- `invitro_lineage_summary.tsv` — 175 rows, of which 20 have finite observed
  and predicted mean karyotype and are used in Panel A.
- `joint_soft_coupling.tsv` — 14 rows. Panel B uses all 84 winner-parameter
  rows, verifies that every saved solution is feasible and unprojected, and
  displays `log2(ratio_vivo_to_vitro)` on one common axis.
- `viz/{invivo,invitro}/functional_curve_oxygen_multi_ploidy.tsv` — saved
  oxygen-response curves. Panels C-D retain 2N and 4N and linearly interpolate
  each saved curve to `seq(0, 5, 0.025)` solely to align traces and compute a
  pointwise median. Every source curve spans 0-5% oxygen.
- `viz/{invivo,invitro}/functional_curve_ploidy.tsv` — 133 rows per context.
  Panel E uses the saved chromosome-number coordinate and
  `viability_after_ms` directly, without interpolation.

Pale points in A are all six fitted predictions for an observation; filled
points are the across-winner median. Vertical ranges are shown only for the
sparser terminal-ploidy and karyotype summaries. Thin curves in C-E are
individual selected solutions and thick curves are pointwise medians.

## Prior-code fidelity

Figure 5A is registered as `novel_no_prior`. The identifiable predecessor
registered for F5B-E,
`oxygen/code/O2_supply_demand_MAP/util/o2_supply_demand_map_joint_parameter_plot.R`,
was copied unchanged before adaptation into `prior_code/F5B/` through
`prior_code/F5E/`. The source and all four copies have SHA-256
`f87d32c9cc231b26081aada85d5ec134136a9b9baa3f7746ac297b82df28f074`.

The prior ratio utility represented one selected fit. Panel B preserves its
natural in-vivo/in-vitro ratio and log2-axis semantics while replacing the
single-fit display with all six values, their range, and median. Historical
single-winner function panels were regenerated from the corresponding saved
curves for the common six-winner universe. Focused adaptation records are in
`code_diffs/F5A.diff` through `code_diffs/F5E.diff`.

## Checks and visual review

- The generator asserts exactly six selected winners and the common seed-10
  in-vitro anchor.
- Panel A asserts exactly six predictions per displayed observation.
- Panel B asserts 84 rows, 14 unique parameters, six values per parameter,
  feasibility at the saved solution, and no projection.
- Panels C-E assert six solutions at every displayed context/state coordinate.
- The required runner completed without warnings after the final layout
  revision.
- The composite is rendered natively at the manuscript size rather than
  raster-resized: the recommended PNG is 2,130 × 2,700 pixels at 300 dpi
  (7.1 × 9.0 inches), and the one-page vector PDF is 511 × 648 points.
- All displayed text is at least 7 points at final size. Panel A is a 2 × 2
  small-multiple block, Panel B remains full-width, and Panels C-E form the
  compact bottom row.
- The final 2,130 × 2,700 PNG was inspected at full resolution after the last
  runner pass. Titles, axes, 2N/4N facet labels, separated oxygen-grid ticks,
  legends, all 14 ratio rows and extremes, thin solution traces, thick medians,
  and the bottom guardrail text are visible without clipping or overlap.

## Scientific caveats

- All four Panel A summaries are in-sample fit checks. There is no held-out
  validation, bootstrap, posterior, Hessian/Jacobian uncertainty product, or
  confidence interval.
- The six winners are competitive optimization solutions sharing the same
  in-vitro seed-10 anchor, not independent biological replicates. Ranges,
  points, and traces across them show solution sensitivity only.
- Within-series burden normalization makes temporal pattern fit visible but
  does not assess absolute burden scale; day 0 is omitted from that compact
  comparison.
- Reducing terminal ploidy and in-vitro karyotype to mean chromosome number
  compresses distributional discrepancies.
- Panel B's common scale deliberately retains the extreme `p_mis_base` ratios;
  broad ranges imply solution dependence and must not be read as mechanistic
  certainty.
- Panels C-E are fitted/model-implied functions. Their cross-context
  differences do not establish unique mechanism, causal necessity, or
  parameter identifiability.
- The stored validation verifies selected-winner provenance and objective
  arithmetic, not predictive validation. Reported DEoptim convergence flags
  and unaccepted local refinements do not supply inferential uncertainty.

## Intentionally skipped outputs

- Necrosis is omitted because all six joint exports contain zero finite
  predicted necrosis values.
- No objective-distribution, optimizer-diagnostic, validation-statistic,
  uncertainty-band, ablation, alternate winner, refit, simulation, or new
  analytical-grid output was produced.
- During isolated panel generation, no CSV/TSV/RDS output, TeX change, or
  production-code change was made. Package-level provenance ledgers and the
  review report are assembled separately.

## Blockers

None. All required saved tables were present, and the final PNG and PDF
rendered successfully.
