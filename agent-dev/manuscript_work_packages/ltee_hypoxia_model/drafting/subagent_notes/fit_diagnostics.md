# Supplementary fit-diagnostics note

## Purpose

Implements FD17 outside the main five figures: direct observed–predicted blocks
remain compact in Figures 3–5, while optimizer-start distributions and
run-level diagnostics are shown here.

## Sources

The active generator reads byte-identical copies under
`source_tables/frozen_inputs/diagnostics/`, plus the frozen Figure 5 selection
table. `source_tables/frozen_input_manifest.csv` records their upstream paths,
byte counts, and SHA-256 values.

- Separate in-vitro and in-vivo
  `run/extra_results/seed_objective_simple.tsv`,
  `seed_summary.tsv`, and `convergence_summary.tsv`.
- All six approved joint-pair
  `run/extra_results/seed_objective_simple.tsv`, `seed_summary.tsv`, and
  `convergence_summary.tsv`.
- `oxygen/results/fitting_output_bundle_20260722/selected_results.tsv`.

The local generator is
`drafting/scripts/make_fit_diagnostics.R`. It reads saved tables only and does
not rerun or refine any optimizer. Counts, selected-run local codes,
convergence totals, and boundary summaries in panel C are computed directly
from those saved tables rather than entered as figure constants.

## Interpretation

Objective differences describe optimization landscapes. Across-start and
across-winner variation describes solution multiplicity/sensitivity; it is not
replicate variation, posterior uncertainty, a confidence interval, or held-out
validation. Optimizer flags and local convergence codes are reported as stored
without treating them as statistical evidence.

Panels A and B use a log10 objective-difference axis. Exact zero differences
cannot appear on that axis and are displayed at a declared floor of `1e-4`.
The six joint winners all have rank 1 and zero objective difference; separate
warm-start-pair facets keep all six open markers visible without displacing
their coordinates.

## Visual QC

The 7.1-inch PNG was inspected after rendering. Selected solutions, axes,
legends, table values, log/floor disclosure, and the uncertainty disclaimer are
visible without clipping.
