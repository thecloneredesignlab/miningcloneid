# Figure 3 drafting notes

## Delivered scope

Figure 3 contains exactly four panels:

- A: observed growth rates and seed-10 fitted growth rates over lineage passage.
- B: seed-10 fitted chromosome-state fractions with individual observed karyotypes overlaid.
- C: the saved fitted post-missegregation survival function.
- D: the saved model-implied nonviable-daughter fraction versus per-chromosome missegregation rate.

There is no panel E.

The single generator is `scripts/make_figure3.R`. It was run with
`scripts/agentRrunner.sh` and does not source the model, read an RDS, optimize,
simulate, refit, or write derived data tables.

## Directive coverage

- FD03: Panel A is the growth-rate comparison only. Its x coordinate is lineage
  depth derived from `segment_id`, labeled “Lineage passage”; the saved
  interleaved `passage_index` and calendar time are not used as the plotted x
  coordinate.
- FD04: Panels A and B retain both control and deprived branches, separately
  faceted for 2N and 4N. Control lineage depths are 1–14; deprived depths are
  3–25 (2N) and 3–24 (4N).
- FD05: Panel E, ablations, and refits are omitted.
- FD06: The figure leads with the direct observed-versus-fitted growth-rate
  comparison. Points are observations and lineage-connected lines are fitted
  values from the same saved run.
- FD17: The main fit block is compact. Optimizer and objective diagnostics were
  not inserted into the main figure.
- FD18: Work is restricted to the approved Figure 3 scope and saved seed-10
  results; no TeX or common workflow files were changed.

## Immutable input tables

The active generator reads byte-identical copies under
`source_tables/frozen_inputs/F3/`. Their upstream directory is
`oxygen/results/fitting_output_bundle_20260722/separate_fits/invitro/selected_seeds/seed10/`:

- `invitro_lineage_summary.tsv` — 131 rows; Panel A.
- `invitro_distribution_summary.tsv` — 9,975 rows; Panel B fitted fractions.
- `invitro_observed_kary.tsv` — 220 rows; Panel B observed overlays.
- `viz/functional_curve_ploidy.tsv` — 133 rows; Panel C.
- `viz/functional_curve_oxygen_multi_ploidy.tsv` — 6,928 rows; Panel D.

`source_tables/frozen_input_manifest.csv` records the original and frozen paths,
byte counts, and SHA-256 values, so rendering does not depend on the ignored
local results bundle.

The identifiable prior generator was copied before adaptation into
`prior_code/F3A/` through `prior_code/F3D/`. All four copies have SHA-256
`20e3bd03cb3bdc077a1d6c57866d9ada0a2359eadc6a1acf5171b14a0f4518d1`.
Focused per-panel adaptation records are in `code_diffs/F3A.diff` through
`code_diffs/F3D.diff`.

## Checks and visual review

- The generator asserts all four 2N/4N × control/deprived strata in Panels A and
  B and all eight 1.5N–5N reference states in Panel D.
- Finite observed growth counts are 12, 12, 46, and 44 for 2N control, 4N
  control, 2N deprived, and 4N deprived, respectively.
- Observed karyotype counts in those strata are 40, 36, 74, and 70.
- As a descriptive check only, the row-level growth-rate RMSE is 0.204 day⁻¹.
  It is not presented as a validation statistic because controls were fitted.
- The saved survival values at the 2N and 4N reference chromosome counts
  (`N = 44` and `N = 88`) are 0.204 and 0.837.
- Final PNG and PDF are 7.1 × 8.8 inches; the PNG is 2,130 × 2,640 pixels at
  300 dpi with a white background. The PDF is predominantly vector; the
  heatmap gradient is embedded as a raster strip.
- The final PNG was inspected at full resolution. An initial collision between
  panel tags and titles was corrected; the revised tags, axes, facet labels,
  heatmap overlay points, and curve legend are legible.

## Scientific caveats

- Seed 10 is the selected rank-1 run among 500 starts, but close alternative
  starts, local convergence code 1, five parameters at bounds, and missing
  Hessian/Jacobian uncertainty diagnostics preclude strong confidence or
  identifiability claims. No uncertainty bands are available.
- Control branches are shorter and were included in fitting; they are not
  held-out validation data.
- “Lineage passage” is a branch-depth coordinate, not elapsed time.
- Shared saved segment predictions are collapsed only after checking that
  replicated rows agree. All finite observed growth measurements remain as
  points.
- Panel B retains the full saved chromosome-state support. Horizontal offsets
  for observed cells are deterministic display jitter only.
- Panels C and D are fitted/model-implied functions, not independent
  measurements and not evidence that either mechanism is uniquely necessary.
- Panel D uses `misseg_nonviable_daughter_fraction`; it does not include
  boundary-dropped state loss.

## Intentionally skipped outputs

- No Figure 3E.
- No alternate seed, refit, ablation, parameter grid, or uncertainty product.
- No optimizer/objective panel in the main figure; any such diagnostic belongs
  in supplementary material under FD17.
- During isolated panel generation, no new CSV/TSV/RDS output, TeX change, or
  shared-ledger/review-report edit was made. Package-level provenance ledgers
  and the review report are assembled separately.

## Blockers

None. The requested saved tables were present and the final render completed.
