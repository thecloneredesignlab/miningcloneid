# Figure 4 drafting notes

## Delivered scope

Figure 4 contains the five approved panels in reading order:

- A: compact seed25 in-sample tumor-burden and terminal chromosome-number fit
  adequacy; necrosis is excluded.
- B: target versus latent model-implied effective-O2 trajectories for all eight
  exposure schedules.
- C: one strongest univariate separator at each of O2 = 0%, 1%, and 5%, with
  explicit noncausal wording.
- D: the approved pooled in-vivo/in-vitro embedding, including its full context
  point universe, both objective overlays, saved geometry, `vi_C01`–`vi_C03`,
  and `vt_C01`–`vt_C02`.
- E: `log10(p_mis_base)` and `n_O` distributions across all three formal
  in-vivo regions.

The single active generator is `scripts/make_figure4.R`. It was run with
`scripts/agentRrunner.sh`. It reads saved tables and a pinned historical PNG;
it does not source or run the model, optimize parameters, simulate new data,
refit, or recompute an embedding.

## Directive coverage

- FD06 and FD17: Panel A leads with direct observed-versus-fitted evidence and
  leaves optimizer/objective diagnostics outside the main figure.
- FD08: the fixed-O2 triptych and pooled landscape remain in the main figure.
- FD09: the exact reproduced terms are used: `mu_hp` at 0% O2 and
  `p_mis_base` at 5% O2. AUC is labeled as univariate discrimination, not
  influence.
- FD15: `vi_C01`, `vi_C02`, and `vi_C03` are retained without remapping them to
  two regimes.
- FD16: the pooled in-vivo/in-vitro point universe, context markers, objective
  overlays, coordinates, and plot geometry are retained. There is no
  in-vivo-only embedding and no re-embedding.
- FD18: no refit, new analytical grid, necrosis reconstruction, TeX edit, or
  Figure 6 work was performed. Package-level provenance ledgers are maintained
  separately.

## Immutable inputs and checks

The active generator reads byte-identical copies under
`source_tables/frozen_inputs/F4/`. The paths described below are the upstream
sources; `source_tables/frozen_input_manifest.csv` records their frozen
counterparts, byte counts, and SHA-256 values.

Panel A reads:

- `separate_fits/invivo/selected_seeds/seed25/burden_fit.tsv` (120 finite
  observed/fitted burden rows across eight scenarios).
- `separate_fits/invivo/selected_seeds/seed25/viz/terminal_ploidy_observed_vs_predicted.tsv`
  (7,566 positive-weight rows across all eight endpoints).

Panel B reads
`separate_fits/invivo/selected_seeds/seed25/viz/o2_lag_timecourse.tsv`.
The blue latent trace is rendered first and the orange dashed target trace
second so coincident trajectories remain visibly encoded.

Panel C reads the saved analytical-attractor tables
`fixed_o2_attractor_mode_by_seed_o2.tsv` and `parameter_values_long.tsv`.
The generator reproduces rank AUC, selects the best feature per O2 level, and
asserts the approved low/high results:

- 0% O2: `mu_hp`, AUC 0.849146, larger in the lower-ploidy mode.
- 1% O2: `p_mis_base`, AUC 0.741 (displayed to three decimals), larger in the
  higher-ploidy mode.
- 5% O2: `p_mis_base`, AUC 0.903459, larger in the higher-ploidy mode.

The input materializer resolves the exact tracked blob
`7e72dca88caf784fc61271d87a1c0dfb564b8303:oxygen/figures/iteration1/fig4f_landscape_tsne_clusters.png`.
The preserved source is 1,628 × 1,430 pixels and 884,007 bytes, with MD5
`14cecff29ab4884823b84d83f0379119` and SHA-256
`e85bfd6116e28d82b5c9c0df276f8c0148f9282409306b5f7fe5315da90921f2`.
The active figure generator reads its frozen byte-identical copy. For the
7.1-inch manuscript canvas, only external whitespace is cropped; the original
legend column is cropped, vertically padded without rescaling, and placed
beside the unchanged square plot crop. No plot coordinate, point, hull, label,
context marker, or objective overlay is regenerated or moved relative to the plot.
The exact untouched source remains in
`initial_subpanels/F4D/figure4D_pooled_embedding_preserved_source.png`.

Panel E uses the 500 saved in-vivo best-fit t-SNE coordinates from
`pooled_invivo_invitro_initial_vs_best_tsne_best_points_curve_class.csv`.
The historical deterministic coordinate-clustering procedure (k = 2,...,8,
mean silhouette selection, seed 123, k-means `nstart = 10`) selects k = 3 and
reproduces the formal counts exactly: `vi_C01` = 99, `vi_C02` = 385, and
`vi_C03` = 16. All 500 rows match finite `p_mis_base` and `n_O` values in
`separate_fits/invivo/run/extra_results/seed_summary.tsv`. This regenerates the
formal labels from saved coordinates; it does not alter or recompute the
embedding.

## Prior-code fidelity

Prior generators were copied before adaptation into:

- `prior_code/F4B/viz_invivo_model_O2_supply_demand_MAP_results.R`
- `prior_code/F4C/mode_parameter_contribution_analysis.R`
- `prior_code/F4D/plot_pooled_embedding_curve_class.R`
- `prior_code/F4E/clustering_report.R`

Panel A had no identifiable prior direct-fit plot and is recorded as a new
saved-table view. Focused adaptation records are in `code_diffs/F4A.diff`
through `code_diffs/F4E.diff`.

## Scientific caveats

- Seed25 has total weighted-MAP objective 14.1193 and ranks first among 500
  starts, but 29 fits are within 1% and 241 within 5% of the best objective.
  Good in-sample agreement is therefore not parameter certainty or evidence
  for a unique mechanism.
- Necrosis observations exist, but all exported necrosis predictions are NA;
  the endpoint cannot support a direct fit-quality comparison.
- Latent effective O2 is inferred by the model and is not a direct tumor oxygen
  measurement.
- Fixed-O2 AUCs are one-feature descriptive separators between model-defined
  attractor modes. They are not causal effects, ablations, or necessity tests.
- The pooled t-SNE is a descriptive projection. Its formal regions represent
  fitted solution families, not biological subtypes, and embedding distance is
  not a calibrated biological distance.
- Panel E has no bootstrap, posterior, or parameter-uncertainty intervals.

## Dimensions and visual QC

- Final PNG: 2,130 × 2,700 pixels, 7.1 × 9.0 inches at 300 ppi, white sRGB
  background.
- Final PDF: one 511 × 648 point page (7.1 × 9.0 inches).
- Final output hashes are recorded in `source_manifest.csv`.
- The final PNG was inspected at full resolution and at intended physical size.
  Tags A–E, titles, subtitles, axes, facet strips, C direction labels, D
  overlays/legends, E region labels, and the three-line footer are inside the
  canvas. Programmatic text is at least 7 pt; the preserved D source remains
  legible at the final scale.
- The final reflow is A above B/C above D/E. Panel D is large at lower left;
  Panel E uses two stacked horizontal distribution facets at lower right.
  Panel A uses sparse 10/100/1k burden ticks, and Panel B draws the dashed
  target trace above the blue latent trace.

## Intentionally skipped outputs

- No necrosis panel or reconstructed necrosis prediction.
- No objective distribution or optimizer diagnostic in the main figure.
- No refit, ablation, alternate seed, causal feature analysis, uncertainty
  product, or new embedding.
- No edits to TeX or production manuscript files. Package-level provenance
  ledgers and the review report are assembled separately.

## Blockers

None.
