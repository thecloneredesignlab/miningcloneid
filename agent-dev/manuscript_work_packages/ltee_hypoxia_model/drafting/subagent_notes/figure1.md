# Figure 1 drafting note

## Purpose and directive coverage

- FD01 addressed: the historical matched-design view is retained in simplified
  form and now includes the chromosome-number observations that motivate the
  study.
- FD18 respected: this is a replot of existing observations, not a fit, refit,
  manuscript edit, or new analysis grid.

## Generation and provenance

Local generator:
`agent-dev/manuscript_work_packages/ltee_hypoxia_model/drafting/scripts/make_figure1.R`.

The active generator reads byte-identical frozen copies under
`source_tables/frozen_inputs/F1/`. Their upstream paths are:

- `figures/invitro_lineage_timeline.tsv`
- `figures/invitro_passage_observations.tsv`
- `figures/invitro_kary_cells.tsv`
- `figures/invivo_burden_long.tsv`
- `figures/invivo_ploidy_cells.tsv`
- `figures/invivo_harvest_catalog.tsv`

`source_tables/frozen_input_manifest.csv` records the upstream-to-frozen
mapping, byte count, and SHA-256 for every file.

The prior design generator was copied to
`drafting/prior_code/F1A/plot_optimization_data_streams.R`; the mechanical
comparison is `drafting/code_diffs/F1A.diff`.

## Panel decisions

- A uses the longest tracked control and oxygen-deprived path for each starting
  ploidy to keep the timeline legible. Known target oxygen is encoded only for
  culture. The in-vivo strip shows burden-measurement times and terminal
  chromosome sampling; it does not portray measured tumor oxygen.
- B shows cell-level observed karyotypes, medians, and interquartile ranges.
  The shared starting reference is used to begin both summary trajectories.
  Control trajectories end earlier than oxygen-deprived trajectories.
- C compares starting cell-line reference karyotypes with terminal tumor-cell
  distributions. Each terminal diamond is one tumor median; violins summarize
  pooled cells. This is a start/end comparison, not longitudinal karyotyping of
  the same tumor.

## Visual QC

Inspected `final_figures/recommended/figure1.png` at its 7.1-inch draft width.
Panel labels, legends, axis labels, and measurement symbols are visible and not
clipped. 2N/4N identity uses the repository blue/vermilion encoding; lineage
condition also uses line type and shape. Low oxygen remains the darkest blue.
The large cell-level terminal distributions are summarized rather than obscured
by overplotting thousands of points.

## Caveats

- The design strip deliberately omits secondary branch variants; it is a
  schematic of the longest tracked paths, not a complete lineage tree.
- In-vitro chromosome measurements are sparse time-point distributions.
- Starting cell-line and terminal tumor measurements are different samples.
- No fitted mechanism is assigned in this figure.
