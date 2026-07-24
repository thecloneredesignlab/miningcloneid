# Gate B drafting panel record

## Package status

This package contains manuscript-ready review drafts for the five-figure
architecture approved at Gate A. It is a drafting package, not production
manuscript assembly: `ltee_hypoxia_model.tex` was used as read-only context and
was not edited. No optimizer was rerun, no necrosis export was repaired, and no
Figure 6 analysis was revived.

The recommended figures are under `final_figures/recommended/`. Each has a
300-dpi PNG and a one-page PDF generated from a local R script under
`scripts/`; PDFs are vector or hybrid where a preserved raster is part of the
approved content. Initial panel renders and refined composites are retained
for review transparency. Figures 1 and 3–5 and the optimizer supplement read
only the 98 immutable inputs under `source_tables/frozen_inputs/` (51.6 MB at
materialization); `source_tables/frozen_input_manifest.csv` maps every copy to
its upstream path and SHA-256 value. Figure 2 is generated from the tracked
model implementation. The manuscript, style reference, and Gate A records used
for scope are frozen under `source_context/`.

## Panel registry

| panel | purpose | generation and provenance | drafting rationale and caveat |
|---|---|---|---|
| 1A | Establish the matched in-vitro and in-vivo experimental design. | `scripts/make_figure1.R`; frozen observed-design tables under `source_tables/frozen_inputs/F1/`; copied prior generator and diff under `prior_code/F1A/` and `code_diffs/F1A.diff`. | Retains the historical design identity while simplifying to the longest tracked paths. Tumor oxygen is not portrayed as measured. |
| 1B | Show the observed in-vitro chromosome trajectories that motivate the result. | Same generator; frozen karyotype and lineage metadata. | Adds the observation requested in FD01. Cell-level distributions, medians, and interquartile ranges are descriptive. |
| 1C | Compare starting reference and terminal tumor chromosome distributions. | Same generator; frozen in-vivo ploidy and harvest tables. | Makes the in-vivo chromosome change visible. Starting and terminal cells are distinct samples, not longitudinal karyotyping of the same tumor. |
| 2A–E | Connect external input and latent resource state, fitted functions, state transitions, survival, and outputs in one schematic. | `scripts/make_figure2.R`; semantic sources `oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.R` and `.cpp`. | The A–E mapping follows the approved architecture. Measured burden is a fit target; simulated live-state demand feeds tumor resource state. Missegregation uses a multi-copy `±ΔN` shift, and the output histogram is schematic. |
| 3A | Lead with the direct in-vitro growth fit. | `scripts/make_figure3.R`; saved seed-10 lineage summary. | Only growth over lineage passage is shown, with control and deprived branches retained. Controls were fitted and are not a held-out negative-control set. |
| 3B | Compare fitted chromosome-state fractions with observed karyotypes. | Same generator; saved seed-10 distribution and observed-karyotype tables. | Preserves both control and deprived conditions and the full saved chromosome-state support. |
| 3C | Show the fitted post-missegregation survival function. | Same generator; saved seed-10 ploidy function table. | A fitted function, not an independent measurement or proof of mechanistic necessity. |
| 3D | Show the model-implied nonviable-daughter fraction. | Same generator; saved seed-10 oxygen/ploidy function table. | Conditional model consequence only. Figure 3E and any restricted-model comparison are intentionally absent. |
| 4A | Lead with direct in-vivo burden and terminal chromosome fit adequacy. | `scripts/make_figure4.R`; saved seed-25 burden and terminal-distribution outputs. | In-sample fit evidence precedes interpretation. No held-out validation or inferential uncertainty is claimed. |
| 4B | Show target and latent model-implied effective-O2 trajectories. | Same generator; saved seed-25 O2 lag time courses. | Effective O2 is a fitted latent state, not a tumor oxygen measurement. |
| 4C | Show the fixed-O2 univariate separator triptych at 0%, 1%, and 5% O2. | Same generator; saved fixed-O2 attractor modes and parameter table. | AUC describes discrimination between fitted solution modes; it is not a causal-effect estimate. |
| 4D | Preserve the pooled in-vivo/in-vitro solution landscape. | Same generator; pinned git raster plus saved pooled coordinates. The byte-identical source copy and layout crop are retained in `initial_subpanels/F4D/`. | Both contexts, saved geometry, point universe, overlays, and all three in-vivo regions are unchanged. Only outer whitespace/legend-column layout is adapted. |
| 4E | Compare `p_mis_base` and `n_O` across `vi_C01`, `vi_C02`, and `vi_C03`. | Same generator; exact formal-region membership on saved coordinates and saved seed summaries. | Regions are fitted solution families, not tumor subtypes. The panel does not promote the separators as causal. |
| 5A | Lead the joint figure with direct fit adequacy in both contexts. | `scripts/make_figure5.R`; frozen outputs for the exact six `joint_pair_best` records in `source_tables/frozen_inputs/F5/selected_results.tsv`. | Pale points show six predictions and filled points their median. All comparisons are in-sample. |
| 5B | Show all 14 in-vivo/in-vitro parameter ratios for all six winners. | Same generator; all six `joint_soft_coupling.tsv` tables. | Preserves the full parameter view and common winner universe. Across-winner range is solution sensitivity, not confidence. |
| 5C | Compare fitted proliferation functions across contexts. | Same generator; saved in-vivo and in-vitro oxygen-response curves for all six winners. | Thin curves are individual solutions and thick curves are pointwise medians; differences are not causal identification. |
| 5D | Compare fitted stress-linked missegregation functions across contexts. | Same generator and six-winner oxygen-response exports. | Uses the identical winner universe as 5A–E. |
| 5E | Compare fitted post-missegregation survival across contexts. | Same generator and six-winner ploidy-function exports. | Necrosis remains absent because finite predicted necrosis values are unavailable in the saved joint exports. |
| S1A–C | Report objective landscapes, stored optimizer diagnostics, and solution multiplicity outside the main figures. | `scripts/make_fit_diagnostics.R`; frozen objective, selected-seed summary, and convergence tables for the separate fits and all six joint-pair runs. | Implements FD17 without crowding Figures 3–5. Optimizer flags are diagnostics, not confidence intervals or held-out validation. |

Detailed input counts, assertions, and scientific caveats are recorded in
`subagent_notes/figure1.md` through `subagent_notes/figure5.md` and
`subagent_notes/fit_diagnostics.md`. File-level SHA-256 provenance is in
`source_manifest.csv`; upstream-to-frozen mappings are in
`source_tables/frozen_input_manifest.csv`.

## Visual quality control

| artifact | final size | visual review |
|---|---|---|
| Figure 1 | 7.1 × 6.05 in | Panel tags, symbols, trajectories, distributions, axes, and legends inspected without clipping. |
| Figure 2 | 7.1 × 4.8 in | Callouts, arrows, formulas, type labels, simulated-live-state feedback, and output labels inspected without overlap; annotations are at least 7 pt. |
| Figure 3 | 7.1 × 8.8 in | Full-resolution inspection confirmed readable axes, facets, overlay points, curves, and panel tags. |
| Figure 4 | 7.1 × 9.0 in | Full-resolution inspection performed after manuscript-size reflow; all A–E labels and the preserved embedding remain within the canvas. |
| Figure 5 | 7.1 × 9.0 in | Regenerated at final size with a 2×2 fit block, full-width ratio panel, and compact C–E row; no clipping detected. |
| Supplementary fit diagnostics | 7.1 × 7.7 in | Selected points, objective axes, diagnostic cells, and uncertainty guardrail inspected without clipping. |

The executable package validator is
`scripts/validate_drafting_package.R`. It checks the PNG/PDF pairs, final-size
envelope, required ledgers, report coverage, source-manifest hashes, and
explicit Figure 3E/Figure 6 exclusions.

## Revision history

- 2026-07-24: created initial subpanels for Figures 1–5 from observed or saved
  fit outputs and assembled the first recommended composites.
- 2026-07-24: corrected Figure 1 start/terminal ordering and Figure 2 context
  colors and labels after full-resolution visual review.
- 2026-07-24: refined Figure 3 tag/title spacing while retaining its approved
  four-panel boundary.
- 2026-07-24: reflowed Figure 4 at manuscript width and height while preserving
  the pinned pooled embedding.
- 2026-07-24: reflowed Figure 5 natively from 16.4 × 13.5 in to 7.1 × 9.0 in;
  no scientific content or source changed.
- 2026-07-24: linked the supplementary optimizer table directly to saved
  seed-summary and convergence tables, replacing manually entered display
  constants with computed values.
- 2026-07-24: froze all active generator inputs plus the manuscript/style/Gate
  context so figure and manifest regeneration do not depend on ignored source
  bundles. `scripts/materialize_frozen_inputs.R` remains a provenance-refresh
  utility and is not part of the clean package build.

No conservative or exploratory final variant was produced because Gate A
resolved the figure architecture and left no competing visual decision that
required a second final draft.
