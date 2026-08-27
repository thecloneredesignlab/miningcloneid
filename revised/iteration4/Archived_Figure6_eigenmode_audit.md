# Archived Figure 6 top-10 eigenmode audit

Status: completed and validated on 2026-08-19; retained as a mathematical-spectrum provenance audit. On 2026-08-20, the manuscript-facing analysis was replaced by the rank-1 weak-gap regime-robustness analysis now documented as `Supplementary_Figure6-3_weak_gap_regime_robustness.md`, because signed or complex higher modes cannot be interpreted as alternative population distributions. All files remain confined to `revised/iteration2`.

## Scientific question

Figure 6A uses the dominant eigenmode and the classical gap between the first two real eigenvalues. This archived audit asks what a small gap means for ploidy: is the nearest competing mode localized near the same chromosome-number state, or can a nearly tied mode be localized in a different ploidy regime?

## Analysis population and grid

- Displayed joint-fit pairs: C01Sc01, C02Sc01, and C03Sc02, abbreviated C01, C02, and C03.
- Eligibility: the pair-specific lowest full-MAP objective decile, 50 numerical seeds per pair.
- Exact endpoint deduplication: 50, 15, and 50 unique 14-parameter endpoints for C01, C02, and C03. Each unique endpoint was calculated once, and its original q10 seed multiplicity was restored before pointwise summaries.
- Grid: 201 oxygen concentrations from 0% to 5% and 60 log-spaced fixed effective per-chromosome missegregation probabilities from 0.005 to 0.500.
- Calculation volume: 115 unique endpoints x 12,060 grid positions = 1,386,900 matrix eigendecompositions. The first 10 modes were retained at every condition, producing 13,869,000 endpoint-level rank-condition records in checkpoint matrices.

## Definitions and interpretation boundary

Eigenmodes are ordered by decreasing real eigenvalue. Rank 1 is oriented and nonnegative-normalized using the same rule as Figure 6D, so its mean ploidy remains a biological dominant-distribution summary.

Ranks 2-10 can be signed or complex and are not chromosome-number distributions. Their primary localization diagnostic is

`P_k = sum_N N * abs(v_Nk) / (N_UNIT * sum_N abs(v_Nk))`.

The gap to rank `k` is

`Delta_lambda_1k = Re(lambda_1) - Re(lambda_k)`.

An L2-amplitude localization was retained as a sensitivity diagnostic. The maximum pointwise median L1-versus-L2 localization difference was 2.323 ploidy units. Therefore, non-dominant localization is used only to diagnose the eigenspectrum geometry and is not interpreted as a second cell distribution.

Two derived quantities were calculated independently for each endpoint before seed-weighted aggregation:

- `J`: maximum `abs(P_k - P_1)` among ranks 2-10 with `Delta_lambda_1k < 0.005`; defined as zero if no rank meets the gap threshold.
- `G`: minimum `Delta_lambda_1k` among ranks 2-10 with `abs(P_k - P_1) >= 1`; undefined if no top-10 rank meets the localization threshold.

## Main findings

- A differently localized top-10 mode was present in at least half of retained endpoints at 59.3%, 61.9%, and 47.7% of the C01, C02, and C03 grids, respectively.
- These differently localized modes were usually not near-degenerate. Within grid points where at least half of retained endpoints contained such a competitor, the median competitor gaps were 0.0742, 0.0753, and 0.0897 for C01, C02, and C03.
- Narrow boundary-like regions did contain nearly tied modes with markedly different localization. Median `J` exceeded one ploidy unit at 1.64%, 1.52%, and 1.45% of the three grids.
- Maximum median `J` was 2.366 for C01 at O2 = 4.450% and fixed `p_miss,eff = 0.2895`; 1.907 for C02 at O2 = 4.825% and `p_miss,eff = 0.3130`; and 3.387 for C03 at O2 = 2.825% and `p_miss,eff = 0.5000`.
- The minimum pointwise median gap to a differently localized mode was 0.00116, 0.000133, and 0.000216 for C01, C02, and C03.

Supported interpretation: weak dominant-mode separation can coincide with a genuine alternative ploidy localization near response boundaries, but differently localized alternatives are usually separated by gaps well above 0.005. The audit does not establish biological switching probabilities or finite-time occupancy of non-dominant modes.

## Validation

- 115/115 unique endpoint checkpoints passed operator QC.
- Maximum relative eigen-equation residual: `8.6743e-15`.
- Top-10 summary: 361,800 rows; derived jump/gap summary: 36,180 rows.
- All 36,180 rank-1 pointwise medians matched the corresponding Figure 6D surface values; maximum absolute difference: `5.3291e-15`.
- The retired atlas passed 7/7 rendering checks at creation; PNG dimensions were 3960 x 5460 and its publication and manuscript copies matched. The live `si6_figure_validation.tsv` now validates the replacement weak-gap figure (9/9), while the atlas files remain unchanged as provenance.

## Delivered files

- Figure: `Figures/archive_fig6_eigenmode_competition.png` and `.pdf`
- Archived manuscript copies: `manuscript/Figures/archive_fig6_eigenmode_competition.png` and `.pdf`
- Historical data function: `si6_data()` in
  `Code/Figures/util/analysis/si_figure6_eigenmodes.R`. The current
  `data_Supp_Figure6_3.R` entry point does not run this archived audit; it
  reuses the fresh Figure 6 q10/q20 rank-1 caches needed by the displayed
  weak-gap robustness panels.
- Archived drawing function: `si6_draw_top10_audit()` in `Code/Figures/util/analysis/si_figure6_eigenmodes.R`
- Analysis and plotting implementation: `Code/Figures/util/analysis/si_figure6_eigenmodes.R`
- Endpoint manifest: `data/Figures/Supp_Figure6_3/archive_figure6_top10_endpoint_manifest.tsv`
- Endpoint QC: `data/Figures/Supp_Figure6_3/archive_figure6_top10_endpoint_qc.tsv`
- Top-10 summary: `data/Figures/Supp_Figure6_3/archive_figure6_top10_eigenmode_summary.tsv`
- Jump/gap summary: `data/Figures/Supp_Figure6_3/archive_figure6_eigenmode_competition_summary.tsv`
- Data validation: `data/Figures/Supp_Figure6_3/supp_figure6-3_data_validation.tsv`
- Current manuscript-facing figure validation: `data/Figures/Supp_Figure6_3/supp_figure6-3_figure_validation.tsv`; the historical atlas validation remains archived.
- Source provenance: `data/Figures/Supp_Figure6_3/supp_figure6-3_source_file_provenance.tsv`
- Archived output manifest: `data/Figures/Supp_Figure6_3/archive_figure6_output_manifest.tsv`

The historical checkpoint directory
`data/Figures/Supp_Figure6_3/top10_endpoint_cache/` is resumable when
`si6_data()` is called explicitly. It is not created by the current
manuscript-facing Figure 6 data-only workflow.
