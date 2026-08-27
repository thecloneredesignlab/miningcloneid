# Supplementary Figure 6-3 weak-gap regime-robustness audit

Status: completed and validated on 2026-08-20. All calculations, code, tables, figures, and manuscript edits are confined to `revised/iteration2`.

## Scientific question

Figure 6A marks conditions where at least half of retained joint-fit endpoints have dominant spectral gap below 0.005. Supplementary Figure 6-3 tests whether those weak-gap cells also lose the qualitative low/intermediate/high dominant-mean-ploidy conclusion, or instead retain regime agreement while the asymptotic dominant mode separates slowly.

The analysis deliberately uses only the nonnegative rank-1 population distribution for biological summaries. Signed or complex higher eigenmodes remain available in the archived top-10 eigenspectrum audit, but they are not treated as alternative chromosome-number distributions.

## Analysis population and grid

- Displayed joint-fit pairs: C01Sc01, C02Sc01, and C03Sc02, abbreviated C01, C02, and C03.
- Eligibility: the pair-specific lowest full-MAP objective decile, 50 numerical endpoints per pair.
- Exact endpoint deduplication: 50, 15, and 50 unique 14-parameter endpoints for C01, C02, and C03. Each unique endpoint was calculated once, then its original q10 optimizer-seed multiplicity was restored.
- Grid: 201 oxygen concentrations from 0% to 5% and 60 log-spaced fixed effective per-chromosome missegregation probabilities from 0.005 to 0.500.
- Weak-gap cell: at least 50% of retained endpoints have `spectral_gap < 0.005`, matching the Figure 6A mask.

## Definitions

Each optimizer-seed endpoint is assigned one of three mutually exclusive classes from its rank-1 dominant mean ploidy:

- Low: `dominant_mean_ploidy <= 2`.
- Intermediate: `2 < dominant_mean_ploidy < 4`.
- High: `dominant_mean_ploidy >= 4`.

At each pair and grid cell, `q_low`, `q_intermediate`, and `q_high` are the corresponding fractions among the 50 retained optimizer-seed endpoints. Regime consensus is their maximum. A stable class requires its fraction to be at least 0.90; cells not meeting that requirement are mixed.

- Across-fit exact-value spread: 90th minus 10th percentile of dominant mean ploidy.

For each endpoint and focal grid cell, the available immediate neighbors are oxygen plus or minus one grid step and fixed missegregation plus or minus one grid step, for at most four neighbors.

- Local regime switch: at least one neighbor changes among the low, intermediate, and high classes.
- Local ploidy change: maximum absolute dominant-mean-ploidy difference to those neighbors.

Endpoint-level indicators and changes are calculated before restoring multiplicity and aggregation. These are one-grid-step numerical sensitivity measures, not derivatives, biological transitions, or switching probabilities.

## Main findings

- Weak-gap regions contain 3,287, 2,838, and 2,118 grid cells in C01, C02, and C03, representing 27.26%, 23.53%, and 17.56% of each complete surface.
- Every weak-gap cell has at least 90% endpoint agreement on the three-level dominant-mean-ploidy regime. Stable-high cells account for 57.80%, 55.78%, and 13.60%; stable-intermediate cells account for 36.32%, 37.21%, and 68.98%; stable-low cells account for 5.87%, 7.01%, and 17.42%; mixed cells account for 0% in C01, C02, and C03, respectively.
- A local three-class switch occurs for at least half of endpoints at 8.67%, 8.17%, and 4.86% of weak-gap cells. At least one endpoint switches at 8.76%, 8.35%, and 5.85%, respectively.
- Across-fit 90th-minus-10th percentile dominant-ploidy spread has medians of 0.000035, 0.000317, and 0.00672 and 90th percentiles of 0.000733, 0.0237, and 0.0107 ploidy units for C01, C02, and C03. Maxima are 0.0481, 0.653, and 0.142.
- The endpoint-median maximum local ploidy change has 90th percentiles of 0.891, 1.147, and 0.385 ploidy units. It reaches at least one unit at 9.13%, 11.31%, and 4.39% of weak-gap cells, showing that exact response magnitude and categorical boundary crossing are related but distinct diagnostics.

Supported interpretation: the Figure 6A weak-gap mask primarily identifies slow dominant-mode separation. It does not generally remove the ensemble-supported qualitative low/intermediate/high dominant-mean-ploidy conclusion. Narrow grid-boundary regions can remain locally sensitive in exact ploidy magnitude.

Unsupported interpretations: these endpoint fractions are not posterior probabilities, biological replicate frequencies, or the fraction of cells in a low, intermediate, or high-ploidy state. Establishing the last quantities would require integrating probability mass over the corresponding ranges within each rank-1 population distribution; the current statistic classifies the rank-1 distribution's mean ploidy.

## Figure panels

- A: stable-low, stable-intermediate, stable-high, and mixed regime maps restricted to weak-gap cells.
- B: fraction of retained optimizer-seed endpoints with a local four-neighbor three-class switch.
- C: across-fit 90th-minus-10th percentile dominant-mean-ploidy spread.
- D: empirical cumulative distribution of endpoint-median maximum local ploidy change.

The dotted purple contour in panels A-C is the same 50%-endpoint weak-gap boundary used to define the masked region.

## Validation

- 115/115 cached unique endpoint calculations passed the original operator QC.
- Weak-gap robustness table: 36,180 rows, comprising 3 pairs x 201 oxygen values x 60 fixed-missegregation values.
- All rank-1 medians and 10th and 90th percentiles match the Figure 6D source surface to numerical tolerance; the largest checked absolute difference is `6.22e-15`. The legacy `dominant_mean_ploidy >= 2` fraction is retained only as an independent source-surface crosscheck. The new low/intermediate/high endpoint fractions are mutually exclusive and sum to one at every grid cell.
- Data validation: 27/27 checks passed.
- Figure validation: 9/9 checks passed; PNG dimensions are 3960 x 4320; publication and manuscript copies have matching MD5 checksums.

## Delivered files

- Figure: `Figures/supp_fig6-3_weak_gap_regime_robustness.png` and `.pdf`
- Manuscript copies: `manuscript/Figures/supp_fig6-3_weak_gap_regime_robustness.png` and `.pdf`
- Data entry point: `Code/Figures/data_Supp_Figure6_3.R`
- Drawing entry point: `Code/Figures/draw_Supp_Figure6_3.R`
- Analysis and plotting implementation: `Code/Figures/util/analysis/si_figure6_eigenmodes.R`
- Grid-level robustness table: `data/Figures/Supp_Figure6_3/supp_figure6-3_weak_gap_regime_robustness.tsv`
- Pair-level summary: `data/Figures/Supp_Figure6_3/supp_figure6-3_weak_gap_pair_summary.tsv`
- Data validation: `data/Figures/Supp_Figure6_3/supp_figure6-3_data_validation.tsv`
- Figure validation: `data/Figures/Supp_Figure6_3/supp_figure6-3_figure_validation.tsv`
- Source provenance: `data/Figures/Supp_Figure6_3/supp_figure6-3_source_file_provenance.tsv`
- Output manifest: `data/Figures/Supp_Figure6_3/supp_figure6-3_output_manifest.tsv`

The original `top10_endpoint_cache/`, `archive_figure6_top10_eigenmode_summary.tsv`, and `archive_figure6_eigenmode_competition_summary.tsv` remain as a reproducibility and mathematical-spectrum audit. Rerunning `data_Supp_Figure6_3.R` without `--rebuild=true` validates and reuses those endpoint caches.
