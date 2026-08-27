# Supplementary Figure 6-3 paired-context weak-gap audit

Status: completed and validated on 2026-08-25. All writes are confined to `revised/iteration2`.

## Scientific question

Figure 6A marks conditions where at least half of retained joint-fit endpoints have dominant spectral gap below 0.005. Supplementary Figure 6-3 asks, separately for the paired *in vivo* and *in vitro* parameter vectors, whether those cells also lose the qualitative low/intermediate/high dominant-mean-ploidy conclusion or instead retain regime agreement while the asymptotic dominant mode separates slowly.

Only the nonnegative rank-1 population distribution is used for biological summaries. Signed or complex higher eigenmodes remain archived for mathematical audit and are not treated as alternative chromosome-number distributions.

## Analysis design

- Displayed pairs: C01Sc01, C02Sc01, and C03Sc02, abbreviated C01, C02, and C03.
- Eligibility: the pair-specific lowest joint full-MAP objective decile, 50 endpoints per pair.
- Context pairing: both parameter vectors come from each retained joint endpoint; the *in vitro* row does not use a separate-fit seed-10 vector.
- Exact endpoint deduplication: C01, C02, and C03 contain 50, 15, and 50 unique parameter endpoints in each context. Each is calculated once and its original optimizer-seed multiplicity restored.
- Grid: 201 oxygen concentrations and 60 fixed effective per-chromosome missegregation probabilities, or 12,060 cells per context--pair panel.
- Weak-gap cell: at least 50% of retained endpoints have `spectral_gap < 0.005`, matching Figure 6A.

Dominant mean ploidy is classified as low at `<=2`, intermediate at `>2 and <4`, and high at `>=4`. A stable class requires at least 90% of the 50 seed-weighted endpoints to agree. Local sensitivity compares a focal cell with its available immediately adjacent oxygen and fixed-missegregation cells.

## Main findings

### *In vivo*

- C01, C02, and C03 contain 3,287, 2,838, and 2,118 weak-gap cells.
- Every weak-gap cell meets the 90% three-class consensus rule.
- Stable-high fractions are 57.80%, 55.78%, and 13.60%; stable-intermediate fractions are 36.32%, 37.21%, and 68.98%; stable-low fractions are 5.87%, 7.01%, and 17.42%.
- A majority-endpoint local three-class switch occurs at 8.67%, 8.17%, and 4.86% of weak-gap cells.
- The endpoint-median maximum local ploidy change reaches at least one unit at 9.13%, 11.31%, and 4.39% of weak-gap cells.

### *In vitro*

- C01, C02, and C03 contain 1,527, 1,526, and 1,526 weak-gap cells.
- Stable-high fractions are 97.38%, 97.38%, and 97.44%; stable-low fractions are 2.55%, 2.56%, and 2.49%. No cell is stable intermediate, and one cell per pair fails the 90% consensus rule.
- A majority-endpoint local three-class switch occurs at 5.30%, 5.24%, and 5.24% of weak-gap cells.
- Median across-endpoint 90th-minus-10th percentile ploidy spreads are 0.000043, 0.000395, and 0.000083; the corresponding 90th percentiles are 0.000061, 0.00121, and 0.000109.
- The endpoint-median maximum local ploidy change reaches at least one unit at 5.30%, 5.24%, and 5.24% of weak-gap cells.

Supported interpretation: weak spectral separation generally does not remove the ensemble-supported low/intermediate/high dominant-ploidy classification in either context. The dominant weak-gap regime differs by context, and isolated boundary cells can remain sensitive in exact ploidy magnitude.

Unsupported interpretations: endpoint fractions are not posterior probabilities, biological replicate frequencies, cell-state fractions, or biological switching probabilities.

## Panels

- A: stable-low, stable-intermediate, stable-high, and mixed regime maps.
- B: fraction of retained endpoints with a local four-neighbor three-class switch.
- C: across-endpoint 90th-minus-10th percentile dominant-mean-ploidy spread.
- D: empirical cumulative distribution of the endpoint-median maximum local ploidy change.

Every panel uses C01--C03 as columns, *in vivo* as the first row, and *in vitro* as the second row. The dotted purple contour in panels A--C is the same majority-endpoint weak-gap boundary used to define the masked region.

## Validation and files

- Context validation: 7/7 checks passed.
- Data validation: 27/27 checks passed.
- Figure validation: 11/11 checks passed.
- Output manifest: 10/10 recorded files exist and match their current MD5 checksums.
- Figure: `Figures/supp_fig6-3_weak_gap_regime_robustness.png` and `.pdf`.
- Data: `data/Figures/Supp_Figure6_3/supp_figure6-3_weak_gap_regime_robustness.tsv` and `supp_figure6-3_weak_gap_pair_summary.tsv`.
- Entry points: `Code/Figures/data_Supp_Figure6_3.R` and `Code/Figures/draw_Supp_Figure6_3.R`.
- Implementation: `Code/Figures/util/analysis/si_figure6_eigenmodes.R` and `Code/Figures/util/analysis/figure6_context_extension.R`.
