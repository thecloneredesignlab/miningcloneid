# Supplementary Figure8_6: in-vitro finite-time oxygen–ploidy matrix

This is the time- and oxygen-range extension of Figure8D, not a fixed-p version of
supp_fig8_5. C01 500 and C02 500 raw optimizer seed endpoints are pooled at
each oxygen coordinate before computing median and Q25-Q75. The 201 oxygen
coordinates are 0-20% in steps of 0.1%. Solid blue, pink and orange curves
represent the pure fitted 2N cohort, pure fitted 4N cohort and their 1:1
cell-number mixture respectively. No Q annotations appear in this export.

Each seed retains its complete 33-parameter in vitro final fitted vector,
including its own native p_misseg and Gaussian initialization parameters.
Only the imposed oxygen coordinate changes. Effective missegregation remains
model-derived; neither the native nor effective probability is held at a
common ensemble value. The initial mixture combines normalized fitted 2N/4N
distributions at the same total inoculum as the pure starts and is simulated
directly, never reconstructed from the arithmetic mean of pure terminal
ploidies.

The ten finite-time snapshots are 100, 200, 300, 400, 500 actual days in
the first row, and 600, 700, 800, 900, 1000 days in the second row. All panels
use common axes. Each condition is simulated as one uninterrupted 1000-day
trajectory: the model is not restarted at the panel times. The snapshots
are recorded after any integer-day threshold passage, following the same
convention as the original 1000-day endpoint. The canonical culture protocol
reseeds 569799 live cells on the first integer day reaching 6180783, with
randomized integerization and multivariate hypergeometric sampling without
replacement. There are 200 nominal independently keyed passage repeats per
seed/condition. Ploidy at each snapshot is averaged within seed across those repeats,
then median/IQR are calculated across 1000 optimizer seeds. A proven upper
bound permits one deterministic evaluation if the passage threshold cannot
be reached and all nominal repeats would be identical. No threshold reached
means no artificial passage.

All 0-20% trajectories reuse the unchanged, approved culture-repair engine. Every
seed uses `fit_result.rds$ctx$invitro_cfg`, including the saved DT=0.05-day
Euler update and population-dependent crowding. Division-linked terms are
crowding-scaled; direct death is not. The full natural-scale parameter
vector is authenticated before prediction. No preceding v4/v9/v10 culture
endpoint is reused, and no original source model or fit result is modified.
The matrix contains 6030000 seed/time/condition records and 6030 median/IQR
points. Its day-1000 panel reproduces all 603000 existing repaired endpoints
and their 603 median/IQR coordinates. Figure8D is that panel's exact 0-5%
subset, and the pure-2N endpoints at 0.5%
and 20% oxygen are the culture model anchors in Figure8A.

Shading is the optimizer-seed Q25-Q75 IQR, not a confidence interval,
posterior credible interval or biological replicate variability. Thin or
zero-width ribbons are not artificially enlarged. Chromosome-grid and
boundary-sensitivity limitations of the source model still apply. Model
oxygen percent is not a calibration to reported organ oxygen mmHg.

The authenticated chromosome grid is 22-154 autosomal chromosomes, with
22 chromosomes per N unit. Lower-three-state (22-24) and upper-three-state
(152-154) day-1000 masses were saved with every seed/condition in the source repair and summarized
in `audit/figure78_invitro_repair_20261007/figure8_full_ensemble_diagnostics.tsv`.
The minimum terminal log-population over passage repeats and population-floor
flags are also recorded; they are not trajectory minima. At 0.5% oxygen the
repaired 2N-start and 1:1 day-1000 medians are 1.051131N and 1.068440N;
median lower-three-state masses are 99.8196% and 98.4691%. Pure 4N instead
has a median of 3.589506N and negligible lower-boundary mass. At 20% oxygen
the three medians are approximately 5.9936N, with upper-three-state mass
approximately 1.76%.
No expanded-grid sensitivity run was requested or performed, so neither
plateau is demonstrated to be independent of chromosome-grid truncation.

Run: `figure8_time_snapshots_20261007_v12`. Numerical source and unchanged
engine: `figure8_invitro_repair_20261007_v11`. Snapshot-storage pilots verify
that recording daily states leaves both day-1000 ploidy and passage counts
identical. The former single-panel export is backed up under
`audit/figure8_time_snapshots_20261007/before/`. Figure8 main panels and
supp_fig8_1/2/3/4/5 remain unchanged. The companion in-vivo matrix is
supp_fig8_7; no supp_fig8_8 was created in that 2026-10-07 extension.

## Family-specific companion, 2026-10-08

The separately requested supp_fig8_8 splits the same seed snapshots into
an upper C01 and lower C02 block, each a complete 2-by-5 time matrix with
500 seeds, within-family median and Q25-Q75 IQR. No new simulations are
performed and this pooled supp_fig8_6 is unchanged.
