# Figure8: pooled context-specific finite-time predictions

The October 7 presentation amendment retains square A/B plotting frames.
A restores the original median-point and Q25-Q75 IQR error-bar style;
no violin densities are drawn. Blue in vivo and pink in vitro colors and
the original PCAWG/DepMap/Model point shapes are retained. The 39 displayed
summaries retain the preceding empirical sample selection and model seed counts;
the culture endpoints are replaced by the approved configuration repair below.
Missing empirical ploidy values remain excluded exactly as before.
Model summaries still pool C01 500 and C02 500 full-parameter seed endpoints
for each of the three displayed model conditions. No oxygen or ploidy
observations are synthesized from summary bars. A/B axis tick labels and
legends are enlarged to 12 pt, axis titles to 14 pt, and panel titles to
16 pt. B organ labels are enlarged from 2.5 to 3.5 mm, with minimal
borderless white backdrops to keep SD lines from crossing the text.
B retains the same
scatter observations, reported SD bars, OLS regression and confidence band.
The subsequent approved culture repair recomputes the in vitro values in A,
D and supp_fig8_1/2/5/6. B empirical inputs and analytical layers, all in vivo
values, and supp_fig8_3/4 exports remain unchanged. A retains its median/IQR
presentation, not its previous erroneous culture predictions.
Explanatory prose stays outside figure exports.

C now shows in vivo and D shows in vitro. Each main panel pools the raw
endpoints of C01 500 and C02 500 optimizer seeds. Both show only 2N, 4N and
the 1:1 cell-number mixture, using solid blue (#0072B2), pink (#CC79A7) and
orange (#E69F00) curves with Q25-Q75 ribbons. No Q annotation appears in
the main figure. Main C/D share oxygen and ploidy scales.

supp_fig8_3 retains the preceding fitted-p in vivo family-separated C/D
curves, all seven initial ratios and their original Q labels and numerical
layers. supp_fig8_4 combines the preceding fixed-p=0.01 and fixed-p=0.3
in vivo supplements into two matched panels, each pooling 1000 seeds.
supp_fig8_5 is their culture counterpart: a 2x2 matrix with columns C01/C02
and rows p_misseg=0.01/0.3, seven initial ratios per panel and 500 seeds
per curve point. Colors and ratios agree across the supplementary figures.

## B: reproduction of Beck et al., Cancer Letters Figure5

Source article: https://doi.org/10.1016/j.canlet.2025.217878
Figure caption: https://pubmed.ncbi.nlm.nih.gov/40517938/

The x axis is WGD-positive tumor fraction, not average ploidy. The y axis is
reported organ oxygen partial pressure (mmHg). The 19 original source rows
come from `share_bundle_figures_3_5/originals/organ_lut.csv`, checked against
Supplementary Table3. Vertical bars reproduce reported oxygen mean +/- SD;
one row has no SD and consequently no bar. The line and shaded region follow
the author's `oxygen_plots.Rmd`: ordinary least-squares oxygen versus WGD,
with the default 95% regression-mean confidence band. This band is not an
organ measurement CI and is unrelated to model ensemble IQRs.

Exact reproduction preserves repeated kidney/marrow rows and distinct
brain/lung oxygen summaries sharing broad tumor organ groups. These are not
19 independent organ cohorts. It also retains the author's apparent
suprarenal-to-thyroid mapping error. No primary-only filter, anatomical repair,
oxygen subsite exclusion or organ deduplication is silently applied. SD
endpoints below zero are preserved as reported mean +/- SD; they are not a
claim that measured oxygen can be negative. Local tables lack explicit unit
metadata and underlying measurement-level sample sizes; mmHg follows the
article's oxygen-pressure context rather than a new unit calibration.

Statistics differ across source versions. The saved 19-row Figure5/Table3
input gives Pearson R=0.6164156998, two-sided p=0.0049439394. The caption
reports R=0.62, p=0.0013. The earlier 20-row `nowicka.txt` gives R=0.6670465156,
p=0.0013154067. This is evidence of a version inconsistency, not proof of its
cause. No caption p-value is hard-coded into the reproduced analysis. Exact
comparisons are saved in `author_figure5_statistic_version_audit.tsv`.

## C: pooled in vivo finite-time initial-condition comparison

Model oxygen is 0-5% in steps of 0.1%, including 0.5%; the endpoint is day
1000. Each of the 1000 seeds retains all 26 native final in vivo parameters,
including its own fitted p_misseg. No endpoint deduplication, refitting,
parameter shuffling or common non-coupling parameter vector is used. Only
the imposed oxygen coordinate is specified externally.

Initial cohorts are pure 2N (44 autosomal chromosomes), pure 4N (88), or
a 1:1 cell-number mixture. Continuous fixed-oxygen state dynamics are linear.
The authenticated v8 finite-horizon predictions use a spectrally shifted
1000-day matrix exponential and exact state superposition. Mixture terminal
ploidies are weighted by the terminal cohort populations, not averaged with
the initial fractions. This is a finite-time prediction, not an eigenvector
equilibrium or an assumption of convergence.

The two families' raw seed endpoints are combined before calculating the
median and Q25-Q75. Averaging family medians or family quartiles is not used.

## D: pooled in vitro finite-time predictions with random passage

The same oxygen grid, day-1000 endpoint and three starting-cohort names are
used. Each seed uses all 33 native final in vitro parameters, including its
own fitted p_misseg and its own init_mean/init_sd parameters for 2N and 4N.
Thus culture starting cohorts use their fitted Gaussian chromosome-state
distributions, as in A, rather than replacing them with exact 44/88-state
point masses. The 1:1 start is the equal-cell-number combination of these
two normalized fitted distributions, with the same total inoculum as each
pure-cohort start.

The approved canonical culture protocol is retained: 569799 reseeded live
cells, threshold 6180783 live cells, and a passage on the first actual
integer day reaching that threshold. At passage, systematic randomized
integerization and multivariate hypergeometric sampling reseed cells without
replacement. Population size and composition are propagated across 1000
actual elapsed days; neither a fixed five-day dilution schedule nor continuous
growth without passage replaces this protocol. If the threshold is never
reached, no artificial passage is imposed. A proven no-crossing upper bound
may evaluate the deterministic path once because all 200 nominal repeats
would be identical; executed repeat counts are recorded separately.

Every culture mixture is simulated directly from its mixed initial vector.
The continuous in vivo pure-endpoint combination is not used for culture:
sampling and passage timing depend on the mixed population's state. There
are 200 independently keyed passage repeats per seed/condition. Each seed's
endpoint is the arithmetic mean of its repeats' population-weighted terminal
ploidies, followed by median and IQR across the 1000 optimizer seeds. The
200 repeats are not treated as additional fitted parameter seeds.

The corrected configuration is strictly loaded from each seed's
`fit_result.rds$ctx$invitro_cfg`, not from the in-vivo `fit_config.rds`.
Every seed's full 33 natural-scale parameters are checked against its final
`invitro_effective` parameter table. The saved DT=0.05-day Euler equations
apply the configured live-population crowding factor to division-linked
growth, maternal removal and WGD, but not to direct death. Crowding is
re-evaluated as the population changes, including after passage. This is the
fitting model's culture dynamics, not the earlier constant linear generator
with passage layered on top. No fit result, fitted parameter or original
model source is changed.

All pure and mixed culture starts are recomputed with newly generated
threshold-triggered random passage events. No previous v4/v9/v10 culture
endpoint or passage trajectory is reused. The inoculum and threshold are
the previously approved equal-lineage geometric summaries of observed
maintenance-passage records, not seed-specific fitted parameters. Historical
observed variable-oxygen lineages are replayed only as validation references;
production predictions impose one fixed oxygen for 1000 actual days.

The fitted-p 0.5% pure-2N endpoints are shared by D and A. The 20% culture
comparison stays in A and supp_fig8_6; main C/D and supp_fig8_5 use 0-5%.
D is the exact 0-5% numerical subset of the native-p 0-20% extension.

## Supplementary initial-ratio and fixed-p comparisons

Q is the initial number of 4N-cohort cells divided by the initial number of
2N-cohort cells. Pure 2N is Q=0, pure 4N is the Q=infinity limit, and 2N:4N
ratios 1:1 through 1:5 correspond to Q=1 through Q=5. Fractions are
1/(1+Q) and Q/(1+Q). Q is not the terminal lineage ratio and a mixture is
not a pure cohort at the initial average ploidy.

supp_fig8_3 preserves the v8 within-family fitted-p predictions and original
seven-Q plot layers. supp_fig8_4 preserves and combines the v7 fixed-p=0.01
and v8 fixed-p=0.3 pooled summaries. Their original numerical files remain
immutable. For continuous in vivo mixtures, the exact terminal 2N-origin
weight is `plogis(log_mass_2N - log_mass_4N - log(Q))`.

supp_fig8_5 uses complete culture parameter vectors and directly propagated
seven-Q culture initial distributions. Only the native scalar p_misseg is
overridden to 0.01 or 0.3. Effective missegregation is still derived by the
canonical model from native baseline, oxygen, buffering, death-rate and
coupling parameters; the effective probability itself is not fixed. All
other seed-specific parameters are retained. Its four panels summarize 500
seeds within each family after averaging passage repeats within each seed.

## Provenance and uncertainty

Current numerical run: `figure8_invitro_repair_20261007_v11`. It preserves
the 153 pooled main-C in vivo points and recomputes 1,319,000 culture
seed/condition endpoints in 6100 tasks: 603000 native-p endpoints for
three starts over 201 oxygen coordinates, 714000 fixed-p endpoints for
seven starts over 51 oxygen coordinates, and 2000 fixed-p pure-2N anchors
at 20% oxygen. Each coordinate has all 1000 optimizer seeds, without
endpoint deduplication. Outputs contain 153 pooled main-D points, 603
pooled day-1000 supp_fig8_6 points and 1428 within-family supp_fig8_5 points.

## Finite-time supplementary matrices, 2026-10-07

The current supp_fig8_6 replaces its former day-1000 single panel with a
chronological 2×5 in-vitro matrix, at 100-day intervals through 1000 days,
using 0–20% imposed oxygen. The new supp_fig8_7 is the corresponding in-vivo
matrix at 0–5% oxygen. Both pool 500 C01 and 500 C02 seeds and show the same
three starts (2N, 4N, and 1:1), native fitted p_misseg, median and Q25–Q75
IQR. Culture panels are snapshots of one uninterrupted threshold-passage
simulation; in-vivo panels are finite-time continuous-growth predictions.
The 1000-day subplots reproduce the preceding repaired culture summary
and existing Figure8C summary. Figure8 main panels, all other supplements,
and the source model are unchanged. No supp_fig8_8 was created in that
2026-10-07 extension. The former
single-panel export is retained in the snapshot audit's before/ directory.

Snapshot run: `figure8_time_snapshots_20261007_v12`; see the dedicated
Supp_Figure8_6 and Supp_Figure8_7 notes for numerical and uncertainty details.

## Family-resolved supplementary matrices, 2026-10-08

Supp_fig8_8 (in vitro, 0-20% oxygen) and supp_fig8_9 (in vivo, 0-5%)
regroup the existing v12 seed snapshots without new simulation. C01 is
the upper 2-by-5 block and C02 the lower 2-by-5 block; each block shows
days 100-1000 at 100-day intervals. Each family contributes its own 500
full-parameter seeds, with native fitted p_misseg, to its median and
pointwise Q25-Q75 IQR. The three starts and palette are unchanged, all
curves remain solid, and no Q labels or diagnostic prose are added to
the graph. Both blocks have common axis limits and a shared legend.
The original pooled supp_fig8_6/7 and all main panels are preserved.
Run: figure8_family_time_snapshots_20261008_v13. Reproduction scripts,
source hashes, numerical checks and rendered-figure review are recorded
under audit/figure8_family_time_snapshots_20261008.

The authenticated source-model files remain unchanged. Configuration,
complete parameter passthrough, timestep refinement, exact experimental
lineage/passage replay and stochastic lineage accounting are tested before
production; full endpoint counts, probability overrides, IQR ordering,
zero numerical clipping, D/extended-subset identity and protected empirical/
in-vivo inputs are checked after aggregation. Audit receipts and the previous
PDF/PNG copies are under `audit/figure78_invitro_repair_20261007/`.

All model ribbons are optimizer-seed Q25-Q75 IQRs, not confidence intervals,
biological replicate variability or posterior credible intervals. Thin or
zero-width IQRs are displayed at their actual width, without jitter or
artificial inflation. Previous chromosome-grid/boundary sensitivity limits
still apply. Model oxygen percent is not calibrated to empirical organ mmHg.

Terminal lower/upper chromosome-grid mass, minimum terminal log-population
across passage repeats, final population-floor flags and passage counts are
saved with every culture seed/condition. These final-population diagnostics
are not minima over the full trajectory. Correcting the culture configuration
does not establish agreement with observed hypoxic 2N ploidy increases;
at fitted native p_misseg, the hypoxic 2N-start and 1:1 culture predictions
still approach the model's lower chromosome boundary. At 0.5% oxygen their
pooled medians are 1.051131N and 1.068440N, with median mass in states 22-24
of 99.8196% and 98.4691%, respectively. The pure 4N-start median is instead
3.589506N and does not show that lower-boundary concentration. No expanded-grid
sensitivity run or empirical-direction adjustment is part of this repair.
