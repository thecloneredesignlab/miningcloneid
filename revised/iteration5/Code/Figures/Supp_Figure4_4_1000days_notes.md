# Supp Figure4_4: separate 1000-day projection

Keep the existing short-follow-up figure and its data unchanged. This variant
continues all eight original untreated tumor scenarios (four in each cohort)
for each of the 500 fitted endpoints. It reuses the existing day0–1000 forecast
caches without refitting. Each seed's predictions must reproduce all its original
568 tumor/time observations to relative tolerance1e-10, and the original fitted
inputs/configuration/cache metadata must have unchanged MD5 hashes.

Use the same two-stage statistics: equal-tumor mean within seed/cohort/day, then
R type-7 Q25/Q50/Q75 across all 500 seeds. Retain four model scenarios per cohort
through day1000; do not impose observed tumor withdrawal on future simulations.
The time-window diagnostic also reports the full range and 5th/95th percentiles
of mean viable ploidy, without changing the plotted IQR.

The original follow-up ends at day56 (2N) and day91 (4N). Later values are model
extrapolations, not validated long-term tumor predictions or biological confidence
intervals. Tumor size becomes extremely large in some fits: plot the raw-volume
median/IQR on a pseudo-log10 axis (sigma1mm3), retaining true zero observations.
Ploidy and effectiveO2 remain linear. Never interpret this display choice as a
change to the underlying raw-volume quantiles or as evidence of realistic growth.

Run `Code/hpc/run_supp_figure4_4_1000days_hpc.sh` through Slurm; hpctpa3pc0028 is
the default. The user approved another available compute node when that node
was fully booked. Override both `--nodelist` and `SUPP4_4_EXPECTED_NODE` together;
the wrapper verifies the actual allocated host and otherwise preserves the runtime.
Supply the local source commit and a SHA256 manifest of the eight task source files
through `SUPP4_4_CODE_COMMIT` and `SUPP4_4_CODE_CHECKSUMS`. The approved rsync handoff
is used while GitHub push is unavailable. All outputs remain inside iteration5.
