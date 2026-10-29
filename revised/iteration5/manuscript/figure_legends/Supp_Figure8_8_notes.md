# Supplementary Figure8_8: family-resolved in-vitro finite-time ploidy

The upper C01 and lower C02 blocks each reproduce the 2-by-5 time layout
of Supplementary Figure8_6. Times are 100, 200, 300, 400, 500 days in the
first row and 600, 700, 800, 900, 1000 days in the second row of each block.
Imposed oxygen ranges from 0 to 20% in steps of 0.1%. Blue, pink and
orange solid curves denote 2N, 4N and equal-cell-number 2N+4N starts.
Each curve is the median of the corresponding family's 500 optimizer
seeds; shading is the pointwise Q25-Q75 interval (type-7 quantiles), not
a confidence interval. The two blocks share axis limits and a legend.
No seeds are deduplicated and no interval is artificially widened.

This is a family-specific regrouping of the existing v12 seed snapshots,
not a new simulation. Every seed retains its complete 33-parameter
in-vitro fit, including its own native fitted p_misseg. The source uses
fit_result.rds$ctx$invitro_cfg, fitted Gaussian starts, crowding, the
original DT, and threshold-triggered random passage. Snapshots are taken
after integer-day passage from one continuous 1000-day trajectory; the
simulation is not restarted at each snapshot. Passage-repeat means are
formed within each seed before the across-seed summary. The source
chromosome-grid and imposed-oxygen interpretation limitations remain.

Source: figure8_time_snapshots_20261007_v12/supp_fig8_6_seed_snapshots.tsv.gz.
Output run: figure8_family_time_snapshots_20261008_v13. Existing pooled
Supplementary Figure8_6 and all other pre-existing figures are unchanged.
