# Supplementary Figure 4-4: dynamics across all 500 in-vivo fitted endpoints

This standalone addition does not rerender Figure 4 or other supplementary figures.
Run the Slurm wrapper `Code/hpc/run_supp_figure4_4_hpc.sh` on hpctpa3pc0028 with
`SUPP4_4_CODE_COMMIT` set to the verified local task-source commit and
`SUPP4_4_CODE_CHECKSUMS` set to the transferred eight-file SHA256 manifest under
iteration5/audit. This supports the approved rsync handoff while GitHub push is
unavailable, without changing HPC Git refs or requiring an unpublished commit object.
It uses the established
R 4.4.2 SIF and read-only original fit caches; every generated file is in iteration5.

## Statistical unit and definitions

- Include seeds 1–500 exactly once, without objective/convergence filtering.
- Each seed contributes eight untreated tumor trajectories from its cached fit
  simulation. Cache parameters must agree with both best_params and the packaged
  Figure 4 500-seed table. Check manifest, configuration, grid, finite values,
  component totals, file dates and source hashes. Seed25 must numerically match
  the approved Figure4A time series. A failed check aborts rather than dropping a seed.
- Ploidy is the within-tumor viable-cell-weighted chromosome number divided by22.
  O2 is dynamic effective O2, not target O2 or the Figure4B fixed-O2 scan.
  Tumor volume includes viable and both dead-cell compartments, in mm³.
- First take the arithmetic mean across available tumors within seed/cohort/day.
  Then calculate R type-7 Q25/Q50/Q75 across the 500 cohort means. Do not pool
  seed/tumor observations or combine biological and optimizer variation.
- Retain the original follow-up: 2N day0–56, 4N day0–91. For 4N, four tumors are
  available through day77 and two from day78. Split both curves and ribbons and
  mark the membership change. Experimental observed means use their own original
  available-tumor denominator and are points only; they do not enter seed quantiles.
- The IQR describes optimizer-endpoint variation, not biological confidence
  intervals, Bayesian posterior uncertainty, or variation among stochastic repeats.

## Outputs

`Figures/Supp_Figure4_4.{pdf,png,svg}` is a 15×8-inch two-row, three-column figure.
Data and provenance are in `data/Figures/Supp_Figure4_4/`, including compressed
tumor-level and cohort-level seed traces, median/IQR summaries, seed validation,
source MD5 hashes, and plotting validation. The run audit additionally records
output SHA256 hashes and confirms existing Figure4-related files are unchanged.
