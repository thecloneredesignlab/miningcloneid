# neighborhood10pct second-order Sobol pilot

Scope: eight representative fitted endpoints; natural-span ±10% clipped neighborhoods; five scrambles per endpoint. This is not a 500-endpoint result.

Oxygen: 0.0, 0.5, 2.5, 5.0%; 14 parameters, all 91 pairs; N=1024,2048,4096.
Each neighborhood is averaged over five scrambles before taking the median across eight endpoints.
Input distributions retain the original independent uniform/log-uniform choices.
Indices describe relative variance within this scope; output variances are retained for comparison.
Raw negative estimates and undefined/unstable cells are retained. A negative estimate is sampling error, not antagonism.
SALib 95% bootstrap intervals are estimator diagnostics; randomized-QMC repeat/resolution checks are reported separately.
Diagnostic tolerances: five-repeat range <=0.05 and mean change N2048 to N4096 <=0.02. These are operational checks, not a proof of convergence.
Conditional surfaces hold the other 12 parameters at the indicated optimum. Their nonadditive contrast is not the marginalized Sobol S2.
Pairwise S2 excludes third- and higher-order interactions. S2/ST estimates do not give biological effect direction.

S2 cells passing all applicable checks: 640/728.
For the neighborhood median, this count requires all eight endpoint cells to pass.

## Strongest estimated pair/output/oxygen cells

| Pair | Output | O2 (%) | S2 | Diagnostic |
|---|---|---:|---:|---|
| buffer_smax × buffer_beta | dominant_mean_ploidy | 2.5 | 0.04852 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_growth_rate | 0.0 | 0.0299 | passed |
| buffer_smax × buffer_beta | dominant_mean_ploidy | 5.0 | 0.0223 | mixed_or_not_passed |
| buffer_beta × mu_hp | dominant_mean_ploidy | 0.0 | 0.022 | passed |
| buffer_smax × buffer_beta | dominant_growth_rate | 0.5 | 0.01333 | passed |
| buffer_beta × O2_crit | dominant_mean_ploidy | 0.5 | 0.01267 | mixed_or_not_passed |
| buffer_smax × O2_crit | dominant_mean_ploidy | 2.5 | 0.01115 | mixed_or_not_passed |
| buffer_smax × n_O | dominant_mean_ploidy | 2.5 | 0.008817 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_mean_ploidy | 0.0 | 0.00867 | passed |
| O2_crit × n_O | dominant_mean_ploidy | 2.5 | 0.008024 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_mean_ploidy | 0.5 | 0.007772 | mixed_or_not_passed |
| buffer_smax × alpha_o2 | dominant_growth_rate | 0.0 | 0.006384 | passed |

These trial oxygen points do not estimate a continuous low/high-band AUC. Claims about the full oxygen range require the 201-point extension.
