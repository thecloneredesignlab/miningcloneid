# neighborhood10pct second-order Sobol pilot

Scope: eight representative fitted endpoints; natural-span ±10% clipped neighborhoods; five scrambles per endpoint. This is not a 500-endpoint result.

Oxygen: 0, 0.5, 2.5 and 5%; 14 parameters, all 91 pairs; N=256,512,1024.
Each neighborhood is averaged over five scrambles before taking the median across eight endpoints.
Input distributions retain the original independent uniform/log-uniform choices.
Indices describe relative variance within this scope; output variances are retained for comparison.
Raw negative estimates and undefined/unstable cells are retained. A negative estimate is sampling error, not antagonism.
SALib 95% bootstrap intervals are estimator diagnostics; randomized-QMC repeat/resolution checks are reported separately.
Diagnostic tolerances: five-repeat range <=0.05 and mean change N512 to N1024 <=0.02. These are operational checks, not a proof of convergence.
Conditional surfaces hold the other 12 parameters at the indicated optimum. Their nonadditive contrast is not the marginalized Sobol S2.
Pairwise S2 excludes third- and higher-order interactions. S2/ST estimates do not give biological effect direction.

S2 cells passing all applicable checks: 450/728.
For the neighborhood median, this count requires all eight endpoint cells to pass.

## Strongest estimated pair/output/oxygen cells

| Pair | Output | O2 (%) | S2 | Diagnostic |
|---|---|---:|---:|---|
| buffer_smax × buffer_beta | dominant_mean_ploidy | 2.5 | 0.04471 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_mean_ploidy | 5.0 | 0.04232 | mixed_or_not_passed |
| buffer_beta × O2_crit | dominant_mean_ploidy | 0.5 | 0.03061 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_growth_rate | 0.0 | 0.02864 | mixed_or_not_passed |
| buffer_beta × mu_hp | dominant_mean_ploidy | 0.0 | 0.01899 | mixed_or_not_passed |
| buffer_smax × n_O | dominant_mean_ploidy | 2.5 | 0.01149 | mixed_or_not_passed |
| buffer_beta × buffer_n_exp | dominant_mean_ploidy | 0.5 | 0.01105 | mixed_or_not_passed |
| buffer_smax × O2_crit | dominant_mean_ploidy | 2.5 | 0.01047 | mixed_or_not_passed |
| buffer_smax × buffer_beta | dominant_mean_ploidy | 0.0 | 0.01038 | mixed_or_not_passed |
| buffer_beta × n_O | dominant_mean_ploidy | 0.5 | 0.01012 | mixed_or_not_passed |
| buffer_beta × gamma_growth | dominant_mean_ploidy | 0.5 | 0.01008 | mixed_or_not_passed |
| buffer_beta × alpha_o2 | dominant_mean_ploidy | 0.5 | 0.01001 | mixed_or_not_passed |

These trial oxygen points do not estimate a continuous low/high-band AUC. Claims about the full oxygen range require the 201-point extension.
