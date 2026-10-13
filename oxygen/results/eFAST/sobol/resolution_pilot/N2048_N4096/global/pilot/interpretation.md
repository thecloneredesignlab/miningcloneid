# global second-order Sobol pilot

Scope: full independent in-vivo bounds; five randomized Sobol scrambles.

Oxygen: 0.0, 0.5, 2.5, 5.0%; 14 parameters, all 91 pairs; N=1024,2048,4096.
Each neighborhood is averaged over five scrambles before taking the median across eight endpoints.
Input distributions retain the original independent uniform/log-uniform choices.
Indices describe relative variance within this scope; output variances are retained for comparison.
Raw negative estimates and undefined/unstable cells are retained. A negative estimate is sampling error, not antagonism.
SALib 95% bootstrap intervals are estimator diagnostics; randomized-QMC repeat/resolution checks are reported separately.
Diagnostic tolerances: five-repeat range <=0.05 and mean change N2048 to N4096 <=0.02. These are operational checks, not a proof of convergence.
Conditional surfaces hold the other 12 parameters at the indicated optimum. Their nonadditive contrast is not the marginalized Sobol S2.
Pairwise S2 excludes third- and higher-order interactions. S2/ST estimates do not give biological effect direction.

S2 cells passing all applicable checks: 602/728.
For the neighborhood median, this count requires all eight endpoint cells to pass.

## Strongest estimated pair/output/oxygen cells

| Pair | Output | O2 (%) | S2 | Diagnostic |
|---|---|---:|---:|---|
| lam_max × p_misseg | dominant_growth_rate | 0.0 | 0.2371 | passed |
| lam_max × p_misseg | dominant_growth_rate | 0.5 | 0.1781 | not_passed |
| lam_max × p_mis_base | dominant_growth_rate | 5.0 | 0.1413 | passed |
| lam_max × p_mis_base | dominant_growth_rate | 2.5 | 0.1349 | not_passed |
| lam_max × p_misseg | dominant_mean_ploidy | 0.0 | 0.09402 | passed |
| lam_max × p_mis_base | dominant_growth_rate | 0.5 | 0.09065 | not_passed |
| lam_max × p_misseg | dominant_mean_ploidy | 0.5 | 0.08494 | not_passed |
| lam_max × O2_crit | dominant_growth_rate | 0.5 | 0.083 | not_passed |
| lam_max × mu_hp | dominant_growth_rate | 0.5 | 0.07577 | not_passed |
| lam_max × p_misseg | dominant_growth_rate | 2.5 | 0.06674 | not_passed |
| p_misseg × mu_hp | dominant_mean_ploidy | 0.0 | 0.06068 | passed |
| lam_max × O2_crit | dominant_growth_rate | 2.5 | 0.05734 | not_passed |

These trial oxygen points do not estimate a continuous low/high-band AUC. Claims about the full oxygen range require the 201-point extension.
