# Fit-results orchestration layer

| File | Ordered stages | Compatibility entry |
|---|---|---|
| `run_extra_results.R` | Prediction simulation → fit-table analysis → visualization → HTML report. | `analysis/fit_results/extra_results.R`. |
| `run_joint_sigma_soft_coupled_paired_seeds.R` | Paired-sigma analysis → visualization → report. | `analysis/fit_results/compare_joint_sigma_soft_coupled_paired_seeds.R`. |
| `run_sigma_burden_extra_results.R` | Sigma-burden analysis → visualization. | `analysis/fit_results/compare_sigma_burden_extra_results.R`. |
| `select_invivo_best_long_ploidy_seed.R` | Long-ploidy value materialization → selection analysis. | `analysis/fit_results/select_invivo_best_long_ploidy_seed.R`. |

Runners only sequence canonical entrypoints. Shared subprocess handling is in
`util/o2_supply_demand_map_fit_results_utils.R`.
