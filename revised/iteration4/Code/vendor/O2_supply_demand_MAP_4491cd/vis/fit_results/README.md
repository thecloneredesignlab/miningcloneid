# Fit-results visualization layer

| File | Responsibility | Materialized inputs | Outputs |
|---|---|---|---|
| `plot_extra_results.R` | Draw objective, parameter-boundary, ploidy, burden, and available joint-fit diagnostics. | Simulation and analysis manifests plus their TSV contracts. | PDF figures and `visualization_manifest.tsv`. |
| `plot_extra_results_objective_violin.R` | Draw only the objective-component distribution figure. | `objective_components_long.tsv` and analysis manifest. | `objective_components_violin.pdf`. |
| `plot_joint_sigma_soft_coupled_paired_seeds.R` | Draw paired objective, fitted-value, and vivo/vitro-ratio figures across sigma values. | Joint-sigma analysis tables. | PNG/PDF figures and visualization manifest. |
| `plot_sigma_burden_extra_results.R` | Draw count, objective, parameter, and rank-flow comparisons across sigma-burden caps. | Sigma-burden merged/count/flow analysis tables. | Nine PDF figures and visualization manifest. |

No file in this directory reads `best_params.tsv`, fit RDS objects, model code,
or raw simulation directories.
