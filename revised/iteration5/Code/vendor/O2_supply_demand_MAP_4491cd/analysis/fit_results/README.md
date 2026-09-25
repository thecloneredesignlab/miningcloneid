# Fit-results analysis layer

This directory now contains table/statistical analysis only. Canonical scripts
do not draw figures, render reports, call model code, or invoke simulation.
Historical paths remain as explicitly marked compatibility wrappers so existing
HPC and runner calls keep working.

## Canonical files

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `analyze_extra_results.R` | Build per-seed objective, convergence, parameter-boundary, ranking, and optional joint/in-vitro summary tables. | Completed fit seed directories plus `simulation_manifest.tsv` and `prediction_gate_metrics.tsv`. | `seed_summary.tsv`, `convergence_summary.tsv`, `parameter_boundary_long.tsv`, objective tables, optional joint/in-vitro tables, `analysis_manifest.tsv`. |
| `extra_results_analysis_utils.R` | Pure fit-summary, parameter-transform, boundary-distance, convergence, ranking, and table-shaping helpers. | In-memory fit/parameter tables. | In-memory analysis data frames only. |
| `analyze_joint_sigma_soft_coupled_paired_seeds.R` | Compare paired seeds across soft-coupling sigma runs. | Existing `fit_summary.tsv` and `joint_soft_coupling.tsv` files. | Paired summary, two-sigma contrast, value-long, objective-long, and analysis manifest TSVs. |
| `analyze_sigma_burden_extra_results.R` | Compare already-materialized `extra_results/seed_summary.tsv` tables across sigma-burden caps. | One analysis table per sigma cap. | Merged seed table, group counts, Sankey-flow table, and analysis manifest. |
| `select_invivo_best_long_ploidy_seed_from_metrics.R` | Apply the long-ploidy eligibility threshold and objective tie-breaking to a simulation table. | `invivo_long_ploidy_metrics.tsv`. | Selection TSV, selected directory pointer, and analysis manifest. |
| `select_best_seed_from_summary.R` | Select the best eligible seed from an existing summary table and required-file contract. | `seed_summary.tsv` and seed-directory existence checks. | Selection audit TSV and selected directory pointer. |

## Compatibility files

| Historical file | Canonical destination / behavior |
|---|---|
| `extra_results.R` | Deprecated orchestrator for `runner/fit_results/run_extra_results.R`; runs simulation, analysis, visualization, then report. |
| `compare_joint_sigma_soft_coupled_paired_seeds.R` | Deprecated orchestrator for the joint-sigma runner. |
| `compare_sigma_burden_extra_results.R` | Deprecated orchestrator for the sigma-burden runner. |
| `select_invivo_best_long_ploidy_seed.R` | Deprecated orchestrator for long-ploidy simulation plus selection analysis. |
| `extra_results_report.R` | Deprecated wrapper for `report/fit_results/render_extra_results_report.R`. |
| `plot_extra_results_objective_violin.R` | Deprecated wrapper for `vis/fit_results/plot_extra_results_objective_violin.R`. |
| `collect_best_separate_fit_reports.py` | Deprecated wrapper for `report/fit_results/collect_best_separate_fit_reports.py`. |

Shared CLI, TSV, manifest, and subprocess helpers live in
`util/o2_supply_demand_map_fit_results_utils.R`.

## Regression contract

The split was checked against the repository's ten-seed in-vivo fixture. Nine
legacy TSVs are byte-identical, including the `10 x 329` seed summary, the
`210 x 24` parameter-boundary table, and the four 1000-day prediction tables.
The joint-sigma four-table contract, sigma-burden three-table contract,
long-ploidy selection table/pointer, and Python ranking CSVs are also
byte-identical to their pre-migration producers.
