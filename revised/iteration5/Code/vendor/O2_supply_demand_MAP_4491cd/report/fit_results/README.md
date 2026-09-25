# Fit-results report layer

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `render_extra_results_report.R` | Assemble the existing rich HTML report from completed simulation, analysis, and visualization contracts. Separate layer directories are staged temporarily; no scientific values are recalculated. | Three upstream manifests, analysis tables, and figures. | `extra_results_report.html` and `report_manifest.tsv`. |
| `render_joint_sigma_soft_coupled_paired_seeds_report.R` | Assemble the paired-sigma HTML summary and figure gallery. | Joint-sigma analysis/visualization manifests and artifacts. | `joint_sigma_soft_coupled_paired_seed_comparison.html` and report manifest. |
| `collect_best_separate_fit_reports.py` | Rank completed separate fits, copy top HTML reports, and export raw/transformed parameter ranking CSVs. | In-vivo and in-vitro seed directories with fit summaries, best-parameter tables, and reports. | Two ranking CSVs plus copied top report HTML files. |

Report code does not run fitting, simulation, analysis, or plotting stages.
