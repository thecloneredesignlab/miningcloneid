# Migration: legacy 06 combine workflow

The result root and filenames remain unchanged under `oxygen/results/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/pooled_embedding_curve_class`.

| Legacy entrypoint | Canonical responsibility |
|---|---|
| `analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/calculate_regression_curve_average_slope.R` | `analysis/combined_fixo2_eigen/calculate_regression_curve_average_slope.R` |
| numerical table work formerly inside `plot_fixo2_eigen_attractor_embedding_curve_class.R` | `analysis/combined_fixo2_eigen/prepare_fixo2_eigen_curve_class_tables.R` |
| `analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/plot_fixo2_eigen_attractor_embedding_curve_class.R` | `vis/combined_fixo2_eigen/plot_fixo2_eigen_attractor_embedding_curve_class.R` |
| `analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/render_fixo2_eigen_attractor_embedding_curve_class_report.R` | `report/combined_fixo2_eigen/render_fixo2_eigen_attractor_embedding_curve_class_report.R` |
| sequential calls in the best-fit-parameter-feature runner | `runner/combined_fixo2_eigen/run_fixo2_eigen_curve_class_pipeline.R` |

The three historical files remain as deprecated compatibility wrappers. The visualization now requires the analysis manifest and never reads the dense-grid class or slope inputs directly.
