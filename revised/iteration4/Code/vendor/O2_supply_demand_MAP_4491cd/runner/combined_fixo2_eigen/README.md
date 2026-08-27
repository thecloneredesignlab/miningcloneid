# Combined FixO2 eigen-attractor runner

- `run_fixo2_eigen_curve_class_pipeline.R` executes the four canonical stages in order: slope analysis, annotated-table preparation, visualization, and HTML report. `--run_report=false` stops after visualization for compatibility with the historical best-fit-parameter-feature runner.

All legacy `--key=value` data-path options are forwarded to the owning stage; the established result root and filenames are preserved.
