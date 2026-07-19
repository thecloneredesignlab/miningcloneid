# Joint in-vivo/in-vitro visualization

## `plot_invivo_invitro_comparison.R`

- Purpose: create matched in-vivo versus in-vitro functional-response figures.
- Inputs: `--fit_dir`, or explicit `--invivo_simulation_dir` and
  `--invitro_simulation_dir`; each simulation directory must contain
  `functional_curve_oxygen.tsv`,
  `functional_curve_oxygen_multi_ploidy.tsv`,
  `functional_curve_ploidy.tsv`, and
  `functional_curve_ploidy_by_o2.tsv`.
- Outputs: nine PDF/PNG figure pairs and `visualization_manifest.tsv` under
  `--out_dir` (default `<fit_dir>/viz/invivo_vs_invitro`).
- Failure contract: missing simulation directories or tables are errors. The
  script never runs a model or a simulation producer.
- Called by: `runner/run_postfit_pipeline.R --scope=joint`.
- Consumed by: `report/render_fit_report.R`, which only embeds pre-existing
  figures.
- Regression: real in-vivo seed366 and in-vitro seed10 tables produced the same
  nine PNG files as the former report-embedded implementation; PDF dimensions
  and byte sizes match, with expected PDF metadata differences.

## `plot_joint_parameter_ratios.R`

- Purpose: draw the joint in-vivo/in-vitro fitted-parameter ratio figure.
- Inputs: `joint_parameter_ratio.tsv` and
  `joint_parameter_ratio_metadata.tsv` from
  `analysis/fit_diagnostics/run_joint_parameter_diagnostics.R`.
- Outputs: one PDF/PNG pair plus
  `joint_parameter_visualization_manifest.tsv` under
  `<fit_dir>/viz/joint_parameters` by default.
- Failure contract: missing analysis tables are errors; this entrypoint never
  reads `best_params.tsv`, `fit_result.rds`, or model code.
- Called by: `runner/run_postfit_pipeline.R --scope=joint`.
- Regression: seed10 PNG is byte-identical to the former direct plotting
  implementation; PDF size matches with expected metadata differences.

## `o2_supply_demand_map_joint_parameter_plot_utils.R`

- Purpose: pure ggplot construction and label/color helpers for the joint
  parameter-ratio table.
- Inputs: an already materialized analysis data frame.
- Outputs: a ggplot object only.
- Called by: `plot_joint_parameter_ratios.R`.
