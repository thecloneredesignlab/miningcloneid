# Utility and stable backend layer

`util/` contains reusable libraries and the stable fitting backends required by
the protected optimizer dispatcher. Utility files must not launch workflows,
draw figures, or write scientific results when sourced.

## Stable fitting interfaces

These paths are consumed by protected `model/` or `optimizer/` code and must
remain available:

- `o2_supply_demand_map_shared.R`: common CLI, scalar, path, table, manifest,
  and workflow-independent helpers.
- `o2_supply_demand_map_common_semantics.R`: canonical configuration and
  parameter semantics shared with the protected core.
- `o2_supply_demand_map_fit_invivo_backend.R`: in-vivo fitting backend.
- `o2_supply_demand_map_fit_invitro_backend.R`: in-vitro fitting backend.
- `o2_supply_demand_map_fit_joint_backend.R`: joint fitting backend.

## Shared utility families

| Family | Files and responsibility |
|---|---|
| In-vitro fitting | `o2_supply_demand_map_invitro_utils.R` is the canonical loader for the migrated `ivt_*` API. The `*_io_utils.R`, `*_lineage_utils.R`, `*_lineage_simulation_utils.R`, `*_summary_utils.R`, and `*_objective_utils.R` modules keep loading, lineage construction, candidate simulation, summaries, and fitting objective logic separate. |
| In-vitro post-fit | `o2_supply_demand_map_invitro_postfit_io_utils.R` and `o2_supply_demand_map_invitro_postfit_response_utils.R` provide shared post-fit table/response helpers without drawing figures. |
| Fixed O2 and eigen attractors | `o2_supply_demand_map_fixed_o2_utils.R`, `o2_supply_demand_map_fixed_o2_format_utils.R`, `o2_supply_demand_map_fixed_o2_mode_utils.R`, `o2_supply_demand_map_fixed_o2_table_utils.R`, `o2_supply_demand_map_fixed_o2_validation_utils.R`, and `o2_supply_demand_map_eigen_attractor_utils.R` provide shared numerical contracts, labels, modes, table/validation helpers, and feature helpers. |
| Fitted-result inputs and model probes | `o2_supply_demand_map_postfit_input_utils.R` centralizes reusable fitted-result discovery/loading, while `o2_supply_demand_map_postfit_probe_utils.R` provides reusable fitted-parameter model probes. |
| Fit-result workflows | `o2_supply_demand_map_fit_results_utils.R` provides CLI, manifest, TSV, validation, and subprocess helpers used across staged fit-result workflows. |
| Report support | `o2_supply_demand_map_report_utils.R` centralizes side-effect-free optional table reads, display-value normalization, parameter-table annotation, dependency checks, Ghostscript preview rendering, and data-URI encoding used by report assemblers. It does not draw figures or assemble reports. `o2_supply_demand_map_html_utils.R` provides shared HTML escaping. |
| Joint parameters | `o2_supply_demand_map_joint_parameter_utils.R` builds joint parameter tables. `o2_supply_demand_map_joint_parameter_plot.R` is a historical data-helper loader retained for compatibility; despite its filename it exposes no plot constructor. Plot code lives in `../vis/joint/`. |
| Parameter landscape and multi-warmup | `o2_supply_demand_map_parameter_landscape_io_utils.R` and `o2_supply_demand_map_multi_warmup_utils.R` provide shared table, manifest, and selection helpers. Executable workflows live in their functional layers. |
| Dense grid and combined embeddings | `o2_supply_demand_map_dense_grid_utils.R` provides dense-grid CLI, path, schema, and artifact-manifest contracts. `o2_supply_demand_map_combined_landscape_utils.R` provides pooled-coordinate discovery, path, and tabular I/O contracts. Neither module runs simulation, analysis, or plotting. |
| Combined FixO2/eigen annotations | `o2_supply_demand_map_combined_fixo2_eigen_utils.R` provides the shared table I/O and reduction/variant normalization used by the table-annotation analysis and consume-only visualization stages. |
| Shell runners and HPC launchers | `o2_supply_demand_map_shell_utils.sh` is the source-only canonical implementation of shared boolean, integer-validation, module-loading, command-rendering, warm-start labeling, path/file, and run-provenance helpers used by `runner/` and `hpc/` shell entrypoints. `hpc/util/write_run_provenance.sh` is only a compatibility loader. The canonical util does not launch or submit work when sourced. |
| Perturbations and curve classes | `o2_supply_demand_map_perturbation_utils.R` and `o2_supply_demand_map_curve_classification_utils.R` provide reusable contracts and classifications. |
| Process fingerprints | `o2_supply_demand_map_process_fingerprint_utils.R` is the canonical shared module for CLI/path handling, table/manifest operations, parameter transformations, distances, feature scaling, and clustering helpers used by both simulation and analysis. It does not read fit objects on source, run a producer, draw figures, or assemble reports. |
| Ploidy regimes | `o2_supply_demand_map_ploidy_regime_utils.R` provides the shared trajectory-feature, clustering, concordance, and diagnostic helpers used by process simulation and ploidy-regime analysis. |
| Best-fit feature compatibility | `o2_supply_demand_map_bpf_cli_utils.R`, `o2_supply_demand_map_bpf_path_utils.R`, `o2_supply_demand_map_bpf_table_io_utils.R`, and `o2_supply_demand_map_bpf_report_utils.R` preserve shared APIs for the staged best-fit feature workflows. |

The former top-level `oxygen/code/in-vitro-utils/` directory is not a runtime
or compatibility path. Its 58 historical `ivt_*` functions are loaded from
the canonical util modules above.

## Verification and file-level reference

The util boundary is checked by:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
Rscript oxygen/tests/run_unit_tests.R stage5-util-consolidation
```

Concrete purpose, inputs, outputs, exported functions, and direct tests for
every util file are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`.
