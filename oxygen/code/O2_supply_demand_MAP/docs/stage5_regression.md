# Stage 5 regression: shared utilities and single-purpose layers

Stage 5 consolidated functions used by more than one workflow under `util/` and
left executable numerical, analytical, visualization, report, runner, and HPC
logic in their owning layers. Compatibility files expose established paths but
do not duplicate canonical implementations.

## In-vitro shared API

- The former `oxygen/code/in-vitro-utils/` runtime directory is absent.
- The canonical loader `util/o2_supply_demand_map_invitro_utils.R` exposes all
  58 historical `ivt_*` functions with identical formal arguments.
- The Stage 5 migration baseline used seed10 objective
  `3.8525352626059366`; that historical numerical contract was intentionally
  superseded by the fixed-time, independent-lineage structure correction.
- Under the corrected structure, the seed10 parameter replay objective is
  `4.0074938125984376`, and schema/behavior regression tests replace the old
  byte-identical table requirement.
- Post-fit I/O and response helpers were moved to dedicated util modules and
  can be sourced from an unrelated working directory without launching a
  workflow.

## Other shared utility contracts

- Fixed-O2 shared numerical, formatting, mode assignment, table I/O, and
  validation helpers have canonical util paths; historical loader APIs remain
  available. The checked legacy/nested loaders exposed 45 functions each, and
  the analytical compatibility loader exposed 50.
- Reusable fitted-result discovery/loading and model-probe helpers live in
  `o2_supply_demand_map_postfit_input_utils.R` and
  `o2_supply_demand_map_postfit_probe_utils.R`; simulation modules retain only
  their workflow-specific numerical producers.
- Joint parameter table preparation lives in
  `o2_supply_demand_map_joint_parameter_utils.R`; the historical
  `o2_supply_demand_map_joint_parameter_plot.R` filename is a data-only loader,
  while plot constructors live under `vis/joint/`.
- Common path, TSV, and key/value-manifest operations are provided by
  `o2_supply_demand_map_shared.R` and are used by profile workflows.
- Best-fit-feature CLI, path, table-I/O, report, fixed-O2, curve-classification,
  parameter-landscape, multi-warmup, eigen-attractor, perturbation, and
  fit-result helpers each have separate util modules.
- Report-side optional table loading, value formatting, parameter annotation,
  dependency detection, preview rendering, and data-URI helpers are centralized
  in `util/o2_supply_demand_map_report_utils.R`; HTML escaping is centralized
  in `util/o2_supply_demand_map_html_utils.R`.
- Common runner/HPC boolean parsing, integer validation, command rendering,
  module loading, path checks, and run-provenance helpers are centralized in
  the source-only `util/o2_supply_demand_map_shell_utils.sh`. The historical
  `hpc/util/write_run_provenance.sh` path is a nine-line compatibility loader.
- Shared visualization scales live in
  `vis/o2_supply_demand_map_common_plot_utils.R`; this source-only plot helper
  neither reads fit objects nor writes figures.
- Shared process-fingerprint parsing, path/table/manifest handling,
  transformations, distances, feature scaling, and clustering helpers now live
  in `util/o2_supply_demand_map_process_fingerprint_utils.R`. The historical
  analysis utility is a thin loader, while
  `simulation/process_fingerprints/process_fingerprint_simulation_legacy_utils.R`
  retains only simulation-owned fit-artifact readers and numerical feature
  production. The canonical
  `simulation/process_fingerprints/process_fingerprint_simulation_utils.R`
  handles materialized simulation contracts.
- Shared ploidy-regime trajectory features, clustering, concordance, and
  diagnostics live in `o2_supply_demand_map_ploidy_regime_utils.R`; the
  historical analysis path is a thin loader and the simulation helper retains
  only simulation-owned collection/build functions.
- The real joint parameter analysis table retained SHA-256
  `71aaedc6eb948d8abb19d55fc4fad1225548626fc354961a2f49e4b8aa777329`,
  and its visualization retained SHA-256
  `4499746e19111a97afd7889f8e1b0f3afba1e624b5fa30fbeb94598deff00291`.

## Entrypoint and HPC placement

- Warm-start table builders live under `runner/warm_start/`; a real seed50 plus
  seed350 smoke produced a 28-row start table.
- Slurm submitters, array workers, and HPC-only R tasks live under `hpc/`.
- Warm-up joint result submission files live under
  `hpc/warm_up_joint_fitting_results_extra/`.
- Multi-warmup pair construction is no longer an executable under `util/`.
- Compatibility entrypoints are required to be explicit thin wrappers or
  orchestration-only files and are covered by path/CWD tests.
- The historical combined FixO2/eigen-attractor `06` workflow is split into
  `analysis/combined_fixo2_eigen/`, `vis/combined_fixo2_eigen/`,
  `report/combined_fixo2_eigen/`, and `runner/combined_fixo2_eigen/`; its
  visualization consumes analysis-annotated coordinate tables only.

## Regression gates

Component gates exercised during consolidation include:

- `test-invitro-utils-migration.R`
- `test-invitro-defaults.R`
- `test-stage5-util-consolidation.R`
- `test-fixed-o2-stage34-boundaries.R`
- `test-fit-results-stage-split.R`
- `test-process-fingerprint-stage-split.R`
- `test-profile-likelihood-stage-split.R`
- `test-combined-fixo2-eigen-stage-split.R`
- `test-code-file-registry.R`

The protected 18-file `model/` and `optimizer/` tree passed repeated complete
size, modification-time, and SHA-256 checks during component testing and both
the pre-unit and post-unit final gates. The frozen repository-wide pass was:

```bash
Rscript oxygen/tests/run_o2_reorganization_regression.R
```

Final aggregate result:

- 272/272 R sources parsed.
- 25/25 shell and Slurm sources passed `bash -n`.
- 3/3 Python entrypoints compiled with bytecode outside the source tree.
- 302 source files are documented in the generated per-file registry.
- The complete unit suite passed.
- Layer boundaries and both immutable-core checks passed.
- `git diff --check` passed.
- `O2_REORGANIZATION_REGRESSION=PASS`.
