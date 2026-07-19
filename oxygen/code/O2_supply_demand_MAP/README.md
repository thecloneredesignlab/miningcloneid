# O2 supply-demand MAP code

This directory owns the O2 supply-demand MAP fitting, fitted-parameter
simulation, analysis, visualization, reporting, local orchestration, and HPC
workflows.

## Protected core

The following directories are immutable:

- `model/`
- `optimizer/`

No source, cache, metadata, or README file may be added, removed, renamed, or
edited inside them. The tracked checksum baseline is
`docs/immutable_core_sha256.tsv`; the complete 18-file tree baseline, including
the pre-existing Rcpp cache and `.DS_Store` files, is
`docs/immutable_core_full_manifest.tsv`.

Verify both contracts with:

```bash
Rscript oxygen/tests/check_immutable_o2_core.R
Rscript oxygen/tests/check_immutable_o2_core.R --full
```

The protected dispatcher also requires these stable backend paths:

- `util/o2_supply_demand_map_shared.R`
- `util/o2_supply_demand_map_common_semantics.R`
- `util/o2_supply_demand_map_fit_invivo_backend.R`
- `util/o2_supply_demand_map_fit_invitro_backend.R`
- `util/o2_supply_demand_map_fit_joint_backend.R`

Their paths and fitting-facing interfaces remain available.

## Layer contract

| Folder | Sole responsibility |
|---|---|
| `model/` | Protected numerical model implementation. |
| `optimizer/` | Protected fitting and profile-likelihood dispatchers. |
| `util/` | Shared, reusable functions and stable fitting backends. Libraries must not draw figures or launch a workflow when sourced. |
| `simulation/` | Use completed fitted parameters/results to materialize concrete O2, ploidy, CIN, population, functional-response, perturbation, and related numerical products. |
| `analysis/` | Consume materialized data and create statistical, diagnostic, clustering, ranking, or plot-ready tables. Analysis must not invoke simulation or draw figures. |
| `vis/` | Consume simulation/analysis tables and create figures plus visualization manifests. It must not read fit RDS objects, best-parameter files, or model code. |
| `report/` | Assemble existing tables and figures into reports. It must not simulate, analyze, or create analytical plots. |
| `runner/` | Sequence independent layers and preserve local compatibility workflows. |
| `hpc/` | Slurm submitters, array workers, dependency jobs, and HPC-only execution wrappers. |
| `docs/` | Architecture, migration, regression, and per-file reference material. |

The canonical dependency direction is:

```text
fit outputs / best parameters
              |
              v
          simulation
              |
              v
           analysis
              |
              v
              vis
              |
              v
            report
```

Each layer may use `util/`. Only `runner/` and `hpc/` may orchestrate multiple
executable layers.

## Canonical post-fit workflow

For a completed in-vivo, in-vitro, or joint seed, use:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/run_postfit_pipeline.R \
  --fit_dir=/absolute/path/to/seed \
  --scope=invivo
```

Valid scopes are `invivo`, `invitro`, and `joint`. The runner executes the
enabled stages in simulation -> analysis -> visualization -> report order and
writes one log per invoked stage. Use `--dry_run=TRUE` to inspect the exact
commands.

Specialized workflows have their own runner subfolders, including fixed-O2,
fit-results, parameter-landscape, profile-likelihood, warm-start, and
multi-warmup orchestration.

## In-vitro migration

The former top-level `oxygen/code/in-vitro-utils/` directory has been removed.
Its reusable `ivt_*` functions now live in canonical `util/` modules, loaded by:

```text
util/o2_supply_demand_map_invitro_utils.R
```

Pure in-vitro plot constructors live under `vis/invitro/`. All 58 historical
`ivt_*` functions and their formal arguments remain available, and the seed10
objective/table/figure goldens are preserved.

## Compatibility paths

When an established path had external callers, the old file is retained as a
small, explicitly labelled deprecated compatibility wrapper or orchestrator.
New numerical, analytical, visualization, and report logic lives only at its
canonical layer path. Compatibility mappings are recorded in
`docs/path_migration_table.tsv`.

## Regression gates

Run the full unit suite:

```bash
Rscript oxygen/tests/run_unit_tests.R
```

Or run the complete parse, shell/Python syntax, architecture, unit,
immutable-core, and whitespace gate in one command:

```bash
Rscript oxygen/tests/run_o2_reorganization_regression.R
```

Run the global architecture and immutable-core gates:

```bash
Rscript oxygen/tests/check_o2_layer_boundaries.R
Rscript oxygen/tests/check_immutable_o2_core.R --full
```

Before handoff, all R files are parsed, all shell/Slurm scripts are checked with
`bash -n`, Python entrypoints are compiled with an external bytecode cache, and
`git diff --check` is run. Stage-specific real-data and frozen-fixture evidence
is recorded in:

- `docs/stage0_regression_baseline.md`
- `docs/stage1_regression.md`
- `docs/stage2_regression.md`
- `docs/stage3_regression.md`
- `docs/stage4_regression.md`
- `docs/stage5_regression.md`

## Documentation

- `docs/README.md`: operational workflow guide and command examples.
- `docs/path_migration_table.tsv`: old-to-canonical path map.
- `docs/CODE_FILE_REGISTRY.md`: concrete responsibility, input/output, function,
  and direct-test information for every source file.
- `docs/code_file_registry.tsv`: machine-readable form of the registry.
- Layer and feature subfolders contain focused READMEs for their own files.

Regenerate the per-file registry after structural changes:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/documentation/generate_code_file_registry.R \
  --workflow_root=oxygen/code/O2_supply_demand_MAP
```
