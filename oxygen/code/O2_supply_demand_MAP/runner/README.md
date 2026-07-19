# Runner layer

Runner scripts sequence canonical workflow layers. They may launch fitting,
simulation, analysis, visualization, and report entrypoints, but they do not
implement those scientific operations themselves.

## Main entrypoints

| Entrypoint | Responsibility |
|---|---|
| `run_o2_fit.sh` | Unified local launcher for in-vivo, in-vitro, direct joint, joint-from-best-single-fits, and multi-warmup fitting modes. |
| `run_postfit_pipeline.R` | Standard completed-seed orchestration in simulation → analysis → visualization → report order for `invivo`, `invitro`, or `joint` scope. |
| `run_fit_model_O2_supply_demand_MAP.sh` | Stable low-level/compatibility fitting launcher used by existing automation. |
| `run_fit_joint_model_O2_supply_demand_MAP.sh` | Stable joint-only wrapper around the fitting launcher. |
| `run_multi_warmup_joint.sh` | Stable shell entry for multi-warmup joint fitting. |

## Specialized runner folders

| Folder | Orchestrated workflow |
|---|---|
| `fit_results/` | Extra-results, paired joint-sigma, sigma-burden, and long-ploidy selection stages. |
| `fixed_o2/` | Fixed-O2 simulation, analysis, visualization, and report stages plus legacy API loading. |
| `fixed_o2_eigen/` | Fixed-O2 eigen-attractor feature, clustering, visualization, and report stages. |
| `multi_warmup/` | Landscape, seed-plan, and collected-results workflows. |
| `parameter_landscape/` | Parameter-landscape simulation, analysis, visualization, and reporting. |
| `profile_likelihood/` | Materialized live-effective-`p_ms` comparison, visualization, and reporting sequence. |
| `warm_start/` | Start-table construction from completed separate fits. |
| `best_fit_parameter_feature/` | Compatibility orchestration for the staged best-fit feature workflow. |
| `documentation/` | Regeneration of the human- and machine-readable per-file code registry. |

## Standard post-fit call

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/run_postfit_pipeline.R \
  --fit_dir=/absolute/path/to/seed \
  --scope=joint
```

The runner accepts `--dry_run=TRUE` and per-layer enable flags so the exact
commands can be audited before execution. Layer-owned outputs remain under
their corresponding output directories; runner logs record each invoked step.

Shared subprocess and path helpers belong in `../util/`. Cluster resource and
Slurm dependency orchestration belongs in `../hpc/`.

Concrete purpose, inputs, outputs, functions, and direct tests for every runner
are listed in `../docs/CODE_FILE_REGISTRY.md` and
`../docs/code_file_registry.tsv`.
