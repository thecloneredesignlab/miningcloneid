# HPC layer

`hpc/` is the only location for Slurm submitters, array workers, dependency-job
wrappers, and task code intended specifically for cluster execution. Numerical,
analytical, visualization, and reporting logic remains in its functional layer;
HPC files only configure resources and invoke those canonical entrypoints.

## HPC folders

| Folder | Responsibility |
|---|---|
| `submit/` | Unified fit submission and multi-warmup task-table/dependency orchestration. |
| `array_workers/` | Slurm array workers for in-vivo, in-vitro, joint, fixed-O2, parameter-landscape, eigen-attractor, and multi-warmup tasks. |
| `postprocess/` | Dependent fit-result postprocessing and report rerender jobs. |
| `best_fit_parameter_feature/` | Top-level best-fit feature submission compatibility entry. |
| `dense_grid_monotonicity_classification/` | Dense-grid monotonicity submission. |
| `fix_o2_simulation/` | Fixed-O2 simulation array submission. |
| `fixo2_eigen_attractor/` | Fixed-O2 eigen-attractor task and array submission. |
| `parameter_landscape/` | Parameter-landscape full-workflow submission. |
| `warm_up_joint_fitting_results_extra/` | Joint warm-up curve-array and collection submissions. |
| `util/` | Shell-only HPC provenance helpers. |

## Primary submitter

Use the unified submitter for production in-vivo, in-vitro, and joint fitting:

```bash
bash oxygen/code/O2_supply_demand_MAP/hpc/submit/submit_o2_fit.sh \
  --fitting_mode=invivo \
  --config_path=/absolute/path/to/O2_supply_demand.yaml \
  --out_root=/absolute/path/to/results
```

Use `--dry_run=TRUE` to inspect submissions without calling `sbatch`. Detailed
fit modes and resource arguments remain documented in `../../../README.md`.

## Operational contract

- Launch submitters from a repository checkout visible to compute nodes, or
  pass the submitter's explicit project/config/result-root options.
- Slurm logs belong under the configured result/log root, not in this source
  tree.
- HPC wrappers must resolve canonical workflow paths from their own location or
  explicit arguments; personal workstation paths are not valid defaults.
- A local regression run checks scripts with `bash -n` and exercises dry-run or
  mock submission paths. It does not submit cluster jobs.

Concrete purpose, inputs, outputs, and direct tests for every HPC file are in
`../docs/CODE_FILE_REGISTRY.md` and `../docs/code_file_registry.tsv`.
