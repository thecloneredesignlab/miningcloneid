# Warm-up joint fitting HPC entrypoints

This folder contains the HPC-only launchers for the compatibility entrypoint
`analysis/warm_up_joint_fitting_results_extra/warm_up_joint_fitting_results_extra.R`.
Canonical table analysis is under `analysis/multi_warmup/`, cross-stage execution
is under `runner/multi_warmup/`, and Slurm orchestration/login-shell setup stays here.

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `run_warm_up_joint_fitting_results_extra_hpc.sh` | Run one or more analysis stages on an HPC node, optionally detached. | Joint multi-warmup result root and CLI options. | Analysis tables/reductions plus an HPC log. |
| `submit_warm_up_joint_curve_array_hpc.sh` | Submit dense-grid tasks, a merge dependency, and the downstream curve/summary job. | Prepared seed manifest, dense-grid submitter, array backend. | Slurm jobs, task tables, logs, merged tables, and downstream analysis outputs. |

Both scripts accept `--help`. Run them from any working directory; canonical
project-relative defaults are resolved from the repository root.
