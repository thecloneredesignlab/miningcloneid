# FixO2 eigen-attractor HPC entrypoints

## File registry

- `fixo2_eigen_attractor_hpc_task.R`: Slurm-array R worker for task-table build,
  per-row eigen-attractor evaluation, and task-row merging. It imports the
  canonical simulation builder and honors `ARRAY_TASK_OFFSET` through the submitter.
- `submit_fixo2_eigen_attractor_array.sh`: submits build/array/merge dependencies
  while preserving the cluster array-size contract.
