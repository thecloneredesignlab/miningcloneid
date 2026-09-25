# Multi-warmup runners

## File registry

- `run_multi_warmup_seed_plan.R`: seed-plan analysis followed by seed-plan figures.
- `run_multi_warmup_landscape_pipeline.R`: landscape table analysis followed by
  cluster figures; task-only stages remain table-only.
- `run_multi_warmup_results_pipeline.R`: result collection, optional integrated
  extra-results execution, and collected-summary figures.
- `render_multi_warmup_results.R`: report figure rendering followed by HTML assembly.
- `run_warm_up_joint_results_pipeline.R`: prepare, embedding, dense-grid curve,
  regression, summary, and visualization orchestration.

HPC launchers and workers remain under `hpc/`; runners do not contain Slurm directives.
