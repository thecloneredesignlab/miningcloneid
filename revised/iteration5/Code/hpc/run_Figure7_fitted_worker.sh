#!/usr/bin/env bash
set -euo pipefail
workspace="$1";campaign="$2";phase="$3"
source "${workspace}/Code/config/figure78_hpc_paths.env"
# Each worker stages the unchanged private model into its own task runtime;
# never attach the older shared-cache mount policy to this new condition.
unset FIGURE78_DATA_ROOT FIGURE78_FIGURE
task="${SLURM_ARRAY_TASK_ID:?Manifest task index required}"
exec "$FIGURE78_RSCRIPT" "${workspace}/Code/Figures/run_Figure7_fitted_p_misseg.R" \
  --stage=simulate "--phase=$phase" "--workspace=$workspace" \
  "--campaign-root=$campaign" "--task-index=$task"
