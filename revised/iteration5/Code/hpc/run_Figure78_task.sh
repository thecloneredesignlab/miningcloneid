#!/usr/bin/env bash
# Execute one ALREADY PREPARED manifest task. This script never submits jobs.
# A future approved submission uses xxlarge and the full manifest index range,
# without a %N concurrency limit. Activate the approved R/container environment.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
: "${FIGURE78_DATA_ROOT:?Set the absolute prepared data/Figures root}"
: "${FIGURE78_FIGURE:?Set 7 or 8}"
: "${SLURM_ARRAY_TASK_ID:?This worker requires the real manifest task index}"
case "$FIGURE78_FIGURE" in 7|8) ;; *) echo "FIGURE78_FIGURE must be 7 or 8" >&2; exit 2 ;; esac
[[ "$FIGURE78_DATA_ROOT" == /* ]] || { echo "FIGURE78_DATA_ROOT must be absolute" >&2; exit 2; }
python3 "${script_dir}/../Figures/figure78_entry_tools.py" check-prepared \
  --workspace="$FIGURE78_WORKSPACE" --figure="$FIGURE78_FIGURE" --data-root="$FIGURE78_DATA_ROOT" \
  --task-index="$SLURM_ARRAY_TASK_ID"
exec "$FIGURE78_RSCRIPT" \
  "$FIGURE78_WORKSPACE/Code/Figures/run_Figure${FIGURE78_FIGURE}.R" \
  --stage=simulate --profile=full \
  "--workspace=$FIGURE78_WORKSPACE" "--data-root=$FIGURE78_DATA_ROOT" \
  "--task-index=$SLURM_ARRAY_TASK_ID"
