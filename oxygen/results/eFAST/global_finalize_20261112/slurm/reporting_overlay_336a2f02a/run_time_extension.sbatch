#!/bin/bash
# Scheduler resources are set explicitly by monitor.py.
set -euo pipefail
runtime_root="$1"
target_job="$2"
export PATH="/cm/shared/apps/slurm/current/bin:/cm/local/apps/slurm/current/bin:${PATH}"
exec python3 -u "${runtime_root}/monitor.py" watch "${target_job}"
