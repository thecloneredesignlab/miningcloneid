#!/usr/bin/env bash

set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=docker_runtime.sh
source "${LOCAL_ROOT}/docker_runtime.sh"

if [[ "$#" -eq 0 ]]; then
  cat <<'EOF' >&2
Usage:
  bash run_postfit_docker.sh --fit_dir=/absolute/path/to/seed \
    --scope=invivo|invitro|joint [run_postfit_pipeline.R options]
EOF
  exit 2
fi

o2sd_docker_exec Rscript \
  "${O2SD_WORKFLOW_ROOT}/runner/run_postfit_pipeline.R" "$@"
