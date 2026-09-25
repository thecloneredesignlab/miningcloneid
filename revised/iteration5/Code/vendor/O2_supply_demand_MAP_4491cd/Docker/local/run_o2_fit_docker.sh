#!/usr/bin/env bash

set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=docker_runtime.sh
source "${LOCAL_ROOT}/docker_runtime.sh"

RUNNER="${O2SD_WORKFLOW_ROOT}/runner/run_o2_fit.sh"

if [[ ! -f "${RUNNER}" ]]; then
  echo "Missing canonical fit runner: ${RUNNER}" >&2
  exit 2
fi

has_project_root="FALSE"
for arg in "$@"; do
  case "${arg}" in
    --project_root=*) has_project_root="TRUE" ;;
  esac
done

runner_args=("$@")
if [[ "${has_project_root}" != "TRUE" ]]; then
  runner_args+=("--project_root=${O2SD_REPO_ROOT}")
fi

o2sd_docker_exec bash "${RUNNER}" "${runner_args[@]}"
