#!/usr/bin/env bash

set -euo pipefail

LOCAL_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=docker_runtime.sh
source "${LOCAL_ROOT}/docker_runtime.sh"

usage() {
  cat <<'EOF'
Usage:
  bash run_full_fit_docker.sh [run_o2_fit.sh options]

Runs the complete fitting chain in Docker: in vivo, then in vitro, then the
only supported joint primary-cluster workflow. The generated separate-fit
directories are passed automatically into joint fitting.

Production-like defaults:
  --fitting_mode=all
  --invivo_total_seeds=500
  --invitro_total_seeds=500
  --joint_total_seeds=500

Set O2SD_*_TOTAL_SEEDS or pass explicit runner options to override. Use
--fitting_mode=joint with both source result directories to run only the joint
portion. Use --dry_run=TRUE to inspect the workflow without fitting.
EOF
}

for arg in "$@"; do
  if [[ "${arg}" == "--help" || "${arg}" == "-h" ]]; then
    usage
    exit 0
  fi
done

has_arg() {
  local key="$1"
  shift
  local arg
  for arg in "$@"; do
    [[ "${arg}" == "--${key}="* ]] && return 0
  done
  return 1
}

runner_args=("$@")

has_arg fitting_mode "$@" \
  || runner_args+=("--fitting_mode=all")
has_arg invivo_total_seeds "$@" \
  || runner_args+=("--invivo_total_seeds=${O2SD_INVIVO_TOTAL_SEEDS:-500}")
has_arg invitro_total_seeds "$@" \
  || runner_args+=("--invitro_total_seeds=${O2SD_INVITRO_TOTAL_SEEDS:-500}")
has_arg joint_total_seeds "$@" \
  || runner_args+=("--joint_total_seeds=${O2SD_JOINT_TOTAL_SEEDS:-500}")
has_arg config_path "$@" || has_arg config "$@" \
  || runner_args+=("--config_path=${O2SD_REPO_ROOT}/oxygen/config/O2_supply_demand.yaml")
has_arg out_root "$@" \
  || runner_args+=("--out_root=${O2SD_REPO_ROOT}/oxygen/results")
has_arg append_run_prefix_timestamp "$@" \
  || runner_args+=("--append_run_prefix_timestamp=FALSE")

exec bash "${LOCAL_ROOT}/run_o2_fit_docker.sh" "${runner_args[@]}"
