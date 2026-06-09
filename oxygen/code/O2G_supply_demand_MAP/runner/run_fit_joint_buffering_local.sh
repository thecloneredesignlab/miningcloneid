#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

RUNNER_SCRIPT="${SCRIPT_DIR}/run_fit_joint_model_O2G_supply_demand_MAP.sh"
CONFIG_PATH="${PROJECT_ROOT}/oxygen/config/O2G_supply_demand.yaml"
OUT_ROOT="${PROJECT_ROOT}/oxygen/results"
RUN_PREFIX="fit_joint_O2G_buffering_local"
APPEND_RUN_PREFIX_TIMESTAMP="TRUE"
SEEDS_CSV="1"
N_CORES="9"
AUTO_VIZ="TRUE"

usage() {
  cat <<'EOF'
Usage:
  run_fit_joint_buffering_local.sh [options]

Options:
  --config=PATH                       Config YAML path.
  --out_root=PATH                     Output root directory.
  --run_prefix=NAME                   Run directory prefix.
  --append_run_prefix_timestamp=BOOL  Append timestamp to run prefix. Default: TRUE.
  --seeds_csv=CSV                     Comma-separated seeds. Default: 1.
  --n_cores=N                         Number of DEoptim workers. Default: 9.
  --auto_viz=BOOL                     Generate per-seed viz/report. Default: TRUE.
  -h, --help                          Show this help.

Example:
  run_fit_joint_buffering_local.sh --seeds_csv=1,2,3 --n_cores=9 --run_prefix=fit_joint_test
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config=*) CONFIG_PATH="${1#*=}" ;;
    --config) shift; CONFIG_PATH="${1:?Missing value for --config}" ;;
    --out_root=*) OUT_ROOT="${1#*=}" ;;
    --out_root) shift; OUT_ROOT="${1:?Missing value for --out_root}" ;;
    --run_prefix=*) RUN_PREFIX="${1#*=}" ;;
    --run_prefix) shift; RUN_PREFIX="${1:?Missing value for --run_prefix}" ;;
    --append_run_prefix_timestamp=*) APPEND_RUN_PREFIX_TIMESTAMP="${1#*=}" ;;
    --append_run_prefix_timestamp) shift; APPEND_RUN_PREFIX_TIMESTAMP="${1:?Missing value for --append_run_prefix_timestamp}" ;;
    --seeds_csv=*) SEEDS_CSV="${1#*=}" ;;
    --seeds_csv) shift; SEEDS_CSV="${1:?Missing value for --seeds_csv}" ;;
    --n_cores=*) N_CORES="${1#*=}" ;;
    --n_cores) shift; N_CORES="${1:?Missing value for --n_cores}" ;;
    --auto_viz=*) AUTO_VIZ="${1#*=}" ;;
    --auto_viz) shift; AUTO_VIZ="${1:?Missing value for --auto_viz}" ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

if [[ ! -f "${RUNNER_SCRIPT}" ]]; then
  echo "Missing runner script: ${RUNNER_SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${CONFIG_PATH}" ]]; then
  echo "Missing config file: ${CONFIG_PATH}" >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"

echo "Running local O2G joint fit"
echo "  runner: ${RUNNER_SCRIPT}"
echo "  config: ${CONFIG_PATH}"
echo "  out_root: ${OUT_ROOT}"
echo "  run_prefix: ${RUN_PREFIX}"
echo "  append_timestamp: ${APPEND_RUN_PREFIX_TIMESTAMP}"
echo "  seeds_csv: ${SEEDS_CSV}"
echo "  n_cores: ${N_CORES}"
echo "  auto_viz: ${AUTO_VIZ}"

exec bash "${RUNNER_SCRIPT}" \
  --mode=run \
  --config="${CONFIG_PATH}" \
  --out_root="${OUT_ROOT}" \
  --run_prefix="${RUN_PREFIX}" \
  --append_run_prefix_timestamp="${APPEND_RUN_PREFIX_TIMESTAMP}" \
  --seeds_csv="${SEEDS_CSV}" \
  --n_cores="${N_CORES}" \
  --auto_viz="${AUTO_VIZ}"
