#!/bin/bash
# Submit the joint soft-coupling sigma=0.35 200-seed run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../../../.." && pwd)}"
SUBMIT_SCRIPT="${SUBMIT_SCRIPT:-${SCRIPT_DIR}/submit_o2g_fit.sh}"

CONFIG_PATH="${CONFIG_PATH:-${PROJECT_ROOT}/oxygen/config/O2G_supply_demand_sigma0p35.yaml}"
OUT_ROOT="${OUT_ROOT:-${PROJECT_ROOT}/oxygen/results}"
RUN_PREFIX="${RUN_PREFIX:-fit_joint_sigma0p35_200seed}"

TOTAL_SEEDS="${TOTAL_SEEDS:-200}"
ARRAY_TASKS="${ARRAY_TASKS:-${TOTAL_SEEDS}}"
SEEDS_PER_TASK="${SEEDS_PER_TASK:-1}"
N_CORES="${N_CORES:-22}"
MEM="${MEM:-32G}"
AUTO_VIZ="${AUTO_VIZ:-TRUE}"
R_MODULE="${R_MODULE:-R/4.4}"
DRY_RUN="${DRY_RUN:-FALSE}"

exec bash "${SUBMIT_SCRIPT}" \
  --fitting_mode=joint \
  --joint_fitting_mode=DIRECT \
  --project_root="${PROJECT_ROOT}" \
  --config_path="${CONFIG_PATH}" \
  --out_root="${OUT_ROOT}" \
  --r_module="${R_MODULE}" \
  --dry_run="${DRY_RUN}" \
  --joint_run_prefix="${RUN_PREFIX}" \
  --joint_total_seeds="${TOTAL_SEEDS}" \
  --joint_array_tasks="${ARRAY_TASKS}" \
  --joint_seeds_per_task="${SEEDS_PER_TASK}" \
  --joint_n_cores="${N_CORES}" \
  --joint_mem="${MEM}" \
  --joint_qos=xxlarge \
  --joint_time_limit=12:00:00 \
  --auto_viz="${AUTO_VIZ}" \
  "$@"
