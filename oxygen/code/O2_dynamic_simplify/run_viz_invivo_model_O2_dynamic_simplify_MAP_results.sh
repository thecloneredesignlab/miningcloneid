#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIZ_SCRIPT="${SCRIPT_DIR}/viz_invivo_model_O2_dynamic_simplify_MAP_results.R"

if [[ ! -f "${VIZ_SCRIPT}" ]]; then
  echo "ERROR: viz script not found: ${VIZ_SCRIPT}" >&2
  exit 1
fi

# Optional defaults (can be overridden by CLI args below).
FIT_DIR_DEFAULT=""
REPORT_DT_DEFAULT="1"
TOP_N_DEFAULT="6"

if [[ -n "${FIT_DIR_DEFAULT}" ]]; then
  Rscript "${VIZ_SCRIPT}" \
    --fit_dir="${FIT_DIR_DEFAULT}" \
    --report_dt="${REPORT_DT_DEFAULT}" \
    --top_n="${TOP_N_DEFAULT}" \
    "$@"
else
  Rscript "${VIZ_SCRIPT}" \
    --report_dt="${REPORT_DT_DEFAULT}" \
    --top_n="${TOP_N_DEFAULT}" \
    "$@"
fi

