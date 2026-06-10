#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNIFIED_RUNNER="${SCRIPT_DIR}/run_fit_model_O2_supply_demand_MAP.sh"

has_mode=false
for arg in "$@"; do
  case "${arg}" in
    --mode|--mode=*)
      has_mode=true
      break
      ;;
  esac
done

if "${has_mode}"; then
  exec bash "${UNIFIED_RUNNER}" --fit_joint "$@"
else
  exec bash "${UNIFIED_RUNNER}" --fit_joint --mode=run "$@"
fi
